from abc import ABC, abstractmethod
from copy import copy
from math import sqrt, pi, sin, cos, isclose
import warnings
import operator

import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from matplotlib.path import Path

import openmc
from openmc.checkvalue import (check_greater_than, check_value,
                               check_iterable_type, check_length)


class CompositeSurface(ABC):
    """Multiple primitive surfaces combined into a composite surface"""

    def translate(self, vector, inplace=False):
        surf = self if inplace else copy(self)
        for name in self._surface_names:
            s = getattr(surf, name)
            setattr(surf, name, s.translate(vector, inplace))
        return surf

    def rotate(self, rotation, pivot=(0., 0., 0.), order='xyz', inplace=False):
        surf = copy(self)
        for name in self._surface_names:
            s = getattr(surf, name)
            setattr(surf, name, s.rotate(rotation, pivot, order, inplace))
        return surf

    @property
    def boundary_type(self):
        return getattr(self, self._surface_names[0]).boundary_type

    @boundary_type.setter
    def boundary_type(self, boundary_type):
        # Set boundary type on underlying surfaces, but not for ambiguity plane
        # on one-sided cones
        for name in self._surface_names:
            if name != 'plane':
                getattr(self, name).boundary_type = boundary_type

    def __repr__(self):
        return "<{} at 0x{:x}>".format(type(self).__name__, id(self))

    @property
    @abstractmethod
    def _surface_names(self):
        """Iterable of attribute names corresponding to underlying surfaces."""

    @abstractmethod
    def __pos__(self):
        """Return the positive half-space of the composite surface."""

    @abstractmethod
    def __neg__(self):
        """Return the negative half-space of the composite surface."""


class CylinderSector(CompositeSurface):
    """Infinite cylindrical sector composite surface.

    A cylinder sector is composed of two cylindrical and two planar surfaces.
    The cylindrical surfaces are concentric, and the planar surfaces intersect
    the central axis of the cylindrical surfaces.

    This class acts as a proper surface, meaning that unary `+` and `-`
    operators applied to it will produce a half-space. The negative
    side is defined to be the region inside of the cylinder sector.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    r1 : float
        Inner radius of sector. Must be less than r2.
    r2 : float
        Outer radius of sector. Must be greater than r1.
    theta1 : float
        Clockwise-most bound of sector in degrees. Assumed to be in the
        counterclockwise direction with respect to the first basis axis
        (+y, +z, or +x). Must be less than :attr:`theta2`.
    theta2 : float
        Counterclockwise-most bound of sector in degrees. Assumed to be in the
        counterclockwise direction with respect to the first basis axis
        (+y, +z, or +x). Must be greater than :attr:`theta1`.
    center : iterable of float
       Coordinate for central axes of cylinders in the (y, z), (z, x), or (x, y)
       basis. Defaults to (0,0).
    axis : {'x', 'y', 'z'}
        Central axis of the cylinders defining the inner and outer surfaces of
        the sector. Defaults to 'z'.
    **kwargs : dict
        Keyword arguments passed to the :class:`Cylinder` and
        :class:`Plane` constructors.

    Attributes
    ----------
    outer_cyl : openmc.ZCylinder, openmc.YCylinder, or openmc.XCylinder
        Outer cylinder surface.
    inner_cyl : openmc.ZCylinder, openmc.YCylinder, or openmc.XCylinder
        Inner cylinder surface.
    plane1 : openmc.Plane
        Plane at angle :math:`\\theta_1` relative to the first basis axis.
    plane2 : openmc.Plane
        Plane at angle :math:`\\theta_2` relative to the first basis axis.

    """

    _surface_names = ('outer_cyl', 'inner_cyl', 'plane1', 'plane2')

    def __init__(self,
                 r1,
                 r2,
                 theta1,
                 theta2,
                 center=(0.,0.),
                 axis='z',
                 **kwargs):

        if r2 <= r1:
            raise ValueError(f'r2 must be greater than r1.')

        if theta2 <= theta1:
            raise ValueError(f'theta2 must be greater than theta1.')

        phi1 = pi / 180 * theta1
        phi2 = pi / 180 * theta2

        # Coords for axis-perpendicular planes
        p1 = np.array([0., 0., 1.])

        p2_plane1 = np.array([r1 * cos(phi1), r1 * sin(phi1), 0.])
        p3_plane1 = np.array([r2 * cos(phi1), r2 * sin(phi1), 0.])

        p2_plane2 = np.array([r1 * cos(phi2), r1 * sin(phi2), 0.])
        p3_plane2 = np.array([r2 * cos(phi2), r2 * sin(phi2), 0.])

        points = [p1, p2_plane1, p3_plane1, p2_plane2, p3_plane2]
        if axis == 'z':
            coord_map = [0, 1, 2]
            self.inner_cyl = openmc.ZCylinder(*center, r1, **kwargs)
            self.outer_cyl = openmc.ZCylinder(*center, r2, **kwargs)
        elif axis == 'y':
            coord_map = [1, 2, 0]
            self.inner_cyl = openmc.YCylinder(*center, r1, **kwargs)
            self.outer_cyl = openmc.YCylinder(*center, r2, **kwargs)
        elif axis == 'x':
            coord_map = [2, 0, 1]
            self.inner_cyl = openmc.XCylinder(*center, r1, **kwargs)
            self.outer_cyl = openmc.XCylinder(*center, r2, **kwargs)

        for p in points:
            p[:] = p[coord_map]

        self.plane1 = openmc.Plane.from_points(p1, p2_plane1, p3_plane1,
                                               **kwargs)
        self.plane2 = openmc.Plane.from_points(p1, p2_plane2, p3_plane2,
                                               **kwargs)

    @classmethod
    def from_theta_alpha(cls,
                         r1,
                         r2,
                         theta,
                         alpha,
                         center = (0.,0.),
                         axis='z',
                         **kwargs):
        r"""Alternate constructor for :class:`CylinderSector`. Returns a
        :class:`CylinderSector` object based on a central angle :math:`\theta`
        and an angular offset :math:`\alpha`. Note that
        :math:`\theta_1 = \alpha` and :math:`\theta_2 = \alpha + \theta`.

        Parameters
        ----------
        r1 : float
            Inner radius of sector. Must be less than r2.
        r2 : float
            Outer radius of sector. Must be greater than r1.
        theta : float
            Central angle, :math:`\theta`, of the sector in degrees. Must be
            greater that 0 and less than 360.
        alpha : float
            Angular offset, :math:`\alpha`, of sector in degrees.
            The offset is in the counter-clockwise direction
            with respect to the first basis axis (+y, +z, or +x). Note that
            negative values translate to an offset in the clockwise direction.
        center : iterable of float
            Coordinate for central axes of cylinders in the (y, z), (z, x), or (x, y)
            basis. Defaults to (0,0).
        axis : {'x', 'y', 'z'}
            Central axis of the cylinders defining the inner and outer surfaces
            of the sector. Defaults to 'z'.
        **kwargs : dict
            Keyword arguments passed to the :class:`Cylinder` and
            :class:`Plane` constructors.

        Returns
        -------
        CylinderSector
            CylinderSector with the given central angle at the given
            offset.
        """
        if theta >= 360. or theta <= 0:
            raise ValueError(f'theta must be less than 360 and greater than 0.')

        theta1 = alpha
        theta2 = alpha + theta

        return cls(r1, r2, theta1, theta2, center=center, axis=axis, **kwargs)

    def __neg__(self):
        return -self.outer_cyl & +self.inner_cyl & -self.plane1 & +self.plane2

    def __pos__(self):
        return +self.outer_cyl | -self.inner_cyl | +self.plane1 | -self.plane2


class IsogonalOctagon(CompositeSurface):
    r"""Infinite isogonal octagon composite surface

    An isogonal octagon is composed of eight planar surfaces. The prism is
    parallel to the x, y, or z axis. The remaining two axes (y and z, z and x,
    or x and y) serve as a basis for constructing the surfaces. Two surfaces
    are parallel to the first basis axis, two surfaces are parallel
    to the second basis axis, and the remaining four surfaces intersect both
    basis axes at 45 degree angles.

    This class acts as a proper surface, meaning that unary `+` and `-`
    operators applied to it will produce a half-space. The negative side is
    defined to be the region inside of the octogonal prism.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    center : iterable of float
        Coordinate for the central axis of the octagon in the
        (y, z), (z, x), or (x, y) basis.
    r1 : float
        Half-width of octagon across its basis axis-parallel sides in units
        of cm. Must be less than :math:`r_2\sqrt{2}`.
    r2 : float
        Half-width of octagon across its basis axis intersecting sides in
        units of cm. Must be less than than :math:`r_1\sqrt{2}`.
    axis : {'x', 'y', 'z'}
        Central axis of octagon. Defaults to 'z'
    **kwargs
        Keyword arguments passed to underlying plane classes

    Attributes
    ----------
    top : openmc.ZPlane, openmc.XPlane, or openmc.YPlane
        Top planar surface of octagon
    bottom : openmc.ZPlane, openmc.XPlane, or openmc.YPlane
        Bottom planar surface of octagon
    right : openmc.YPlane, openmc.ZPlane, or openmc.XPlane
        Right planar surface of octagon
    left : openmc.YPlane, openmc.ZPlane, or openmc.XPlane
        Left planar surface of octagon
    upper_right : openmc.Plane
        Upper right planar surface of octagon
    lower_right : openmc.Plane
        Lower right planar surface of octagon
    lower_left : openmc.Plane
        Lower left planar surface of octagon
    upper_left : openmc.Plane
        Upper left planar surface of octagon

    """

    _surface_names = ('top', 'bottom',
                      'upper_right', 'lower_left',
                      'right', 'left',
                      'lower_right', 'upper_left')

    def __init__(self, center, r1, r2, axis='z', **kwargs):
        c1, c2 = center

        # Coords for axis-perpendicular planes
        ctop = c1 + r1
        cbottom = c1 - r1

        cright = c2 + r1
        cleft = c2 - r1

        # Side lengths
        if r2 > r1 * sqrt(2):
            raise ValueError(f'r2 is greater than sqrt(2) * r1. Octagon' +
                             ' may be erroneous.')
        if r1 > r2 * sqrt(2):
            raise ValueError(f'r1 is greater than sqrt(2) * r2. Octagon' +
                             ' may be erroneous.')

        L_basis_ax = (r2 * sqrt(2) - r1)

        # Coords for quadrant planes
        p1_ur = np.array([L_basis_ax, r1, 0.])
        p2_ur = np.array([r1, L_basis_ax, 0.])
        p3_ur = np.array([r1, L_basis_ax, 1.])

        p1_lr = np.array([r1, -L_basis_ax, 0.])
        p2_lr = np.array([L_basis_ax, -r1, 0.])
        p3_lr = np.array([L_basis_ax, -r1, 1.])

        points = [p1_ur, p2_ur, p3_ur, p1_lr, p2_lr, p3_lr]

        # Orientation specific variables
        if axis == 'z':
            coord_map = [0, 1, 2]
            self.top = openmc.YPlane(ctop, **kwargs)
            self.bottom = openmc.YPlane(cbottom, **kwargs)
            self.right = openmc.XPlane(cright, **kwargs)
            self.left = openmc.XPlane(cleft, **kwargs)
        elif axis == 'y':
            coord_map = [1, 2, 0]
            self.top = openmc.XPlane(ctop, **kwargs)
            self.bottom = openmc.XPlane(cbottom, **kwargs)
            self.right = openmc.ZPlane(cright, **kwargs)
            self.left = openmc.ZPlane(cleft, **kwargs)
        elif axis == 'x':
            coord_map = [2, 0, 1]
            self.top = openmc.ZPlane(ctop, **kwargs)
            self.bottom = openmc.ZPlane(cbottom, **kwargs)
            self.right = openmc.YPlane(cright, **kwargs)
            self.left = openmc.YPlane(cleft, **kwargs)

        # Put our coordinates in (x,y,z) order
        for p in points:
            p[:] = p[coord_map]

        self.upper_right = openmc.Plane.from_points(p1_ur, p2_ur, p3_ur,
                                                    **kwargs)
        self.lower_right = openmc.Plane.from_points(p1_lr, p2_lr, p3_lr,
                                                    **kwargs)
        self.lower_left = openmc.Plane.from_points(-p1_ur, -p2_ur, -p3_ur,
                                                   **kwargs)
        self.upper_left = openmc.Plane.from_points(-p1_lr, -p2_lr, -p3_lr,
                                                   **kwargs)

    def __neg__(self):
        return -self.top & +self.bottom & -self.right & +self.left & \
            +self.upper_right & +self.lower_right & -self.lower_left & \
            -self.upper_left

    def __pos__(self):
        return +self.top | -self.bottom | +self.right | -self.left | \
            -self.upper_right | -self.lower_right | +self.lower_left | \
            +self.upper_left


class RightCircularCylinder(CompositeSurface):
    """Right circular cylinder composite surface

    A right circular cylinder is composed of a cylinder and two planar surface
    perpendicular to the axis of the cylinder. This class acts as a proper
    surface, meaning that unary `+` and `-` operators applied to it will produce
    a half-space. The negative side is defined to be the region inside of the
    right circular cylinder.

    .. versionadded:: 0.12

    Parameters
    ----------
    center_base : iterable of float
        Cartesian coordinate of the center of the base of the cylinder
    height : float
        Height of the cylinder
    radius : float
        Radius of the cylinder
    axis : {'x', 'y', 'z'}
        Axis of the cylinder
    **kwargs
        Keyword arguments passed to underlying cylinder and plane classes

    Attributes
    ----------
    cyl : openmc.Cylinder
        Underlying cylinder surface
    bottom : openmc.Plane
        Bottom planar surface of the cylinder
    top : openmc.Plane
        Top planar surface of the cylinder

    """
    _surface_names = ('cyl', 'bottom', 'top')

    def __init__(self, center_base, height, radius, axis='z', **kwargs):
        cx, cy, cz = center_base
        check_greater_than('cylinder height', height, 0.0)
        check_greater_than('cylinder radius', radius, 0.0)
        check_value('cylinder axis', axis, ('x', 'y', 'z'))
        if axis == 'x':
            self.cyl = openmc.XCylinder(y0=cy, z0=cz, r=radius, **kwargs)
            self.bottom = openmc.XPlane(x0=cx, **kwargs)
            self.top = openmc.XPlane(x0=cx + height, **kwargs)
        elif axis == 'y':
            self.cyl = openmc.YCylinder(x0=cx, z0=cz, r=radius, **kwargs)
            self.bottom = openmc.YPlane(y0=cy, **kwargs)
            self.top = openmc.YPlane(y0=cy + height, **kwargs)
        elif axis == 'z':
            self.cyl = openmc.ZCylinder(x0=cx, y0=cy, r=radius, **kwargs)
            self.bottom = openmc.ZPlane(z0=cz, **kwargs)
            self.top = openmc.ZPlane(z0=cz + height, **kwargs)

    def __neg__(self):
        return -self.cyl & +self.bottom & -self.top

    def __pos__(self):
        return +self.cyl | -self.bottom | +self.top


class RectangularParallelepiped(CompositeSurface):
    """Rectangular parallelpiped composite surface

    A rectangular parallelpiped is composed of six planar surfaces. This class
    acts as a proper surface, meaning that unary `+` and `-` operators applied
    to it will produce a half-space. The negative side is defined to be the
    region inside of the rectangular parallelpiped.

    .. versionadded:: 0.12

    Parameters
    ----------
    xmin, xmax : float
        Minimum and maximum x coordinates of the parallelepiped
    ymin, ymax : float
        Minimum and maximum y coordinates of the parallelepiped
    zmin, zmax : float
        Minimum and maximum z coordinates of the parallelepiped
    **kwargs
        Keyword arguments passed to underlying plane classes

    Attributes
    ----------
    xmin, xmax : openmc.XPlane
        Sides of the parallelepiped
    ymin, ymax : openmc.YPlane
        Sides of the parallelepiped
    zmin, zmax : openmc.ZPlane
        Sides of the parallelepiped

    """
    _surface_names = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax')

    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax, **kwargs):
        if xmin >= xmax:
            raise ValueError('xmin must be less than xmax')
        if ymin >= ymax:
            raise ValueError('ymin must be less than ymax')
        if zmin >= zmax:
            raise ValueError('zmin must be less than zmax')
        self.xmin = openmc.XPlane(x0=xmin, **kwargs)
        self.xmax = openmc.XPlane(x0=xmax, **kwargs)
        self.ymin = openmc.YPlane(y0=ymin, **kwargs)
        self.ymax = openmc.YPlane(y0=ymax, **kwargs)
        self.zmin = openmc.ZPlane(z0=zmin, **kwargs)
        self.zmax = openmc.ZPlane(z0=zmax, **kwargs)

    def __neg__(self):
        return -self.xmax & +self.xmin & -self.ymax & +self.ymin & -self.zmax & +self.zmin

    def __pos__(self):
        return +self.xmax | -self.xmin | +self.ymax | -self.ymin | +self.zmax | -self.zmin


class XConeOneSided(CompositeSurface):
    """One-sided cone parallel the x-axis

    A one-sided cone is composed of a normal cone surface and an "ambiguity"
    surface that eliminates the ambiguity as to which region of space is
    included. This class acts as a proper surface, meaning that unary `+` and
    `-` operators applied to it will produce a half-space. The negative side is
    defined to be the region inside of the cone.

    .. versionadded:: 0.12

    Parameters
    ----------
    x0 : float, optional
        x-coordinate of the apex. Defaults to 0.
    y0 : float, optional
        y-coordinate of the apex. Defaults to 0.
    z0 : float, optional
        z-coordinate of the apex. Defaults to 0.
    r2 : float, optional
        Parameter related to the aperature. Defaults to 1.
    up : bool
        Whether to select the side of the cone that extends to infinity in the
        positive direction of the coordinate axis (the positive half-space of
        the ambiguity plane)
    **kwargs
        Keyword arguments passed to underlying plane classes

    Attributes
    ----------
    cone : openmc.XCone
        Regular two-sided cone
    plane : openmc.XPlane
        Ambiguity surface
    up : bool
        Whether to select the side of the cone that extends to infinity in the
        positive direction of the coordinate axis (the positive half-space of
        the ambiguity plane)

    """
    _surface_names = ('cone', 'plane')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., up=True, **kwargs):
        check_greater_than('cone R^2', r2, 0.0)
        self.cone = openmc.XCone(x0, y0, z0, r2, **kwargs)
        self.plane = openmc.XPlane(x0)
        self.up = up

    def __neg__(self):
        return -self.cone & (+self.plane if self.up else -self.plane)

    def __pos__(self):
        if self.up:
            return (+self.cone & +self.plane) | -self.plane
        else:
            return (+self.cone & -self.plane) | +self.plane


class YConeOneSided(CompositeSurface):
    """One-sided cone parallel the y-axis

    A one-sided cone is composed of a normal cone surface and an "ambiguity"
    surface that eliminates the ambiguity as to which region of space is
    included. This class acts as a proper surface, meaning that unary `+` and
    `-` operators applied to it will produce a half-space. The negative side is
    defined to be the region inside of the cone.

    .. versionadded:: 0.12

    Parameters
    ----------
    x0 : float, optional
        x-coordinate of the apex. Defaults to 0.
    y0 : float, optional
        y-coordinate of the apex. Defaults to 0.
    z0 : float, optional
        z-coordinate of the apex. Defaults to 0.
    r2 : float, optional
        Parameter related to the aperature. Defaults to 1.
    up : bool
        Whether to select the side of the cone that extends to infinity in the
        positive direction of the coordinate axis (the positive half-space of
        the ambiguity plane)
    **kwargs
        Keyword arguments passed to underlying plane classes

    Attributes
    ----------
    cone : openmc.YCone
        Regular two-sided cone
    plane : openmc.YPlane
        Ambiguity surface
    up : bool
        Whether to select the side of the cone that extends to infinity in the
        positive direction of the coordinate axis (the positive half-space of
        the ambiguity plane)

    """
    _surface_names = ('cone', 'plane')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., up=True, **kwargs):
        check_greater_than('cone R^2', r2, 0.0)
        self.cone = openmc.YCone(x0, y0, z0, r2, **kwargs)
        self.plane = openmc.YPlane(y0)
        self.up = up

    __neg__ = XConeOneSided.__neg__
    __pos__ = XConeOneSided.__pos__


class ZConeOneSided(CompositeSurface):
    """One-sided cone parallel the z-axis

    A one-sided cone is composed of a normal cone surface and an "ambiguity"
    surface that eliminates the ambiguity as to which region of space is
    included. This class acts as a proper surface, meaning that unary `+` and
    `-` operators applied to it will produce a half-space. The negative side is
    defined to be the region inside of the cone.

    .. versionadded:: 0.12

    Parameters
    ----------
    x0 : float, optional
        x-coordinate of the apex. Defaults to 0.
    y0 : float, optional
        y-coordinate of the apex. Defaults to 0.
    z0 : float, optional
        z-coordinate of the apex. Defaults to 0.
    r2 : float, optional
        Parameter related to the aperature. Defaults to 1.
    up : bool
        Whether to select the side of the cone that extends to infinity in the
        positive direction of the coordinate axis (the positive half-space of
        the ambiguity plane)
    **kwargs
        Keyword arguments passed to underlying plane classes

    Attributes
    ----------
    cone : openmc.ZCone
        Regular two-sided cone
    plane : openmc.ZPlane
        Ambiguity surface
    up : bool
        Whether to select the side of the cone that extends to infinity in the
        positive direction of the coordinate axis (the positive half-space of
        the ambiguity plane)

    """
    _surface_names = ('cone', 'plane')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., up=True, **kwargs):
        check_greater_than('cone R^2', r2, 0.0)
        self.cone = openmc.ZCone(x0, y0, z0, r2, **kwargs)
        self.plane = openmc.ZPlane(z0)
        self.up = up

    __neg__ = XConeOneSided.__neg__
    __pos__ = XConeOneSided.__pos__


class Polygon(CompositeSurface):
    """Create a polygon composite surface from a path of closed points.

    Parameters
    ----------
    points : np.ndarray
        An Nx2 array of points defining the vertices of the polygon.
    basis : {'rz', 'xy', 'yz', 'xz'}, optional
        2D basis set for the polygon. The polygon is two dimensional and has
        infinite extent in the third (unspecified) dimension. For example, the
        'xy' basis produces a polygon with infinite extent in the +/- z
        direction. For the 'rz' basis the phi extent is infinite, thus forming
        an axisymmetric surface.

    Attributes
    ----------
    points : np.ndarray
        An Nx2 array of points defining the vertices of the polygon.
    basis : {'rz', 'xy', 'yz', 'xz'}
        2D basis set for the polygon.
    regions : list of openmc.Region
        A list of :class:`openmc.Region` objects, one for each of the convex polygons
        formed during the decomposition of the input polygon.
    region : openmc.Union
        The union of all the regions comprising the polygon.
    """

    def __init__(self, points, basis='rz'):
        check_value('basis', basis, ('xy', 'yz', 'xz', 'rz'))
        self._basis = basis
        points = np.asarray(points, dtype=float)
        check_iterable_type('points', points, float, min_depth=2, max_depth=2)
        check_length('points', points[0, :], 2, 2)

        # If the last point is the same as the first, remove it and make sure
        # there are still at least 3 points for a valid polygon.
        if np.allclose(points[0, :], points[-1, :]):
            points = points[:-1, :]
        check_length('points', points, 3)

        # Order the points counter-clockwise (necessary for offset method)
        self._points = self._make_ccw(points)

        # Create a triangulation of the points.
        self._tri = Delaunay(self._points, qhull_options='QJ')

        # Decompose the polygon into groups of simplices forming convex subsets
        # and get the sets of (surface, operator) pairs defining the polygon
        self._surfsets = self._decompose_polygon_into_convex_sets()

        # Set surface names as required by CompositeSurface protocol
        surfnames = []
        i = 0
        for surfset in self._surfsets:
            for surf, op, on_boundary in surfset:
                if on_boundary:
                    setattr(self, f'surface_{i}', surf)
                    surfnames.append(f'surface_{i}')
                    i += 1
        self._surfnames = tuple(surfnames)

        # Generate a list of regions whose union represents the polygon.
        regions = []
        for surfs_ops in self._surfsets:
            regions.append([op(surf) for surf, op, _ in surfs_ops])
        self._regions = [openmc.Intersection(regs) for regs in regions]

        # Create the union of all the convex subsets
        self._region = openmc.Union(self._regions)

    def __neg__(self):
        return self._region

    def __pos__(self):
        return ~self._region

    @property
    def _surface_names(self):
        return self._surfnames

    @CompositeSurface.boundary_type.setter
    def boundary_type(self, boundary_type):
        if boundary_type != 'transmission':
            warnings.warn("Setting boundary_type to a value other than "
                          "'transmission' on Polygon composite surfaces can "
                          "result in unintended behavior. Please use the "
                          "regions property of the Polygon to generate "
                          "individual openmc.Cell objects to avoid unwanted "
                          "behavior.")
        for name in self._surface_names:
            getattr(self, name).boundary_type = boundary_type

    @property
    def points(self):
        return self._points

    @property
    def basis(self):
        return self._basis

    @property
    def _normals(self):
        """Generate the outward normal unit vectors for the polygon."""
        # Rotation matrix for 90 degree clockwise rotation (-90 degrees about z
        # axis for an 'xy' basis).
        rotation = np.array([[0., 1.], [-1., 0.]])
        # Get the unit vectors that point from one point in the polygon to the
        # next given that they are ordered counterclockwise and that the final
        # point is connected to the first point
        tangents = np.diff(self._points, axis=0, append=[self._points[0, :]])
        tangents /= np.linalg.norm(tangents, axis=-1, keepdims=True)
        # Rotate the tangent vectors clockwise by 90 degrees, which for a
        # counter-clockwise ordered polygon will produce the outward normal
        # vectors.
        return rotation.dot(tangents.T).T

    @property
    def _equations(self):
        normals = self._normals
        equations = np.empty((normals.shape[0], 3))
        equations[:, :2] = normals
        equations[:, 2] = -np.sum(normals*self.points, axis=-1)
        return equations

    @property
    def regions(self):
        return self._regions

    @property
    def region(self):
        return self._region

    def _make_ccw(self, points):
        """Order a set of points counter-clockwise.

        Parameters
        ----------
        points : np.ndarray (Nx2)
            An Nx2 array of coordinate pairs describing the vertices.

        Returns
        -------
        ordered_points : the input points ordered counter-clockwise
        """
        # Calculates twice the signed area of the polygon using the "Shoelace
        # Formula" https://en.wikipedia.org/wiki/Shoelace_formula
        xpts, ypts = points.T

        # If signed area is positive the curve is oriented counter-clockwise
        if np.sum(ypts*(np.roll(xpts, 1) - np.roll(xpts, -1))) > 0:
            return points

        return points[::-1, :]

    def _group_simplices(self, neighbor_map, group=None):
        """Generate a convex grouping of simplices.

        Parameters
        ----------
        neighbor_map : dict
            A map whose keys are simplex indices for simplices inside the polygon
            and whose values are a list of simplex indices that neighbor this
            simplex and are also inside the polygon.
        group : list
            A list of simplex indices that comprise the current convex group.

        Returns
        -------
        group : list
            The list of simplex indices that comprise the complete convex group.
        """
        # If neighbor_map is empty there's nothing left to do
        if not neighbor_map:
            return group
        # If group is empty, grab the next simplex in the dictionary and recurse
        if group is None:
            sidx = next(iter(neighbor_map))
            return self._group_simplices(neighbor_map, group=[sidx])
        # Otherwise use the last simplex in the group 
        else:
            sidx = group[-1]
            # Remove current simplex from dictionary since it is in a group
            neighbors = neighbor_map.pop(sidx, [])
            # For each neighbor check if it is part of the same convex
            # hull as the rest of the group. If yes, recurse. If no, continue on.
            for n in neighbors:
                if n in group or neighbor_map.get(n, None) is None:
                    continue
                test_group = group + [n]
                test_point_idx = np.unique(self._tri.simplices[test_group, :])
                test_points = self._tri.points[test_point_idx]
                # If test_points are convex keep adding to this group
                if len(test_points) == len(ConvexHull(test_points).vertices):
                    group = self._group_simplices(neighbor_map, group=test_group)
            return group

    def _get_convex_hull_surfs(self, qhull):
        """Generate a list of surfaces given by a set of linear equations

        Parameters
        ----------
        qhull : scipy.spatial.ConvexHull
            A ConvexHull object representing the sub-region of the polygon.

        Returns
        -------
        surfs_ops : list of (surface, operator) tuples

        """
        basis = self.basis
        boundary_eqns = self._equations
        # Collect surface/operator pairs such that the intersection of the
        # regions defined by these pairs is the inside of the polygon.
        surfs_ops = []
        # hull facet equation: dx*x + dy*y + c = 0
        for dx, dy, c in qhull.equations:
            # check if this facet is on the boundary of the polygon
            facet_eq = np.array([dx, dy, c])
            on_boundary = any([np.allclose(facet_eq, eq) for eq in boundary_eqns])
            # Check if the facet is horizontal
            if isclose(dx, 0, abs_tol=1e-8):
                if basis in ('xz', 'yz', 'rz'):
                    surf = openmc.ZPlane(z0=-c/dy)
                else:
                    surf = openmc.YPlane(y0=-c/dy)
                # if (0, 1).(dx, dy) < 0 we want positive halfspace instead
                op = operator.pos if dy < 0 else operator.neg
            # Check if the facet is vertical
            elif isclose(dy, 0, abs_tol=1e-8):
                if basis in ('xy', 'xz'):
                    surf = openmc.XPlane(x0=-c/dx)
                elif basis == 'yz':
                    surf = openmc.YPlane(y0=-c/dx)
                else:
                    surf = openmc.ZCylinder(r=-c/dx)
                # if (1, 0).(dx, dy) < 0 we want positive halfspace instead
                op = operator.pos if dx < 0 else operator.neg
            # Otherwise the facet is at an angle
            else:
                op = operator.neg
                if basis == 'xy':
                    surf = openmc.Plane(a=dx, b=dy, d=-c)
                elif basis == 'yz':
                    surf = openmc.Plane(b=dx, c=dy, d=-c)
                elif basis == 'xz':
                    surf = openmc.Plane(a=dx, c=dy, d=-c)
                else:
                    y0 = -c/dy
                    r2 = dy**2 / dx**2
                    # Check if the *slope* of the facet is positive. If dy/dx < 0
                    # then we want up to be True for the one-sided cones.
                    up = dy / dx < 0
                    surf = openmc.model.ZConeOneSided(z0=y0, r2=r2, up=up)
                    # if (1, -1).(dx, dy) < 0 for up cones we want positive halfspace
                    # if (1, 1).(dx, dy) < 0 for down cones we want positive halfspace
                    # otherwise we keep the negative halfspace operator
                    if (up and dx - dy < 0) or (not up and dx + dy < 0):
                        op = operator.pos

            surfs_ops.append((surf, op, on_boundary))

        return surfs_ops

    def _decompose_polygon_into_convex_sets(self):
        """Decompose the Polygon into a set of convex polygons.

        Returns
        -------
        surfsets : a list of lists of surface, operator pairs
        """

        # Get centroids of all the simplices and determine if they are inside
        # the polygon defined by input vertices or not.
        centroids = np.mean(self._points[self._tri.simplices], axis=1)
        in_polygon = Path(self._points).contains_points(centroids)
        self._in_polygon = in_polygon

        # Build a map with keys of simplex indices inside the polygon whose
        # values are lists of that simplex's neighbors also inside the
        # polygon
        neighbor_map = {}
        for i, nlist in enumerate(self._tri.neighbors):
            if not in_polygon[i]:
                continue
            neighbor_map[i] = [n for n in nlist if in_polygon[n] and n >=0]

        # Get the groups of simplices forming convex polygons whose union
        # comprises the full input polygon. While there are still simplices
        # left in the neighbor map, group them together into convex sets.
        groups = []
        while neighbor_map:
            groups.append(self._group_simplices(neighbor_map))
        self._groups = groups

        # Generate lists of (surface, operator) pairs for each convex
        # sub-region.
        surfsets = []
        for group in groups:
            # Find all the unique points in the convex group of simplices,
            # generate the convex hull and find the resulting surfaces and
            # unary operators that represent this convex subset of the polygon.
            idx = np.unique(self._tri.simplices[group, :])
            qhull = ConvexHull(self._tri.points[idx, :])
            surf_ops = self._get_convex_hull_surfs(qhull)
            surfsets.append(surf_ops)
        return surfsets

    def offset(self, distance):
        """Offset this polygon by a set distance

        Parameters
        ----------
        distance : float
            The distance to offset the polygon by. Positive is outward
            (expanding) and negative is inward (shrinking).

        Returns
        -------
        offset_polygon : openmc.model.Polygon
        """
        normals = np.insert(self._normals, 0, self._normals[-1, :], axis=0)
        cos2theta = np.sum(normals[1:, :]*normals[:-1, :], axis=-1, keepdims=True)
        costheta = np.cos(np.arccos(cos2theta) / 2)
        nvec = (normals[1:, :] + normals[:-1, :])
        unit_nvec = nvec / np.linalg.norm(nvec, axis=-1, keepdims=True)
        disp_vec = distance / costheta * unit_nvec

        return type(self)(self.points + disp_vec, basis=self.basis)
