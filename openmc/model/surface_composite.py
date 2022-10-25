from abc import ABC, abstractmethod
from copy import copy
from math import sqrt, pi, sin, cos
import numpy as np

import openmc
from openmc.checkvalue import check_greater_than, check_value


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
