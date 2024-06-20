from abc import ABC, abstractmethod
from collections.abc import Iterable
from copy import copy
from functools import partial
from math import sqrt, pi, sin, cos, isclose
from numbers import Real
import warnings
import operator
from typing import Sequence

import numpy as np
from scipy.spatial import ConvexHull, Delaunay

import openmc
from openmc.checkvalue import (check_greater_than, check_value, check_less_than,
                               check_iterable_type, check_length, check_type)


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
        return f"<{type(self).__name__} at 0x{id(self):x}>"

    @property
    @abstractmethod
    def _surface_names(self):
        """Iterable of attribute names corresponding to underlying surfaces."""

    @abstractmethod
    def __neg__(self) -> openmc.Region:
        """Return the negative half-space of the composite surface."""

    def __pos__(self) -> openmc.Region:
        """Return the positive half-space of the composite surface."""
        return ~(-self)


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
       Coordinate for central axes of cylinders in the (y, z), (x, z), or (x, y)
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
            raise ValueError('r2 must be greater than r1.')

        if theta2 <= theta1:
            raise ValueError('theta2 must be greater than theta1.')

        phi1 = pi / 180 * theta1
        phi2 = pi / 180 * theta2

        # Coords for axis-perpendicular planes
        p1 = np.array([center[0], center[1], 1.])

        p2_plane1 = np.array([r1 * cos(phi1) + center[0], r1 * sin(phi1) + center[1], 0.])
        p3_plane1 = np.array([r2 * cos(phi1) + center[0], r2 * sin(phi1) + center[1], 0.])

        p2_plane2 = np.array([r1 * cos(phi2) + center[0], r1 * sin(phi2)+ center[1], 0.])
        p3_plane2 = np.array([r2 * cos(phi2) + center[0], r2 * sin(phi2)+ center[1], 0.])

        points = [p1, p2_plane1, p3_plane1, p2_plane2, p3_plane2]
        if axis == 'z':
            coord_map = [0, 1, 2]
            self.inner_cyl = openmc.ZCylinder(*center, r1, **kwargs)
            self.outer_cyl = openmc.ZCylinder(*center, r2, **kwargs)
        elif axis == 'y':
            coord_map = [0, 2, 1]
            self.inner_cyl = openmc.YCylinder(*center, r1, **kwargs)
            self.outer_cyl = openmc.YCylinder(*center, r2, **kwargs)
        elif axis == 'x':
            coord_map = [2, 0, 1]
            self.inner_cyl = openmc.XCylinder(*center, r1, **kwargs)
            self.outer_cyl = openmc.XCylinder(*center, r2, **kwargs)

        # Reorder the points to correspond to the correct central axis
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
            Coordinate for central axes of cylinders in the (y, z), (x, z), or
            (x, y) basis. Defaults to (0,0).
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
            raise ValueError('theta must be less than 360 and greater than 0.')

        theta1 = alpha
        theta2 = alpha + theta

        return cls(r1, r2, theta1, theta2, center=center, axis=axis, **kwargs)

    def __neg__(self):
        if isinstance(self.inner_cyl, openmc.YCylinder):
            return -self.outer_cyl & +self.inner_cyl & +self.plane1 & -self.plane2
        else:
            return -self.outer_cyl & +self.inner_cyl & -self.plane1 & +self.plane2

    def __pos__(self):
        if isinstance(self.inner_cyl, openmc.YCylinder):
            return +self.outer_cyl | -self.inner_cyl | -self.plane1 | +self.plane2
        else:
            return +self.outer_cyl | -self.inner_cyl | +self.plane1 | -self.plane2


class IsogonalOctagon(CompositeSurface):
    r"""Infinite isogonal octagon composite surface

    An isogonal octagon is composed of eight planar surfaces. The prism is
    parallel to the x, y, or z axis. The remaining two axes (y and z, x and z,
    or x and y) serve as a basis for constructing the surfaces. Two surfaces
    are parallel to the first basis axis, two surfaces are parallel
    to the second basis axis, and the remaining four surfaces intersect both
    basis axes at 45 degree angles.

    This class acts as a proper surface, meaning that unary `+` and `-`
    operators applied to it will produce a half-space. The negative side is
    defined to be the region inside of the octagonal prism.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    center : iterable of float
        Coordinate for the central axis of the octagon in the
        (y, z), (x, z), or (x, y) basis depending on the axis parameter.
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

        # Coordinates for axis-perpendicular planes
        cright = c1 + r1
        cleft = c1 - r1

        ctop = c2 + r1
        cbottom = c2 - r1

        # Side lengths
        if r2 > r1 * sqrt(2):
            raise ValueError('r2 is greater than sqrt(2) * r1. Octagon' +
                             ' may be erroneous.')
        if r1 > r2 * sqrt(2):
            raise ValueError('r1 is greater than sqrt(2) * r2. Octagon' +
                             ' may be erroneous.')

        L_basis_ax = (r2 * sqrt(2) - r1)

        # Coordinates for quadrant planes
        p1_ur = np.array([L_basis_ax, r1, 0.])
        p2_ur = np.array([r1, L_basis_ax, 0.])
        p3_ur = np.array([r1, L_basis_ax, 1.])

        p1_lr = np.array([r1, -L_basis_ax, 0.])
        p2_lr = np.array([L_basis_ax, -r1, 0.])
        p3_lr = np.array([L_basis_ax, -r1, 1.])

        p1_ll = -p1_ur
        p2_ll = -p2_ur
        p3_ll = -p3_ur

        p1_ul = -p1_lr
        p2_ul = -p2_lr
        p3_ul = -p3_lr

        points = [p1_ur, p2_ur, p3_ur, p1_lr, p2_lr, p3_lr,
                  p1_ll, p2_ll, p3_ll, p1_ul, p2_ul, p3_ul]

        # Orientation specific variables
        if axis == 'z':
            coord_map = [0, 1, 2]
            self.top = openmc.YPlane(ctop, **kwargs)
            self.bottom = openmc.YPlane(cbottom, **kwargs)
            self.right = openmc.XPlane(cright, **kwargs)
            self.left = openmc.XPlane(cleft, **kwargs)
        elif axis == 'y':
            coord_map = [0, 2, 1]
            self.top = openmc.ZPlane(ctop, **kwargs)
            self.bottom = openmc.ZPlane(cbottom, **kwargs)
            self.right = openmc.XPlane(cright, **kwargs)
            self.left = openmc.XPlane(cleft, **kwargs)
        elif axis == 'x':
            coord_map = [2, 0, 1]
            self.top = openmc.ZPlane(ctop, **kwargs)
            self.bottom = openmc.ZPlane(cbottom, **kwargs)
            self.right = openmc.YPlane(cright, **kwargs)
            self.left = openmc.YPlane(cleft, **kwargs)
        self.axis = axis

        # Put our coordinates in (x,y,z) order and add the offset
        for p in points:
            p[0] += c1
            p[1] += c2
            p[:] = p[coord_map]

        self.upper_right = openmc.Plane.from_points(p1_ur, p2_ur, p3_ur,
                                                    **kwargs)
        self.lower_right = openmc.Plane.from_points(p1_lr, p2_lr, p3_lr,
                                                    **kwargs)
        self.lower_left = openmc.Plane.from_points(p1_ll, p2_ll, p3_ll,
                                                   **kwargs)
        self.upper_left = openmc.Plane.from_points(p1_ul, p2_ul, p3_ul,
                                                   **kwargs)

    def __neg__(self):
        if self.axis == 'y':
            region = -self.top & +self.bottom & -self.right & +self.left & \
                -self.upper_right & -self.lower_right & +self.lower_left & \
                +self.upper_left
        else:
            region = -self.top & +self.bottom & -self.right & +self.left & \
                +self.upper_right & +self.lower_right & -self.lower_left & \
                -self.upper_left

        return region

    def __pos__(self):
        if self.axis == 'y':
            region = +self.top | -self.bottom | +self.right | -self.left | \
                +self.upper_right | +self.lower_right | -self.lower_left | \
                -self.upper_left
        else:
            region = +self.top | -self.bottom | +self.right | -self.left | \
                -self.upper_right | -self.lower_right | +self.lower_left | \
                +self.upper_left
        return region


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
    upper_fillet_radius : float
        Upper edge fillet radius in [cm].
    lower_fillet_radius : float
        Lower edge fillet radius in [cm].
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
    upper_fillet_torus : openmc.Torus
        Surface that creates the filleted edge for the upper end of the
        cylinder. Only present if :attr:`upper_fillet_radius` is set.
    upper_fillet_cylinder : openmc.Cylinder
        Surface that bounds :attr:`upper_fillet_torus` radially. Only present
        if :attr:`upper_fillet_radius` is set.
    upper_fillet_plane : openmc.Plane
        Surface that bounds :attr:`upper_fillet_torus` axially. Only present if
        :attr:`upper_fillet_radius` is set.
    lower_fillet_torus : openmc.Torus
        Surface that creates the filleted edge for the lower end of the
        cylinder. Only present if :attr:`lower_fillet_radius` is set.
    lower_fillet_cylinder : openmc.Cylinder
        Surface that bounds :attr:`lower_fillet_torus` radially. Only present
        if :attr:`lower_fillet_radius` is set.
    lower_fillet_plane : openmc.Plane
        Surface that bounds :attr:`lower_fillet_torus` axially. Only present if
        :attr:`lower_fillet_radius` is set.

    """
    _surface_names = ('cyl', 'bottom', 'top')

    def __init__(self, center_base, height, radius, axis='z',
                 upper_fillet_radius=0., lower_fillet_radius=0., **kwargs):
        cx, cy, cz = center_base
        check_greater_than('cylinder height', height, 0.0)
        check_greater_than('cylinder radius', radius, 0.0)
        check_value('cylinder axis', axis, ('x', 'y', 'z'))
        check_type('upper_fillet_radius', upper_fillet_radius, float)
        check_less_than('upper_fillet_radius', upper_fillet_radius,
                        radius, equality=True)
        check_type('lower_fillet_radius', lower_fillet_radius, float)
        check_less_than('lower_fillet_radius', lower_fillet_radius,
                        radius, equality=True)

        if axis == 'x':
            self.cyl = openmc.XCylinder(y0=cy, z0=cz, r=radius, **kwargs)
            self.bottom = openmc.XPlane(x0=cx, **kwargs)
            self.top = openmc.XPlane(x0=cx + height, **kwargs)
            x1, x2 = 'y', 'z'
            axcoord, axcoord1, axcoord2 = 0, 1, 2
        elif axis == 'y':
            self.cyl = openmc.YCylinder(x0=cx, z0=cz, r=radius, **kwargs)
            self.bottom = openmc.YPlane(y0=cy, **kwargs)
            self.top = openmc.YPlane(y0=cy + height, **kwargs)
            x1, x2 = 'x', 'z'
            axcoord, axcoord1, axcoord2 = 1, 0, 2
        elif axis == 'z':
            self.cyl = openmc.ZCylinder(x0=cx, y0=cy, r=radius, **kwargs)
            self.bottom = openmc.ZPlane(z0=cz, **kwargs)
            self.top = openmc.ZPlane(z0=cz + height, **kwargs)
            x1, x2 = 'x', 'y'
            axcoord, axcoord1, axcoord2 = 2, 0, 1

        def _create_fillet_objects(axis_args, height, center_base, radius, fillet_radius, pos='upper'):
            axis, x1, x2, axcoord, axcoord1, axcoord2 = axis_args
            fillet_ext = height / 2 - fillet_radius
            sign = 1
            if pos == 'lower':
                sign = -1
            coord = center_base[axcoord] + (height / 2) + sign * fillet_ext

            # cylinder
            cyl_name = f'{pos}_min'
            cylinder_args = {
                x1 + '0': center_base[axcoord1],
                x2 + '0': center_base[axcoord2],
                'r': radius - fillet_radius
            }
            cls = getattr(openmc, f'{axis.upper()}Cylinder')
            cyl = cls(name=f'{cyl_name} {axis}', **cylinder_args)

            #torus
            tor_name = f'{axis} {pos}'
            tor_args = {
                'a': radius - fillet_radius,
                'b': fillet_radius,
                'c': fillet_radius,
                x1 + '0': center_base[axcoord1],
                x2 + '0': center_base[axcoord2],
                axis + '0': coord
            }
            cls = getattr(openmc, f'{axis.upper()}Torus')
            torus = cls(name=tor_name, **tor_args)

            # plane
            p_name = f'{pos} ext'
            p_args = {axis + '0': coord}
            cls = getattr(openmc, f'{axis.upper()}Plane')
            plane = cls(name=p_name, **p_args)

            return cyl, torus, plane

        if upper_fillet_radius > 0. or lower_fillet_radius > 0.:
            if 'boundary_type' in kwargs:
                if kwargs['boundary_type'] == 'periodic':
                    raise ValueError('Periodic boundary conditions not permitted when '
                                     'rounded corners are used.')

            axis_args = (axis, x1, x2, axcoord, axcoord1, axcoord2)
            if upper_fillet_radius > 0.:
                cylinder, torus, plane = _create_fillet_objects(
                    axis_args, height, center_base, radius, upper_fillet_radius)
                self.upper_fillet_cylinder = cylinder
                self.upper_fillet_torus = torus
                self.upper_fillet_plane = plane
                self._surface_names += ('upper_fillet_cylinder',
                                        'upper_fillet_torus',
                                        'upper_fillet_plane')

            if lower_fillet_radius > 0.:
                cylinder, torus, plane = _create_fillet_objects(
                    axis_args, height, center_base, radius, lower_fillet_radius,
                    pos='lower'
                )
                self.lower_fillet_cylinder = cylinder
                self.lower_fillet_torus = torus
                self.lower_fillet_plane = plane

                self._surface_names += ('lower_fillet_cylinder',
                                        'lower_fillet_torus',
                                        'lower_fillet_plane')

    def _get_fillet(self):
        upper_fillet = self._get_upper_fillet()
        lower_fillet = self._get_lower_fillet()
        has_upper_fillet = upper_fillet is not None
        has_lower_fillet = lower_fillet is not None
        if has_lower_fillet and has_upper_fillet:
            fillet = lower_fillet | upper_fillet
        elif has_upper_fillet and not has_lower_fillet:
            fillet = upper_fillet
        elif not has_upper_fillet and has_lower_fillet:
            fillet = lower_fillet
        else:
            fillet = None
        return fillet

    def _get_upper_fillet(self):
        has_upper_fillet = hasattr(self, 'upper_fillet_plane')
        if has_upper_fillet:
            upper_fillet = +self.upper_fillet_cylinder & +self.upper_fillet_torus & +self.upper_fillet_plane
        else:
            upper_fillet = None
        return upper_fillet

    def _get_lower_fillet(self):
        has_lower_fillet = hasattr(self, 'lower_fillet_plane')
        if has_lower_fillet:
            lower_fillet = +self.lower_fillet_cylinder & +self.lower_fillet_torus & -self.lower_fillet_plane
        else:
            lower_fillet = None
        return lower_fillet

    def __neg__(self):
        prism = -self.cyl & +self.bottom & -self.top
        fillet = self._get_fillet()
        if fillet is not None:
            prism = prism & ~fillet
        return prism

    def __pos__(self):
        prism = +self.cyl | -self.bottom | +self.top
        fillet = self._get_fillet()
        if fillet is not None:
            prism = prism | fillet
        return prism


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

    A one-sided cone is composed of a normal cone surface and a "disambiguation"
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
        Parameter related to the aperture [:math:`\\rm cm^2`].
        It can be interpreted as the increase in the radius squared per cm along
        the cone's axis of revolution.
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
        Disambiguation surface
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

    A one-sided cone is composed of a normal cone surface and a "disambiguation"
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
        Parameter related to the aperture [:math:`\\rm cm^2`].
        It can be interpreted as the increase in the radius squared per cm along
        the cone's axis of revolution.
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
        Disambiguation surface
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

    A one-sided cone is composed of a normal cone surface and a "disambiguation"
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
        Parameter related to the aperture [:math:`\\rm cm^2`].
        It can be interpreted as the increase in the radius squared per cm along
        the cone's axis of revolution.
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
        Disambiguation surface
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
    """Polygon formed from a path of closed points.

    .. versionadded:: 0.13.3

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

        # Create a constrained triangulation of the validated points.
        # The constrained triangulation is set to the _tri attribute
        self._constrain_triangulation(self._validate_points(points))

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
        return self._tri.points

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
        tangents = np.diff(self.points, axis=0, append=[self.points[0, :]])
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

    def _validate_points(self, points):
        """Ensure the closed path defined by points does not intersect and is
        oriented counter-clockwise.

        Parameters
        ----------
        points : np.ndarray (Nx2)
            An Nx2 array of coordinate pairs describing the vertices.

        Returns
        -------
        ordered_points : the input points ordered counter-clockwise
        """
        points = np.asarray(points, dtype=float)
        check_iterable_type('points', points, float, min_depth=2, max_depth=2)
        check_length('points', points[0, :], 2, 2)

        # If the last point is the same as the first, remove it and make sure
        # there are still at least 3 points for a valid polygon.
        if np.allclose(points[0, :], points[-1, :]):
            points = points[:-1, :]
        check_length('points', points, 3)

        if len(points) != len(np.unique(points, axis=0)):
            raise ValueError('Duplicate points were detected in the Polygon input')

        # Order the points counter-clockwise (necessary for offset method)
        # Calculates twice the signed area of the polygon using the "Shoelace
        # Formula" https://en.wikipedia.org/wiki/Shoelace_formula
        # If signed area is positive the curve is oriented counter-clockwise.
        # If the signed area is negative the curve is oriented clockwise.
        xpts, ypts = points.T
        if np.sum(ypts*(np.roll(xpts, 1) - np.roll(xpts, -1))) < 0:
            points = points[::-1, :]

        # Check if polygon is self-intersecting by comparing edges pairwise
        n = len(points)
        for i in range(n):
            p0 = np.append(points[i, :], 0)
            p1 = np.append(points[(i + 1) % n, :], 0)
            for j in range(i + 1, n):
                p2 = np.append(points[j, :], 0)
                p3 = np.append(points[(j + 1) % n, :], 0)
                # Compute orientation of p0 wrt p2->p3 line segment
                cp0 = np.cross(p3-p0, p2-p0)[-1]
                # Compute orientation of p1 wrt p2->p3 line segment
                cp1 = np.cross(p3-p1, p2-p1)[-1]
                # Compute orientation of p2 wrt p0->p1 line segment
                cp2 = np.cross(p1-p2, p0-p2)[-1]
                # Compute orientation of p3 wrt p0->p1 line segment
                cp3 = np.cross(p1-p3, p0-p3)[-1]

                # Group cross products in an array and find out how many are 0
                cross_products = np.array([[cp0, cp1], [cp2, cp3]])
                cps_near_zero = np.isclose(cross_products, 0).astype(int)
                num_zeros = np.sum(cps_near_zero)

                # Topologies of 2 finite line segments categorized by the number
                # of zero-valued cross products:
                #
                # 0: No 3 points lie on the same line
                # 1: 1 point lies on the same line defined by the other line
                # segment, but is not coincident with either of the points
                # 2: 2 points are coincident, but the line segments are not
                # collinear which guarantees no intersection
                # 3: not possible, except maybe floating point issues?
                # 4: Both line segments are collinear, simply need to check if
                # they overlap or not
                # adapted from algorithm linked below and modified to only
                # consider intersections on the interior of line segments as
                # proper intersections: i.e. segments sharing end points do not
                # count as intersections.
                # https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

                if num_zeros == 0:
                    # If the orientations of p0 and p1 have opposite signs
                    # and the orientations of p2 and p3 have opposite signs
                    # then there is an intersection.
                    if all(np.prod(cross_products, axis=-1) < 0):
                        raise ValueError('Polygon cannot be self-intersecting')
                    continue

                elif num_zeros == 1:
                    # determine which line segment has 2 out of the 3 collinear
                    # points
                    idx = np.argwhere(np.sum(cps_near_zero, axis=-1) == 0)
                    if np.prod(cross_products[idx, :]) < 0:
                        raise ValueError('Polygon cannot be self-intersecting')
                    continue

                elif num_zeros == 2:
                    continue

                elif num_zeros == 3:
                    warnings.warn('Unclear if Polygon is self-intersecting')
                    continue

                else:
                    # All 4 cross products are zero
                    # Determine number of unique points, x span and y span for
                    # both line segments
                    xmin1, xmax1 = min(p0[0], p1[0]), max(p0[0], p1[0])
                    ymin1, ymax1 = min(p0[1], p1[1]), max(p0[1], p1[1])
                    xmin2, xmax2 = min(p2[0], p3[0]), max(p2[0], p3[0])
                    ymin2, ymax2 = min(p2[1], p3[1]), max(p2[1], p3[1])
                    xlap = xmin1 < xmax2 and xmin2 < xmax1
                    ylap = ymin1 < ymax2 and ymin2 < ymax1
                    if xlap or ylap:
                        raise ValueError('Polygon cannot be self-intersecting')
                    continue

        return points

    def _constrain_triangulation(self, points, depth=0):
        """Generate a constrained triangulation by ensuring all edges of the
        Polygon are contained within the simplices.

        Parameters
        ----------
        points : np.ndarray (Nx2)
            An Nx2 array of coordinate pairs describing the vertices. These
            points represent a planar straight line graph.

        Returns
        -------
        None
        """
        # Only attempt the triangulation up to 5 times.
        if depth > 4:
            raise RuntimeError('Could not create a valid triangulation after 5'
                               ' attempts')

        tri = Delaunay(points, qhull_options='QJ')
        # Loop through the boundary edges of the polygon. If an edge is not
        # included in the triangulation, break it into two line segments.
        n = len(points)
        new_pts = []
        for i, j in zip(range(n), range(1, n + 1)):
            # If both vertices of any edge are not found in any simplex, insert
            # a new point between them.
            if not any([i in s and j % n in s for s in tri.simplices]):
                newpt = (points[i, :] + points[j % n, :]) / 2
                new_pts.append((j, newpt))

        # If all the edges are included in the triangulation set it, otherwise
        # try again with additional points inserted on offending edges.
        if not new_pts:
            self._tri = tri
        else:
            for i, pt in new_pts[::-1]:
                points = np.insert(points, i, pt, axis=0)
            self._constrain_triangulation(points, depth=depth + 1)

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
            # Start with smallest neighbor lists
            sidx = sorted(neighbor_map.items(), key=lambda item: len(item[1]))[0][0]
            return self._group_simplices(neighbor_map, group=[sidx])
        # Otherwise use the last simplex in the group
        else:
            sidx = group[-1]
            # Remove current simplex from dictionary since it is in a group
            neighbors = neighbor_map.pop(sidx, [])
            # For each neighbor check if it is part of the same convex
            # hull as the rest of the group. If yes, recurse. If no, continue on.
            for n in neighbors:
                if n in group or neighbor_map.get(n) is None:
                    continue
                test_group = group + [n]
                test_point_idx = np.unique(self._tri.simplices[test_group, :])
                test_points = self.points[test_point_idx]
                test_hull = ConvexHull(test_points, qhull_options='Qc')
                pts_on_hull = len(test_hull.vertices) + len(test_hull.coplanar)
                # If test_points are convex (including coplanar) keep adding to
                # this group
                if len(test_points) == pts_on_hull:
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
        from matplotlib.path import Path

        # Get centroids of all the simplices and determine if they are inside
        # the polygon defined by input vertices or not.
        centroids = np.mean(self.points[self._tri.simplices], axis=1)
        in_polygon = Path(self.points).contains_points(centroids)
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
            qhull = ConvexHull(self.points[idx, :])
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


class CruciformPrism(CompositeSurface):
    """Generalized cruciform prism

    This surface represents a prism parallel to an axis formed by planes at
    multiple distances from the center. Equivalent to the 'gcross' derived
    surface in Serpent.

    .. versionadded:: 0.14.0

    Parameters
    ----------
    distances : iterable of float
        A monotonically increasing (or decreasing) iterable of distances in [cm]
        that form the planes of the generalized cruciform.
    center : iterable of float
        The center of the prism in the two non-parallel axes (e.g., (x, y) when
        axis is 'z') in [cm]
    axis : {'x', 'y', 'z'}
        Axis to which the prism is parallel
    **kwargs
        Keyword arguments passed to underlying plane classes

    """

    def __init__(self, distances, center=(0., 0.), axis='z', **kwargs):
        x0, y0 = center
        self.distances = distances

        if axis == 'x':
            cls_horizontal = openmc.YPlane
            cls_vertical = openmc.ZPlane
        elif axis == 'y':
            cls_horizontal = openmc.XPlane
            cls_vertical = openmc.ZPlane
        elif axis == 'z':
            cls_horizontal = openmc.XPlane
            cls_vertical = openmc.YPlane
        else:
            raise ValueError("axis must be 'x', 'y', or 'z'")

        # Create each planar surface
        surfnames = []
        for i, d in enumerate(distances):
            setattr(self, f'hmin{i}', cls_horizontal(x0 - d, **kwargs))
            setattr(self, f'hmax{i}', cls_horizontal(x0 + d, **kwargs))
            setattr(self, f'vmin{i}', cls_vertical(y0 - d, **kwargs))
            setattr(self, f'vmax{i}', cls_vertical(y0 + d, **kwargs))
            surfnames.extend([f'hmin{i}', f'hmax{i}', f'vmin{i}', f'vmax{i}'])

        # Set _surfnames to satisfy CompositeSurface protocol
        self._surfnames = tuple(surfnames)

    @property
    def _surface_names(self):
        return self._surfnames

    @property
    def distances(self):
        return self._distances

    @distances.setter
    def distances(self, values):
        values = np.array(values, dtype=float)
        # check for positive values
        if not (values > 0).all():
            raise ValueError("distances must be positive")
        # Check for monotonicity
        if (values[1:] > values[:-1]).all() or (values[1:] < values[:-1]).all():
            self._distances = values
        else:
            raise ValueError("distances must be monotonic")

    def __neg__(self):
        n = len(self.distances)
        regions = []
        for i in range(n):
            regions.append(
                +getattr(self, f'hmin{i}') &
                -getattr(self, f'hmax{i}') &
                +getattr(self, f'vmin{n-1-i}') &
                -getattr(self, f'vmax{n-1-i}')
            )
        return openmc.Union(regions)


# Define function to create a plane on given axis
def _plane(axis, name, value, boundary_type='transmission', albedo=1.0):
        cls = getattr(openmc, f'{axis.upper()}Plane')
        return cls(value, name=f'{name} {axis}',
                   boundary_type=boundary_type, albedo=albedo)


class RectangularPrism(CompositeSurface):
    """Infinite rectangular prism bounded by four planar surfaces.

    .. versionadded:: 0.14.0

    Parameters
    ----------
    width : float
        Prism width in units of [cm]. The width is aligned with the x, x, or z
        axes for prisms parallel to the x, y, or z axis, respectively.
    height : float
        Prism height in units of [cm]. The height is aligned with the x, y, or z
        axes for prisms parallel to the x, y, or z axis, respectively.
    axis : {'x', 'y', 'z'}
        Axis with which the infinite length of the prism should be aligned.
    origin : Iterable of two floats
        Origin of the prism. The two floats correspond to (y,z), (x,z) or (x,y)
        for prisms parallel to the x, y or z axis, respectively.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}
        Boundary condition that defines the behavior for particles hitting the
        surfaces comprising the rectangular prism.
    albedo : float, optional
        Albedo of the prism's surfaces as a ratio of particle weight after
        interaction with the surface to the initial weight. Values must be
        positive. Only applicable if the boundary type is 'reflective',
        'periodic', or 'white'.
    corner_radius : float
        Prism corner radius in units of [cm].

    """
    _surface_names = ('min_x1', 'max_x1', 'min_x2', 'max_x2')

    def __init__(
            self,
            width: float,
            height: float,
            axis: str = 'z',
            origin: Sequence[float] = (0., 0.),
            boundary_type: str = 'transmission',
            albedo: float = 1.,
            corner_radius: float = 0.
        ):
        check_type('width', width, Real)
        check_type('height', height, Real)
        check_type('albedo', albedo, Real)
        check_type('corner_radius', corner_radius, Real)
        check_value('axis', axis, ('x', 'y', 'z'))
        check_type('origin', origin, Iterable, Real)

        if axis == 'x':
            x1, x2 = 'y', 'z'
        elif axis == 'y':
            x1, x2 = 'x', 'z'
        else:
            x1, x2 = 'x', 'y'

        # Get cylinder class corresponding to given axis
        cyl = getattr(openmc, f'{axis.upper()}Cylinder')

        # Create container for boundary arguments
        bc_args = {'boundary_type': boundary_type, 'albedo': albedo}

        # Create rectangular region
        self.min_x1 = _plane(x1, 'minimum', -width/2 + origin[0], **bc_args)
        self.max_x1 = _plane(x1, 'maximum', width/2 + origin[0], **bc_args)
        self.min_x2 = _plane(x2, 'minimum', -height/2 + origin[1], **bc_args)
        self.max_x2 = _plane(x2, 'maximum', height/2 + origin[1], **bc_args)
        if boundary_type == 'periodic':
            self.min_x1.periodic_surface = self.max_x1
            self.min_x2.periodic_surface = self.max_x2

        # Handle rounded corners if given
        if corner_radius > 0.:
            if boundary_type == 'periodic':
                raise ValueError('Periodic boundary conditions not permitted when '
                                'rounded corners are used.')

            args = {'r': corner_radius, 'boundary_type': boundary_type, 'albedo': albedo}

            args[x1 + '0'] = origin[0] - width/2 + corner_radius
            args[x2 + '0'] = origin[1] - height/2 + corner_radius
            self.x1_min_x2_min = cyl(name=f'{x1} min {x2} min', **args)

            args[x1 + '0'] = origin[0] - width/2 + corner_radius
            args[x2 + '0'] = origin[1] + height/2 - corner_radius
            self.x1_min_x2_max = cyl(name=f'{x1} min {x2} max', **args)

            args[x1 + '0'] = origin[0] + width/2 - corner_radius
            args[x2 + '0'] = origin[1] - height/2 + corner_radius
            self.x1_max_x2_min = cyl(name=f'{x1} max {x2} min', **args)

            args[x1 + '0'] = origin[0] + width/2 - corner_radius
            args[x2 + '0'] = origin[1] + height/2 - corner_radius
            self.x1_max_x2_max = cyl(name=f'{x1} max {x2} max', **args)

            self.x1_min = _plane(x1, 'min', -width/2 + origin[0] + corner_radius,
                                 **bc_args)
            self.x1_max = _plane(x1, 'max', width/2 + origin[0] - corner_radius,
                                 **bc_args)
            self.x2_min = _plane(x2, 'min', -height/2 + origin[1] + corner_radius,
                                 **bc_args)
            self.x2_max = _plane(x2, 'max', height/2 + origin[1] - corner_radius,
                                 **bc_args)
            self._surface_names += (
                'x1_min_x2_min', 'x1_min_x2_max', 'x1_max_x2_min',
                'x1_max_x2_max', 'x1_min', 'x1_max', 'x2_min', 'x2_max'
            )

    def __neg__(self):
        prism = +self.min_x1 & -self.max_x1 & +self.min_x2 & -self.max_x2

        # Cut out corners if a corner radius was given
        if hasattr(self, 'x1_min'):
            corners = (
                (+self.x1_min_x2_min & -self.x1_min & -self.x2_min) |
                (+self.x1_min_x2_max & -self.x1_min & +self.x2_max) |
                (+self.x1_max_x2_min & +self.x1_max & -self.x2_min) |
                (+self.x1_max_x2_max & +self.x1_max & +self.x2_max)
            )
            prism &= ~corners

        return prism


class HexagonalPrism(CompositeSurface):
    """Hexagonal prism comoposed of six planar surfaces

    .. versionadded:: 0.14.0

    Parameters
    ----------
    edge_length : float
        Length of a side of the hexagon in [cm]
    orientation : {'x', 'y'}
        An 'x' orientation means that two sides of the hexagon are parallel to
        the x-axis and a 'y' orientation means that two sides of the hexagon are
        parallel to the y-axis.
    origin : Iterable of two floats
        Origin of the prism.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}
        Boundary condition that defines the behavior for particles hitting the
        surfaces comprising the hexagonal prism.
    albedo : float, optional
        Albedo of the prism's surfaces as a ratio of particle weight after
        interaction with the surface to the initial weight. Values must be
        positive. Only applicable if the boundary type is 'reflective',
        'periodic', or 'white'.
    corner_radius : float
        Prism corner radius in units of [cm].

    """
    _surface_names = ('plane_max', 'plane_min', 'upper_right', 'upper_left',
                      'lower_right', 'lower_left')

    def __init__(
            self,
            edge_length: float = 1.,
            orientation: str = 'y',
            origin: Sequence[float] = (0., 0.),
            boundary_type: str = 'transmission',
            albedo: float = 1.,
            corner_radius: float = 0.
    ):
        check_type('edge_length', edge_length, Real)
        check_type('albedo', albedo, Real)
        check_type('corner_radius', corner_radius, Real)
        check_value('orientation', orientation, ('x', 'y'))
        check_type('origin', origin, Iterable, Real)

        l = edge_length
        x, y = origin

        # Create container for boundary arguments
        bc_args = {'boundary_type': boundary_type, 'albedo': albedo}

        if orientation == 'y':
            # Left and right planes
            self.plane_max = openmc.XPlane(x + sqrt(3.)/2*l, **bc_args)
            self.plane_min = openmc.XPlane(x - sqrt(3.)/2*l, **bc_args)
            c = sqrt(3.)/3.

            # y = -x/sqrt(3) + a
            self.upper_right = openmc.Plane(a=c, b=1., d=l+x*c+y, **bc_args)

            # y = x/sqrt(3) + a
            self.upper_left = openmc.Plane(a=-c, b=1., d=l-x*c+y, **bc_args)

            # y = x/sqrt(3) - a
            self.lower_right = openmc.Plane(a=-c, b=1., d=-l-x*c+y, **bc_args)

            # y = -x/sqrt(3) - a
            self.lower_left = openmc.Plane(a=c, b=1., d=-l+x*c+y, **bc_args)

        elif orientation == 'x':
            self.plane_max = openmc.YPlane(y + sqrt(3.)/2*l, **bc_args)
            self.plane_min = openmc.YPlane(y - sqrt(3.)/2*l, **bc_args)
            c = sqrt(3.)

            # Upper-right surface: y = -sqrt(3)*(x - a)
            self.upper_right = openmc.Plane(a=c, b=1., d=c*l+x*c+y, **bc_args)

            # Lower-right surface: y = sqrt(3)*(x + a)
            self.lower_right = openmc.Plane(a=-c, b=1., d=-c*l-x*c+y, **bc_args)

            # Lower-left surface: y = -sqrt(3)*(x + a)
            self.lower_left = openmc.Plane(a=c, b=1., d=-c*l+x*c+y, **bc_args)

            # Upper-left surface: y = sqrt(3)*(x + a)
            self.upper_left = openmc.Plane(a=-c, b=1., d=c*l-x*c+y, **bc_args)

        # Handle periodic boundary conditions
        if boundary_type == 'periodic':
            self.plane_min.periodic_surface = self.plane_max
            self.upper_right.periodic_surface = self.lower_left
            self.lower_right.periodic_surface = self.upper_left

        # Handle rounded corners if given
        if corner_radius > 0.:
            if boundary_type == 'periodic':
                raise ValueError('Periodic boundary conditions not permitted '
                                 'when rounded corners are used.')

            c = sqrt(3.)/2
            t = l - corner_radius/c

            # Cylinder with corner radius and boundary type pre-applied
            cyl1 = partial(openmc.ZCylinder, r=corner_radius, **bc_args)
            cyl2 = partial(openmc.ZCylinder, r=corner_radius/(2*c), **bc_args)

            if orientation == 'x':
                self.x_min_y_min_in = cyl1(name='x min y min in', x0=x-t/2, y0=y-c*t)
                self.x_min_y_max_in = cyl1(name='x min y max in', x0=x+t/2, y0=y-c*t)
                self.x_max_y_min_in = cyl1(name='x max y min in', x0=x-t/2, y0=y+c*t)
                self.x_max_y_max_in = cyl1(name='x max y max in', x0=x+t/2, y0=y+c*t)
                self.min_in = cyl1(name='x min in', x0=x-t, y0=y)
                self.max_in = cyl1(name='x max in', x0=x+t, y0=y)

                self.x_min_y_min_out = cyl2(name='x min y min out', x0=x-l/2, y0=y-c*l)
                self.x_min_y_max_out = cyl2(name='x min y max out', x0=x+l/2, y0=y-c*l)
                self.x_max_y_min_out = cyl2(name='x max y min out', x0=x-l/2, y0=y+c*l)
                self.x_max_y_max_out = cyl2(name='x max y max out', x0=x+l/2, y0=y+c*l)
                self.min_out = cyl2(name='x min out', x0=x-l, y0=y)
                self.max_out = cyl2(name='x max out', x0=x+l, y0=y)

            elif orientation == 'y':
                self.x_min_y_min_in = cyl1(name='x min y min in', x0=x-c*t, y0=y-t/2)
                self.x_min_y_max_in = cyl1(name='x min y max in', x0=x-c*t, y0=y+t/2)
                self.x_max_y_min_in = cyl1(name='x max y min in', x0=x+c*t, y0=y-t/2)
                self.x_max_y_max_in = cyl1(name='x max y max in', x0=x+c*t, y0=y+t/2)
                self.min_in = cyl1(name='y min in', x0=x, y0=y-t)
                self.max_in = cyl1(name='y max in', x0=x, y0=y+t)

                self.x_min_y_min_out = cyl2(name='x min y min out', x0=x-c*l, y0=y-l/2)
                self.x_min_y_max_out = cyl2(name='x min y max out', x0=x-c*l, y0=y+l/2)
                self.x_max_y_min_out = cyl2(name='x max y min out', x0=x+c*l, y0=y-l/2)
                self.x_max_y_max_out = cyl2(name='x max y max out', x0=x+c*l, y0=y+l/2)
                self.min_out = cyl2(name='y min out', x0=x, y0=y-l)
                self.max_out = cyl2(name='y max out', x0=x, y0=y+l)

            # Add to tuple of surface names
            for s in ('in', 'out'):
                self._surface_names += (
                    f'x_min_y_min_{s}', f'x_min_y_max_{s}',
                    f'x_max_y_min_{s}', f'x_max_y_max_{s}',
                    f'min_{s}', f'max_{s}')

    def __neg__(self) -> openmc.Region:
        prism = (
            -self.plane_max & +self.plane_min &
            -self.upper_right & -self.upper_left &
            +self.lower_right & +self.lower_left
        )

        # Cut out corners if a corner radius was given
        if hasattr(self, 'min_in'):
            corners = (
                +self.x_min_y_min_in & -self.x_min_y_min_out |
                +self.x_min_y_max_in & -self.x_min_y_max_out |
                +self.x_max_y_min_in & -self.x_max_y_min_out |
                +self.x_max_y_max_in & -self.x_max_y_max_out |
                +self.min_in & -self.min_out |
                +self.max_in & -self.max_out
            )
            prism &= ~corners

        return prism
