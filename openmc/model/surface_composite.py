from abc import ABC, abstractmethod
from copy import copy

import openmc
from openmc.checkvalue import check_greater_than, check_value
from numpy import sqrt

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

class IsogonalOctagon(CompositeSurface):
    """Infinite isogonal octagon composite surface

    An isogonal octagon is composed of eight planar surfaces. The prism is
    parallel to the x, y, or z axis. The remaining two axes (y and z, x and z,
    or x and y) serve as a basis for constructing the surfaces. Two surfaces
    are parallel to the first basis axis, two surfaces are parallel
    to the second basis axis, and the remaining four surfaces intersect both
    basis axes at 45 degree angles.

    This class acts as a proper surface, meaning that unary `+` and `-`
    operators applied to it will produce a half-space. The negative side is
    defined to be the region inside of the octogonal prism.

    Parameters
    ----------
    center : iterable of float
        Coordinate for the central axis of the octagon in the
        (y, z), (x, z), or (x, y) basis.
    r1 : float
        Half-width of octagon across its basis axis-parallel sides in units
        of cm
    r2 : float
        Half-width of octagon across its basis axis intersecting sides in
        units of cm. Must be less than :math:`\sqrt{2} r_1`.
    axis : {'x', 'y', 'z'}
        Central axis of octagon. Defaults to 'z'
    **kwargs
        Keyword arguments passed to underlying plane classes

    Attributes
    ----------
    top : openmc.ZPlane or openmc.YPlane
        Top planar surface of octagon
    bottom : openmc.ZPlane or openmc.YPlane
        Bottom planar surface of octagon
    right: openmc.YPlane or openmc.XPlane
        Right planaer surface of octagon
    left: openmc.YPlane or openmc.XPlane
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
        self._axis = axis
        q1c, q2c = center

        # Coords for axis-perpendicular planes
        ctop = q1c+r1
        cbottom = q1c-r1

        cright= q2c+r1
        cleft = q2c-r1

        # Side lengths
        if (r2 > r1 * sqrt(2)):
            raise ValueError(f'r2 must be less than {sqrt(2) * r1}')

        L_perp_ax1 = (r2 * sqrt(2) - r1) * 2
        L_perp_ax2 = (r1 * sqrt(2) - r2) * 2

        # Coords for quadrant planes
        p1_ur = [L_perp_ax1/2, r1, 0]
        p2_ur = [r1, L_perp_ax1/2, 0]
        p3_ur = [r1, L_perp_ax1/2, 1]

        p1_lr = [r1, -L_perp_ax1/2, 0]
        p2_lr = [L_perp_ax1/2, -r1, 0]
        p3_lr = [L_perp_ax1/2, -r1, 1]

        points = [p1_ur, p2_ur, p3_ur, p1_lr, p2_lr, p3_lr]

        # Orientation specific variables
        if axis == 'z':
            coord_map = [0,1,2]
            self.top = openmc.YPlane(y0=ctop, **kwargs)
            self.bottom = openmc.YPlane(y0=cbottom, **kwargs)
            self.right = openmc.XPlane(x0=cright, **kwargs)
            self.left = openmc.XPlane(x0=cleft, **kwargs)
        elif axis == 'y':
            coord_map = [0,2,1]
            self.top = openmc.ZPlane(z0=ctop, **kwargs)
            self.bottom = openmc.ZPlane(z0=cbottom, **kwargs)
            self.right = openmc.XPlane(x0=cright, **kwargs)
            self.left = openmc.XPlane(x0=cleft, **kwargs)
        elif axis == 'x':
            coord_map = [2,0,1]
            self.top = openmc.ZPlane(z0=ctop, **kwargs)
            self.bottom = openmc.ZPlane(z0=cbottom, **kwargs)
            self.right = openmc.YPlane(y0=cright, **kwargs)
            self.left = openmc.YPlane(y0=cleft, **kwargs)

        # Put our coordinates in (x,y,z) order
        calibrated_points = []
        for p in points:
            p_temp = []
            for i in coord_map:
                p_temp += [p[i]]
            calibrated_points += [p_temp]

        p1_ur, p2_ur, p3_ur, p1_lr, p2_lr, p3_lr = calibrated_points

        self.upper_right = openmc.Plane.from_points(p1_ur, p2_ur, p3_ur,
                                                    **kwargs)
        self.lower_right = openmc.Plane.from_points(p1_lr, p2_lr, p3_lr,
                                                    **kwargs)
        self.lower_left = openmc.Plane.from_points(-p1_ur, -p2_ur, -p3_ur,
                                                   **kwargs)
        self.upper_left = openmc.Plane.from_points(-p1_lr, -p2_lr, -p3_lr,
                                                   **kwargs)


    def __neg__(self):
        region1 = -self.top & +self.bottom & -self.right & +self.left
        if self._axis == 'y':
            region2 = -self.upper_right & -self.lower_right & \
                +self.lower_left & +self.upper_left
        else:
            region2 = +self.upper_right & +self.lower_right & \
                -self.lower_left & -self.upper_left

        return region1 & region2


    def __pos__(self):
        region1 = +self.top | -self.bottom | +self.right | -self.left
        if self._axis == 'y':
            region2 = +self.upper_right | +self.lower_right | \
                -self.lower_left | -self.upper_left
        else:
            region2 = -self.upper_right | -self.lower_right | \
                +self.lower_left | +self.upper_left

        return region1 | region2


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
        return +self.xmin & -self.xmax & +self.ymin & -self.ymax & +self.zmin & -self.zmax

    def __pos__(self):
        return -self.xmin | +self.xmax | -self.ymin | +self.ymax | -self.zmin | +self.zmax


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
