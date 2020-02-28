from copy import copy

import openmc


class CompositeMixin:
    def evaluate(self, point):
        raise NotImplementedError('Composite surfaces do not have a surface equation.')

    def _get_base_coeffs(self):
        raise NotImplementedError('Composite surfaces do not have base coefficients.')

    def translate(self, vector, inplace=False):
        surf = self if inplace else copy(self)
        for name in self._surface_names:
            s = getattr(surf, name)
            setattr(surf, name, s.translate(vector, inplace))
        return surf

    def rotate(self, rotation):
        surf = copy(self)
        for name in self._surface_names:
            s = getattr(surf, name)
            setattr(surf, name, s.rotate(rotation))
        return surf

    def __repr__(self):
        return "<{} at 0x{:x}>".format(type(self).__name__, id(self))


class RightCircularCylinder(CompositeMixin, openmc.Surface):
    _surface_names = ('cyl', 'bottom', 'top')

    def __init__(self, center_base, height, radius, axis='z',
                 boundary_type='transmission'):
        kwargs = {'boundary_type': boundary_type}
        cx, cy, cz = center_base
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


class RectangularParallelepiped(CompositeMixin, openmc.Surface):
    _surface_names = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax')

    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax, boundary_type='transmission'):
        kwargs = {'boundary_type': boundary_type}
        self.xmin = openmc.XPlane(x0=xmin, **kwargs)
        self.xmax = openmc.XPlane(x0=xmax, **kwargs)
        self.ymin = openmc.YPlane(y0=ymin, **kwargs)
        self.ymax = openmc.YPlane(y0=ymax, **kwargs)
        self.zmin = openmc.ZPlane(z0=zmin, **kwargs)
        self.zmax = openmc.ZPlane(z0=zmax, **kwargs)

    def __neg__(self):
        return +self.xmin & -self.xmax & +self.ymin & -self.ymax & +self.zmin & -self.zmax

    def __pos__(self):
        return -self.xmin | +self.ymax | -self.ymin | +self.ymax | -self.zmin | +self.zmax


class Box(CompositeMixin, openmc.Surface):
    _surface_names = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax')

    def __init__(self, v, a1, a2, a3, boundary_type='transmission'):
        kwargs = {'boundary_type': boundary_type}
        vx, vy, vz = v
        a1x, a1y, a1z = a1
        a2x, a2y, a2z = a2
        a3x, a3y, a3z = a3

        # Only support boxes with axis-aligned vectors
        if any(x != 0.0 for x in (a1y, a1z, a2x, a2z, a3x, a3y)):
            raise NotImplementedError('Box macrobody with non-axis-aligned '
                                      'vector not supported.')

        # Determine each side of the box
        if a1x > 0:
            xmin, xmax = vx, vx + a1x
        else:
            xmin, xmax = vx + a1x, vx
        if a2y > 0:
            ymin, ymax = vy, vy + a2y
        else:
            ymin, ymax = vy + a2y, vy
        if a3z > 0:
            zmin, zmax = vz, vz + a3z
        else:
            zmin, zmax = vz + a3z, vz

        # Create surfaces
        self.xmin = openmc.XPlane(x0=xmin, **kwargs)
        self.xmax = openmc.XPlane(x0=xmax, **kwargs)
        self.ymin = openmc.YPlane(y0=ymin, **kwargs)
        self.ymax = openmc.YPlane(y0=ymax, **kwargs)
        self.zmin = openmc.ZPlane(z0=zmin, **kwargs)
        self.zmax = openmc.ZPlane(z0=zmax, **kwargs)

    def __neg__(self):
        return (+self.xmin & -self.xmax &
                +self.ymin & -self.ymax &
                +self.zmin & -self.zmax)

    def __pos__(self):
        return (-self.xmin | +self.xmax |
                -self.ymin | +self.ymax |
                -self.zmin | +self.zmax)


class XConeOneSided(CompositeMixin, openmc.Surface):
    _surface_names = ('cone', 'plane')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., up=True, **kwargs):
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


class YConeOneSided(CompositeMixin, openmc.Surface):
    _surface_names = ('cone', 'plane')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., up=True, **kwargs):
        self.cone = openmc.YCone(x0, y0, z0, r2, **kwargs)
        self.plane = openmc.YPlane(y0)
        self.up = up

    __neg__ = XConeOneSided.__neg__
    __pos__ = XConeOneSided.__pos__


class ZConeOneSided(CompositeMixin, openmc.Surface):
    _surface_names = ('cone', 'plane')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., up=True, **kwargs):
        self.cone = openmc.ZCone(x0, y0, z0, r2, **kwargs)
        self.plane = openmc.ZPlane(z0)
        self.up = up

    __neg__ = XConeOneSided.__neg__
    __pos__ = XConeOneSided.__pos__
