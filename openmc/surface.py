from abc import ABCMeta
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

from openmc.checkvalue import check_type, check_value, check_greater_than
from openmc.constants import BC_TYPES

if sys.version_info[0] >= 3:
    basestring = str

# A static variable for auto-generated Surface IDs
AUTO_SURFACE_ID = 10000


def reset_auto_surface_id():
    global AUTO_SURFACE_ID
    AUTO_SURFACE_ID = 10000


class Surface(object):
    """A two-dimensional surface that can be used define regions of space with an
    associated boundary condition.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the surface. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface, e.g. 'x-plane'
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coeffs : dict
        Dictionary of surface coefficients

    """

    def __init__(self, surface_id=None, boundary_type='transmission', name=''):
        # Initialize class attributes
        self.id = surface_id
        self.name = name
        self._type = ''
        self.boundary_type = boundary_type

        # A dictionary of the quadratic surface coefficients
        # Key        - coefficeint name
        # Value    - coefficient value
        self._coeffs = {}

        # An ordered list of the coefficient names to export to XML in the
        # proper order
        self._coeff_keys = []

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def type(self):
        return self._type

    @property
    def boundary_type(self):
        return self._boundary_type

    @property
    def coeffs(self):
        return self._coeffs

    @id.setter
    def id(self, surface_id):
        if surface_id is None:
            global AUTO_SURFACE_ID
            self._id = AUTO_SURFACE_ID
            AUTO_SURFACE_ID += 1
        else:
            check_type('surface ID', surface_id, Integral)
            check_greater_than('surface ID', surface_id, 0)
            self._id = surface_id

    @name.setter
    def name(self, name):
        check_type('surface name', name, basestring)
        self._name = name

    @boundary_type.setter
    def boundary_type(self, boundary_type):
        check_type('boundary type', boundary_type, basestring)
        check_value('boundary type', boundary_type, BC_TYPES.values())
        self._boundary_type = boundary_type

    def __repr__(self):
        string = 'Surface\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
        string += '{0: <16}{1}{2}\n'.format('\tBoundary', '=\t', self._boundary_type)

        coeffs = '{0: <16}'.format('\tCoefficients') + '\n'

        for coeff in self._coeffs:
            coeffs += '{0: <16}{1}{2}\n'.format(coeff, '=\t', self._coeffs[coeff])

        string += coeffs

        return string

    def create_xml_subelement(self):
        element = ET.Element("surface")
        element.set("id", str(self._id))

        if len(self._name) > 0:
            element.set("name", str(self._name))

        element.set("type", self._type)
        element.set("boundary", self._boundary_type)
        element.set("coeffs", ' '.join([str(self._coeffs[key])
                                        for key in self._coeff_keys]))

        return element


class Plane(Surface):
    """An arbitrary plane of the form :math:`Ax + By + Cz = D`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    A : float
        The 'A' parameter for the plane
    B : float
        The 'B' parameter for the plane
    C : float
        The 'C' parameter for the plane
    D : float
        The 'D' parameter for the plane
    name : str
        Name of the plane. If not specified, the name will be the empty string.

    Attributes
    ----------
    a : float
        The 'A' parameter for the plane
    b : float
        The 'B' parameter for the plane
    c : float
        The 'C' parameter for the plane
    d : float
        The 'D' parameter for the plane

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 A=None, B=None, C=None, D=None, name=''):
        # Initialize Plane class attributes
        super(Plane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'plane'
        self._coeff_keys = ['A', 'B', 'C', 'D']

        if A is not None:
            self.a = A

        if B is not None:
            self.b = B

        if C is not None:
            self.c = C

        if D is not None:
            self.d = D

    @property
    def a(self):
        return self.coeffs['A']

    @property
    def b(self):
        return self.coeffs['B']

    @property
    def c(self):
        return self.coeffs['C']

    @property
    def d(self):
        return self.coeffs['D']

    @a.setter
    def a(self, A):
        check_type('A coefficient', A, Real)
        self._coeffs['A'] = A

    @b.setter
    def b(self, B):
        check_type('B coefficient', B, Real)
        self._coeffs['B'] = B

    @c.setter
    def c(self, C):
        check_type('C coefficient', C, Real)
        self._coeffs['C'] = C

    @d.setter
    def d(self, D):
        check_type('D coefficient', D, Real)
        self._coeffs['D'] = D


class XPlane(Plane):
    """A plane perpendicular to the x axis, i.e. a surface of the form :math:`x -
    x_0 = 0`

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        Location of the plane
    name : str
        Name of the plane. If not specified, the name will be the empty string.

    Attributes
    ----------
    x0 : float
        Location of the plane

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, name=''):
        # Initialize XPlane class attributes
        super(XPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'x-plane'
        self._coeff_keys = ['x0']

        if x0 is not None:
            self.x0 = x0

    @property
    def x0(self):
        return self.coeff['x0']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coeffs['x0'] = x0


class YPlane(Plane):
    """A plane perpendicular to the y axis, i.e. a surface of the form :math:`y -
    y_0 = 0`

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    y0 : float
        Location of the plane
    name : str
        Name of the plane. If not specified, the name will be the empty string.

    Attributes
    ----------
    y0 : float
        Location of the plane

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 y0=None, name=''):
        # Initialize YPlane class attributes
        super(YPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'y-plane'
        self._coeff_keys = ['y0']

        if y0 is not None:
            self.y0 = y0

    @property
    def y0(self):
        return self.coeffs['y0']

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coeffs['y0'] = y0


class ZPlane(Plane):
    """A plane perpendicular to the z axis, i.e. a surface of the form :math:`z -
    z_0 = 0`

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    z0 : float
        Location of the plane
    name : str
        Name of the plane. If not specified, the name will be the empty string.

    Attributes
    ----------
    z0 : float
        Location of the plane

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 z0=None, name=''):
        # Initialize ZPlane class attributes
        super(ZPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'z-plane'
        self._coeff_keys = ['z0']

        if z0 is not None:
            self.z0 = z0

    @property
    def z0(self):
        return self.coeffs['z0']

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coeffs['z0'] = z0


class Cylinder(Surface):
    """A cylinder whose length is parallel to the x-, y-, or z-axis.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    R : float
        Radius of the cylinder
    name : str
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    r : float
        Radius of the cylinder

    """

    __metaclass__ = ABCMeta

    def __init__(self, surface_id=None, boundary_type='transmission',
                 R=None, name=''):
        # Initialize Cylinder class attributes
        super(Cylinder, self).__init__(surface_id, boundary_type, name=name)

        self._coeff_keys = ['R']

        if R is not None:
            self.r = R

    @property
    def r(self):
        return self.coeffs['R']

    @r.setter
    def r(self, R):
        check_type('R coefficient', R, Real)
        self._coeffs['R'] = R


class XCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the x-axis. This is a
    quadratic surface of the form :math:`(y - y_0)^2 + (z - z_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    y0 : float
        y-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder
    R : float
        Radius of the cylinder
    name : str
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    y0 : float
        y-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 y0=None, z0=None, R=None, name=''):
        # Initialize XCylinder class attributes
        super(XCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'x-cylinder'
        self._coeff_keys = ['y0', 'z0', 'R']

        if y0 is not None:
            self.y0 = y0

        if z0 is not None:
            self.z0 = z0

    @property
    def y0(self):
        return self.coeffs['y0']

    @property
    def z0(self):
        return self.coeffs['z0']

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coeffs['y0'] = y0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coeffs['z0'] = z0


class YCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the y-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (z - z_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        x-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder
    R : float
        Radius of the cylinder
    name : str
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, z0=None, R=None, name=''):
        # Initialize YCylinder class attributes
        super(YCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'y-cylinder'
        self._coeff_keys = ['x0', 'z0', 'R']

        if x0 is not None:
            self.x0 = x0

        if z0 is not None:
            self.z0 = z0

    @property
    def x0(self):
        return self.coeffs['x0']

    @property
    def z0(self):
        return self.coeffs['z0']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coeffs['x0'] = x0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coeffs['z0'] = z0


class ZCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the z-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (y - y_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        x-coordinate of the center of the cylinder
    y0 : float
        y-coordinate of the center of the cylinder
    R : float
        Radius of the cylinder
    name : str
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the center of the cylinder
    y0 : float
        y-coordinate of the center of the cylinder

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, R=None, name=''):
        # Initialize ZCylinder class attributes
        super(ZCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'z-cylinder'
        self._coeff_keys = ['x0', 'y0', 'R']

        if x0 is not None:
            self.x0 = x0

        if y0 is not None:
            self.y0 = y0

    @property
    def x0(self):
        return self.coeffs['x0']

    @property
    def y0(self):
        return self.coeffs['y0']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coeffs['x0'] = x0

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coeffs['y0'] = y0


class Sphere(Surface):
    """A sphere of the form :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        x-coordinate of the center of the sphere
    y0 : float
        y-coordinate of the center of the sphere
    z0 : float
        z-coordinate of the center of the sphere
    R : float
        Radius of the sphere
    name : str
        Name of the sphere. If not specified, the name will be the empty string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the center of the sphere
    y0 : float
        y-coordinate of the center of the sphere
    z0 : float
        z-coordinate of the center of the sphere
    R : float
        Radius of the sphere

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R=None, name=''):
        # Initialize Sphere class attributes
        super(Sphere, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'sphere'
        self._coeff_keys = ['x0', 'y0', 'z0', 'R']

        if x0 is not None:
            self.x0 = x0

        if y0 is not None:
            self.y0 = y0

        if z0 is not None:
            self.z0 = z0

        if R is not None:
            self.r = R

    @property
    def x0(self):
        return self.coeffs['x0']

    @property
    def y0(self):
        return self.coeffs['y0']

    @property
    def z0(self):
        return self.coeffs['z0']

    @property
    def r(self):
        return self.coeffs['R']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coeffs['x0'] = x0

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coeffs['y0'] = y0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coeffs['z0'] = z0

    @r.setter
    def r(self, R):
        check_type('R coefficient', R, Real)
        self._coeffs['R'] = R


class Cone(Surface):
    """A conical surface parallel to the x-, y-, or z-axis.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature
    name : str
        Name of the cone. If not specified, the name will be the empty string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature

    """

    __metaclass__ = ABCMeta

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):
        # Initialize Cone class attributes
        super(Cone, self).__init__(surface_id, boundary_type, name=name)

        self._coeff_keys = ['x0', 'y0', 'z0', 'R2']

        if x0 is not None:
            self.x0 = x0

        if y0 is not None:
            self.y0 = y0

        if z0 is not None:
            self.z0 = z0

        if R2 is not None:
            self.r2 = R2

    @property
    def x0(self):
        return self.coeffs['x0']

    @property
    def y0(self):
        return self.coeffs['y0']

    @property
    def z0(self):
        return self.coeffs['z0']

    @property
    def r2(self):
        return self.coeffs['r2']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coeffs['x0'] = x0

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coeffs['y0'] = y0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coeffs['z0'] = z0

    @r2.setter
    def r2(self, R2):
        check_type('R^2 coefficient', R2, Real)
        self._coeffs['R2'] = R2


class XCone(Cone):
    """A cone parallel to the x-axis of the form :math:`(y - y_0)^2 + (z - z_0)^2 =
    R^2 (x - x_0)^2`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature
    name : str
        Name of the cone. If not specified, the name will be the empty string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):
        # Initialize XCone class attributes
        super(XCone, self).__init__(surface_id, boundary_type, x0, y0,
                                    z0, R2, name=name)

        self._type = 'x-cone'


class YCone(Cone):
    """A cone parallel to the y-axis of the form :math:`(x - x_0)^2 + (z - z_0)^2 =
    R^2 (y - y_0)^2`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature
    name : str
        Name of the cone. If not specified, the name will be the empty string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):
        # Initialize YCone class attributes
        super(YCone, self).__init__(surface_id, boundary_type, x0, y0, z0,
                                    R2, name=name)

        self._type = 'y-cone'


class ZCone(Cone):
    """A cone parallel to the x-axis of the form :math:`(x - x_0)^2 + (y - y_0)^2 =
    R^2 (z - z_0)^2`.

    Parameters
    ----------
    surface_id : int
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature
    name : str
        Name of the cone. If not specified, the name will be the empty string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    R2 : float
        Parameter related to the aperature

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):
        # Initialize ZCone class attributes
        super(ZCone, self).__init__(surface_id, boundary_type, x0, y0, z0,
                                    R2, name=name)

        self._type = 'z-cone'
