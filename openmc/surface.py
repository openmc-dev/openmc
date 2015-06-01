from xml.etree import ElementTree as ET

from openmc.checkvalue import *
from openmc.constants import BC_TYPES


# A static variable for auto-generated Surface IDs
AUTO_SURFACE_ID = 10000

def reset_auto_surface_id():
    global AUTO_SURFACE_ID
    AUTO_SURFACE_ID = 10000



class Surface(object):

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

        # Check that the ID is an integer and wasn't already used
        elif not is_integer(surface_id):
            msg = 'Unable to set a non-integer Surface ' \
                  'ID {0}'.format(surface_id)
            raise ValueError(msg)

        elif surface_id < 0:
            msg = 'Unable to set Surface ID to {0} since it must be a ' \
                  'non-negative integer'.format(surface_id)
            raise ValueError(msg)

        else:
            self._id = surface_id


    @name.setter
    def name(self, name):

        if not is_string(name):
            msg = 'Unable to set name for Surface ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name


    @boundary_type.setter
    def boundary_type(self, boundary_type):

        if not is_string(boundary_type):
            msg = 'Unable to set boundary type for Surface ID={0} with a ' \
                  'non-string value {1}'.format(self._id, boundary_type)
            raise ValueError(msg)

        elif not boundary_type in BC_TYPES.values():
            msg = 'Unable to set boundary type for Surface ID={0} to ' \
                  '{1} which is not trasmission, vacuum or ' \
                  'reflective'.format(boundary_type)
            raise ValueError(msg)

        else:
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

        coeffs = ''

        for coeff in self._coeff_keys:
            coeffs += '{0} '.format(self._coeffs[coeff])

        element.set("coeffs", coeffs.rstrip(' '))

        return element



class Plane(Surface):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 A=None, B=None, C=None, D=None, name='',):

        # Initialize Plane class attributes
        super(Plane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'plane'
        self._coeff_keys = ['A', 'B', 'C', 'D']

        if not A is None:
            self.a = A

        if not B is None:
            self.b = B

        if not C is None:
            self.c = C

        if not D is None:
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

        if not is_integer(A) and not is_float(A):
            msg = 'Unable to set A coefficient for Plane ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, A)
            raise ValueError(msg)

        self._coeffs['A'] = A


    @b.setter
    def b(self, B):

        if not is_integer(B) and not is_float(B):
            msg = 'Unable to set B coefficient for Plane ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, B)
            raise ValueError(msg)

        self._coeffs['B'] = B


    @c.setter
    def c(self, C):

        if not is_integer(C) and not is_float(C):
            msg = 'Unable to set C coefficient for Plane ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, C)
            raise ValueError(msg)

        self._coeffs['C'] = C


    @d.setter
    def d(self, D):

        if not is_integer(D) and not is_float(D):
            msg = 'Unable to set D coefficient for Plane ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, D)
            raise ValueError(msg)

        self._coeffs['D'] = D



class XPlane(Plane):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, name=''):

        # Initialize XPlane class attributes
        super(XPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'x-plane'
        self._coeff_keys = ['x0']

        if not x0 is None:
            self.x0 = x0


    @property
    def x0(self):
        return self.coeff['x0']


    @x0.setter
    def x0(self, x0):

        if not is_integer(x0) and not is_float(x0):
            msg = 'Unable to set x0 coefficient for XPlane ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, x0)
            raise ValueError(msg)

        self._coeffs['x0'] = x0



class YPlane(Plane):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 y0=None, name=''):

        # Initialize YPlane class attributes
        super(YPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'y-plane'
        self._coeff_keys = ['y0']

        if not y0 is None:
            self.y0 = y0


    @property
    def y0(self):
        return self.coeffs['y0']


    @y0.setter
    def y0(self, y0):

        if not is_integer(y0) and not is_float(y0):
            msg = 'Unable to set y0 coefficient for XPlane ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, y0)
            raise ValueError(msg)

        self._coeffs['y0'] = y0



class ZPlane(Plane):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 z0=None, name=''):

        # Initialize ZPlane class attributes
        super(ZPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'z-plane'
        self._coeff_keys = ['z0']

        if not z0 is None:
            self.z0 = z0


    @property
    def z0(self):
        return self.coeffs['z0']


    @z0.setter
    def z0(self, z0):

        if not is_integer(z0) and not is_float(z0):
            msg = 'Unable to set z0 coefficient for ZPlane ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, z0)
            raise ValueError(msg)

        self._coeffs['z0'] = z0



class Cylinder(Surface):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 R=None, name=''):

        # Initialize Cylinder class attributes
        super(Cylinder, self).__init__(surface_id, boundary_type, name=name)

        self._coeff_keys = ['R']

        if not R is None:
            self.r = R


    @property
    def r(self):
        return self.coeffs['R']


    @r.setter
    def r(self, R):

        if not is_integer(R) and not is_float(R):
            msg = 'Unable to set R coefficient for Cylinder ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, R)
            raise ValueError(msg)

        self._coeffs['R'] = R



class XCylinder(Cylinder):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 y0=None, z0=None, R=None, name=''):

        # Initialize XCylinder class attributes
        super(XCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'x-cylinder'
        self._coeff_keys = ['y0', 'z0', 'R']

        if not y0 is None:
            self.y0 = y0

        if not z0 is None:
            self.z0 = z0


    @property
    def y0(self):
        return self.coeffs['y0']


    @property
    def z0(self):
        return self.coeffs['z0']


    @y0.setter
    def y0(self, y0):

        if not is_integer(y0) and not is_float(y0):
            msg = 'Unable to set y0 coefficient for XCylinder ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, y0)
            raise ValueError(msg)

        self._coeffs['y0'] = y0


    @z0.setter
    def z0(self, z0):

        if not is_integer(z0) and not is_float(z0):
            msg = 'Unable to set z0 coefficient for XCylinder ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, z0)
            raise ValueError(msg)

        self._coeffs['z0'] = z0



class YCylinder(Cylinder):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, z0=None, R=None, name=''):

        # Initialize YCylinder class attributes
        super(YCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'y-cylinder'
        self._coeff_keys = ['x0', 'z0', 'R']

        if not x0 is None:
            self.x0 = x0

        if not z0 is None:
            self.z0 = z0


    @property
    def x0(self):
        return self.coeffs['x0']


    @property
    def z0(self):
        return self.coeffs['z0']


    @x0.setter
    def x0(self, x0):

        if not is_integer(x0) and not is_float(x0):
            msg = 'Unable to set x0 coefficient for YCylinder ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, x0)
            raise ValueError(msg)

        self._coeffs['x0'] = x0


    @z0.setter
    def z0(self, z0):

        if not is_integer(z0) and not is_float(z0):
            msg = 'Unable to set z0 coefficient for YCylinder ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, z0)
            raise ValueError(msg)

        self._coeffs['z0'] = z0


class ZCylinder(Cylinder):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, R=None, name=''):

        # Initialize ZCylinder class attributes
        # Initialize YPlane class attributes
        super(ZCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'z-cylinder'
        self._coeff_keys = ['x0', 'y0', 'R']

        if not x0 is None:
            self.x0 = x0

        if not y0 is None:
            self.y0 = y0


    @property
    def x0(self):
        return self.coeffs['x0']


    @property
    def y0(self):
        return self.coeffs['y0']


    @x0.setter
    def x0(self, x0):

        if not is_integer(x0) and not is_float(x0):
            msg = 'Unable to set x0 coefficient for ZCylinder ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, x0)
            raise ValueError(msg)

        self._coeffs['x0'] = x0


    @y0.setter
    def y0(self, y0):

        if not is_integer(y0) and not is_float(y0):
            msg = 'Unable to set y0 coefficient for ZCylinder ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, y0)
            raise ValueError(msg)

        self._coeffs['y0'] = y0



class Sphere(Surface):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R=None, name=''):

        # Initialize Sphere class attributes
        super(Sphere, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'sphere'
        self._coeff_keys = ['x0', 'y0', 'z0', 'R']

        if not x0 is None:
            self.x0 = x0

        if not y0 is None:
            self.y0 = y0

        if not z0 is None:
            self.z0 = z0

        if not R is None:
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

        if not is_integer(x0) and not is_float(x0):
            msg = 'Unable to set x0 coefficient for Sphere ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, x0)
            raise ValueError(msg)

        self._coeffs['x0'] = x0


    @y0.setter
    def y0(self, y0):

        if not is_integer(y0) and not is_float(y0):
            msg = 'Unable to set y0 coefficient for Sphere ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, y0)
            raise ValueError(msg)

        self._coeffs['y0'] = y0


    @z0.setter
    def z0(self, z0):

        if not is_integer(z0) and not is_float(z0):
            msg = 'Unable to set z0 coefficient for Sphere ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, z0)
            raise ValueError(msg)

        self._coeffs['z0'] = z0


    @r.setter
    def r(self, R):

        if not is_integer(R) and not is_float(R):
            msg = 'Unable to set R coefficient for Sphere ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, R)
            raise ValueError(msg)

        self._coeffs['R'] = R



class Cone(Surface):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):

        # Initialize Cone class attributes
        super(Cone, self).__init__(surface_id, boundary_type, name=name)

        self._coeff_keys = ['x0', 'y0', 'z0', 'R2']

        if not x0 is None:
            self.x0 = x0

        if not y0 is None:
            self.y0 = y0

        if not z0 is None:
            self.z0 = z0

        if not R2 is None:
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

        if not is_integer(x0) and not is_float(x0):
            msg = 'Unable to set x0 coefficient for Cone ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, x0)
            raise ValueError(msg)

        self._coeffs['x0'] = x0


    @y0.setter
    def y0(self, y0):

        if not is_integer(y0) and not is_float(y0):
            msg = 'Unable to set y0 coefficient for Cone ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, y0)
            raise ValueError(msg)

        self._coeffs['y0'] = y0


    @z0.setter
    def z0(self, z0):

        if not is_integer(z0) and not is_float(z0):
            msg = 'Unable to set z0 coefficient for Cone ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, z0)
            raise ValueError(msg)

        self._coeffs['z0'] = z0


    @r2.setter
    def r2(self, R2):

        if not is_integer(R2) and not is_float(R2):
            msg = 'Unable to set R^2 coefficient for Cone ID={0} to a ' \
                  'non-integer value {1}'.format(self._id, R2)
            raise ValueError(msg)

        self._coeffs['R2'] = R2



class XCone(Cone):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):

        # Initialize XCone class attributes
        super(XCone, self).__init__(surface_id, boundary_type, x0, y0,
                                    z0, R2, name=name)

        self._type = 'x-cone'



class YCone(Cone):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):

        # Initialize YCone class attributes
        super(YCone, self).__init__(surface_id, boundary_type, x0, y0, z0,
                                    R2, name=name)

        self._type = 'y-cone'



class ZCone(Cone):

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=None, y0=None, z0=None, R2=None, name=''):

        # Initialize ZCone class attributes
        super(ZCone, self).__init__(surface_id, boundary_type, x0, y0, z0,
                                    R2, name=name)

        self._type = 'z-cone'
