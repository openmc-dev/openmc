from abc import ABCMeta
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

import numpy as np

from openmc.checkvalue import check_type, check_value, check_greater_than
from openmc.region import Region

if sys.version_info[0] >= 3:
    basestring = str

# A static variable for auto-generated Surface IDs
AUTO_SURFACE_ID = 10000

_BC_TYPES = ['transmission', 'vacuum', 'reflective', 'periodic']


def reset_auto_surface_id():
    global AUTO_SURFACE_ID
    AUTO_SURFACE_ID = 10000


class Surface(object):
    """An implicit surface with an associated boundary condition.

    An implicit surface is defined as the set of zeros of a function of the
    three Cartesian coordinates. Surfaces in OpenMC are limited to a set of
    algebraic surfaces, i.e., surfaces that are polynomial in x, y, and z.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface. Note that periodic boundary conditions
        can only be applied to x-, y-, and z-planes, and only axis-aligned
        periodicity is supported.
    name : str, optional
        Name of the surface. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission', name=''):
        self.id = surface_id
        self.name = name
        self._type = ''
        self.boundary_type = boundary_type

        # A dictionary of the quadratic surface coefficients
        # Key        - coefficeint name
        # Value    - coefficient value
        self._coefficients = {}

        # An ordered list of the coefficient names to export to XML in the
        # proper order
        self._coeff_keys = []

    def __neg__(self):
        return Halfspace(self, '-')

    def __pos__(self):
        return Halfspace(self, '+')

    def __repr__(self):
        string = 'Surface\n'
        string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
        string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
        string += '{0: <16}{1}{2}\n'.format('\tBoundary', '=\t', self._boundary_type)

        coefficients = '{0: <16}'.format('\tCoefficients') + '\n'

        for coeff in self._coefficients:
            coefficients += '{0: <16}{1}{2}\n'.format(
                coeff, '=\t', self._coefficients[coeff])

        string += coefficients

        return string

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
    def coefficients(self):
        return self._coefficients

    @id.setter
    def id(self, surface_id):
        if surface_id is None:
            global AUTO_SURFACE_ID
            self._id = AUTO_SURFACE_ID
            AUTO_SURFACE_ID += 1
        else:
            check_type('surface ID', surface_id, Integral)
            check_greater_than('surface ID', surface_id, 0, equality=True)
            self._id = surface_id

    @name.setter
    def name(self, name):
        if name is not None:
            check_type('surface name', name, basestring)
            self._name = name
        else:
            self._name = ''

    @boundary_type.setter
    def boundary_type(self, boundary_type):
        check_type('boundary type', boundary_type, basestring)
        check_value('boundary type', boundary_type, _BC_TYPES)
        self._boundary_type = boundary_type

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. If the half-space is
        unbounded in a particular direction, numpy.inf is used to represent
        infinity.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        return (np.array([-np.inf, -np.inf, -np.inf]),
                np.array([np.inf, np.inf, np.inf]))

    def create_xml_subelement(self):
        element = ET.Element("surface")
        element.set("id", str(self._id))

        if len(self._name) > 0:
            element.set("name", str(self._name))

        element.set("type", self._type)
        if self.boundary_type != 'transmission':
            element.set("boundary", self.boundary_type)
        element.set("coeffs", ' '.join([str(self._coefficients.setdefault(key, 0.0))
                                        for key in self._coeff_keys]))

        return element


class Plane(Surface):
    """An arbitrary plane of the form :math:`Ax + By + Cz = D`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    A : float, optional
        The 'A' parameter for the plane. Defaults to 1.
    B : float, optional
        The 'B' parameter for the plane. Defaults to 0.
    C : float, optional
        The 'C' parameter for the plane. Defaults to 0.
    D : float, optional
        The 'D' parameter for the plane. Defaults to 0.
    name : str, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    periodic_surface : openmc.Surface
        If a periodic boundary condition is used, the surface with which this
        one is periodic with
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 A=1., B=0., C=0., D=0., name=''):
        super(Plane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'plane'
        self._coeff_keys = ['A', 'B', 'C', 'D']
        self._periodic_surface = None
        self.a = A
        self.b = B
        self.c = C
        self.d = D

    @property
    def a(self):
        return self.coefficients['A']

    @property
    def b(self):
        return self.coefficients['B']

    @property
    def c(self):
        return self.coefficients['C']

    @property
    def d(self):
        return self.coefficients['D']

    @property
    def periodic_surface(self):
        return self._periodic_surface

    @a.setter
    def a(self, A):
        check_type('A coefficient', A, Real)
        self._coefficients['A'] = A

    @b.setter
    def b(self, B):
        check_type('B coefficient', B, Real)
        self._coefficients['B'] = B

    @c.setter
    def c(self, C):
        check_type('C coefficient', C, Real)
        self._coefficients['C'] = C

    @d.setter
    def d(self, D):
        check_type('D coefficient', D, Real)
        self._coefficients['D'] = D

    @periodic_surface.setter
    def periodic_surface(self, periodic_surface):
        check_type('periodic surface', periodic_surface, Plane)
        self._periodic_surface = periodic_surface
        periodic_surface._periodic_surface = self

    def create_xml_subelement(self):
        element = super(Plane, self).create_xml_subelement()

        # Add periodic surface pair information
        if self.boundary_type == 'periodic':
            if self.periodic_surface is not None:
                element.set("periodic_surface_id", str(self.periodic_surface.id))
        return element


class XPlane(Plane):
    """A plane perpendicular to the x axis of the form :math:`x - x_0 = 0`

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface. Only axis-aligned periodicity is
        supported, i.e., x-planes can only be paired with x-planes.
    x0 : float, optional
        Location of the plane. Defaults to 0.
    name : str, optional
        Name of the plane. If not specified, the name will be the empty string.

    Attributes
    ----------
    x0 : float
        Location of the plane
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    periodic_surface : openmc.Surface
        If a periodic boundary condition is used, the surface with which this
        one is periodic with
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., name=''):
        super(XPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'x-plane'
        self._coeff_keys = ['x0']
        self.x0 = x0

    @property
    def x0(self):
        return self.coefficients['x0']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coefficients['x0'] = x0

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. For the x-plane surface, the
        half-spaces are unbounded in their y- and z- directions. To represent
        infinity, numpy.inf is used.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        if side == '-':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([self.x0, np.inf, np.inf]))
        elif side == '+':
            return (np.array([self.x0, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))


class YPlane(Plane):
    """A plane perpendicular to the y axis of the form :math:`y - y_0 = 0`

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface. Only axis-aligned periodicity is
        supported, i.e., x-planes can only be paired with x-planes.
    y0 : float, optional
        Location of the plane
    name : str, optional
        Name of the plane. If not specified, the name will be the empty string.

    Attributes
    ----------
    y0 : float
        Location of the plane
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    periodic_surface : openmc.Surface
        If a periodic boundary condition is used, the surface with which this
        one is periodic with
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 y0=0., name=''):
        # Initialize YPlane class attributes
        super(YPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'y-plane'
        self._coeff_keys = ['y0']
        self.y0 = y0

    @property
    def y0(self):
        return self.coefficients['y0']

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coefficients['y0'] = y0

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. For the y-plane surface, the
        half-spaces are unbounded in their x- and z- directions. To represent
        infinity, numpy.inf is used.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        if side == '-':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, self.y0, np.inf]))
        elif side == '+':
            return (np.array([-np.inf, self.y0, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))


class ZPlane(Plane):
    """A plane perpendicular to the z axis of the form :math:`z - z_0 = 0`

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface. Only axis-aligned periodicity is
        supported, i.e., x-planes can only be paired with x-planes.
    z0 : float, optional
        Location of the plane. Defaults to 0.
    name : str, optional
        Name of the plane. If not specified, the name will be the empty string.

    Attributes
    ----------
    z0 : float
        Location of the plane
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    periodic_surface : openmc.Surface
        If a periodic boundary condition is used, the surface with which this
        one is periodic with
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 z0=0., name=''):
        # Initialize ZPlane class attributes
        super(ZPlane, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'z-plane'
        self._coeff_keys = ['z0']
        self.z0 = z0

    @property
    def z0(self):
        return self.coefficients['z0']

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coefficients['z0'] = z0

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. For the z-plane surface, the
        half-spaces are unbounded in their x- and y- directions. To represent
        infinity, numpy.inf is used.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        if side == '-':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, self.z0]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, self.z0]),
                    np.array([np.inf, np.inf, np.inf]))


class Cylinder(Surface):
    """A cylinder whose length is parallel to the x-, y-, or z-axis.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    R : float, optional
        Radius of the cylinder. Defaults to 1.
    name : str, optional
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    r : float
        Radius of the cylinder
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    __metaclass__ = ABCMeta

    def __init__(self, surface_id=None, boundary_type='transmission',
                 R=1., name=''):
        super(Cylinder, self).__init__(surface_id, boundary_type, name=name)

        self._coeff_keys = ['R']
        self.r = R

    @property
    def r(self):
        return self.coefficients['R']

    @r.setter
    def r(self, R):
        check_type('R coefficient', R, Real)
        self._coefficients['R'] = R


class XCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the x-axis of the form
    :math:`(y - y_0)^2 + (z - z_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    y0 : float, optional
        y-coordinate of the center of the cylinder. Defaults to 0.
    z0 : float, optional
        z-coordinate of the center of the cylinder. Defaults to 0.
    R : float, optional
        Radius of the cylinder. Defaults to 0.
    name : str, optional
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    y0 : float
        y-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 y0=0., z0=0., R=1., name=''):
        super(XCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'x-cylinder'
        self._coeff_keys = ['y0', 'z0', 'R']
        self.y0 = y0
        self.z0 = z0

    @property
    def y0(self):
        return self.coefficients['y0']

    @property
    def z0(self):
        return self.coefficients['z0']

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coefficients['y0'] = y0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coefficients['z0'] = z0

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. For the x-cylinder surface,
        the negative half-space is unbounded in the x- direction and the
        positive half-space is unbounded in all directions. To represent
        infinity, numpy.inf is used.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        if side == '-':
            return (np.array([-np.inf, self.y0 - self.r, self.z0 - self.r]),
                    np.array([np.inf, self.y0 + self.r, self.z0 + self.r]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))


class YCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the y-axis of the form
    :math:`(x - x_0)^2 + (z - z_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the center of the cylinder. Defaults to 0.
    z0 : float, optional
        z-coordinate of the center of the cylinder. Defaults to 0.
    R : float, optional
        Radius of the cylinder. Defaults to 1.
    name : str, optional
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., z0=0., R=1., name=''):
        super(YCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'y-cylinder'
        self._coeff_keys = ['x0', 'z0', 'R']
        self.x0 = x0
        self.z0 = z0

    @property
    def x0(self):
        return self.coefficients['x0']

    @property
    def z0(self):
        return self.coefficients['z0']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coefficients['x0'] = x0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coefficients['z0'] = z0

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. For the y-cylinder surface,
        the negative half-space is unbounded in the y- direction and the
        positive half-space is unbounded in all directions. To represent
        infinity, numpy.inf is used.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        if side == '-':
            return (np.array([self.x0 - self.r, -np.inf, self.z0 - self.r]),
                    np.array([self.x0 + self.r, np.inf, self.z0 + self.r]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))


class ZCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the z-axis of the form
    :math:`(x - x_0)^2 + (y - y_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the center of the cylinder. Defaults to 0.
    y0 : float, optional
        y-coordinate of the center of the cylinder. Defaults to 0.
    R : float, optional
        Radius of the cylinder. Defaults to 1.
    name : str, optional
        Name of the cylinder. If not specified, the name will be the empty
        string.

    Attributes
    ----------
    x0 : float
        x-coordinate of the center of the cylinder
    y0 : float
        y-coordinate of the center of the cylinder
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., y0=0., R=1., name=''):
        super(ZCylinder, self).__init__(surface_id, boundary_type, R, name=name)

        self._type = 'z-cylinder'
        self._coeff_keys = ['x0', 'y0', 'R']
        self.x0 = x0
        self.y0 = y0

    @property
    def x0(self):
        return self.coefficients['x0']

    @property
    def y0(self):
        return self.coefficients['y0']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coefficients['x0'] = x0

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coefficients['y0'] = y0

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. For the z-cylinder surface,
        the negative half-space is unbounded in the z- direction and the
        positive half-space is unbounded in all directions. To represent
        infinity, numpy.inf is used.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        if side == '-':
            return (np.array([self.x0 - self.r, self.y0 - self.r, -np.inf]),
                    np.array([self.x0 + self.r, self.y0 + self.r, np.inf]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))


class Sphere(Surface):
    """A sphere of the form :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 = R^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the center of the sphere. Defaults to 0.
    y0 : float, optional
        y-coordinate of the center of the sphere. Defaults to 0.
    z0 : float, optional
        z-coordinate of the center of the sphere. Defaults to 0.
    R : float, optional
        Radius of the sphere. Defaults to 1.
    name : str, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., y0=0., z0=0., R=1., name=''):
        super(Sphere, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'sphere'
        self._coeff_keys = ['x0', 'y0', 'z0', 'R']
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r = R

    @property
    def x0(self):
        return self.coefficients['x0']

    @property
    def y0(self):
        return self.coefficients['y0']

    @property
    def z0(self):
        return self.coefficients['z0']

    @property
    def r(self):
        return self.coefficients['R']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coefficients['x0'] = x0

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coefficients['y0'] = y0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coefficients['z0'] = z0

    @r.setter
    def r(self, R):
        check_type('R coefficient', R, Real)
        self._coefficients['R'] = R

    def bounding_box(self, side):
        """Determine an axis-aligned bounding box.

        An axis-aligned bounding box for surface half-spaces is represented by
        its lower-left and upper-right coordinates. The positive half-space of a
        sphere is unbounded in all directions. To represent infinity, numpy.inf
        is used.

        Parameters
        ----------
        side : {'+', '-'}
            Indicates the negative or positive half-space

        Returns
        -------
        numpy.ndarray
            Lower-left coordinates of the axis-aligned bounding box for the
            desired half-space
        numpy.ndarray
            Upper-right coordinates of the axis-aligned bounding box for the
            desired half-space

        """

        if side == '-':
            return (np.array([self.x0 - self.r, self.y0 - self.r,
                              self.z0 - self.r]),
                    np.array([self.x0 + self.r, self.y0 + self.r,
                              self.z0 + self.r]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))


class Cone(Surface):
    """A conical surface parallel to the x-, y-, or z-axis.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the apex. Defaults to 0.
    y0 : float
        y-coordinate of the apex. Defaults to 0.
    z0 : float
        z-coordinate of the apex. Defaults to 0.
    R2 : float
        Parameter related to the aperature. Defaults to 1.
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
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    __metaclass__ = ABCMeta

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., y0=0., z0=0., R2=1., name=''):
        super(Cone, self).__init__(surface_id, boundary_type, name=name)

        self._coeff_keys = ['x0', 'y0', 'z0', 'R2']
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r2 = R2

    @property
    def x0(self):
        return self.coefficients['x0']

    @property
    def y0(self):
        return self.coefficients['y0']

    @property
    def z0(self):
        return self.coefficients['z0']

    @property
    def r2(self):
        return self.coefficients['r2']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coefficients['x0'] = x0

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coefficients['y0'] = y0

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coefficients['z0'] = z0

    @r2.setter
    def r2(self, R2):
        check_type('R^2 coefficient', R2, Real)
        self._coefficients['R2'] = R2


class XCone(Cone):
    """A cone parallel to the x-axis of the form :math:`(y - y_0)^2 + (z - z_0)^2 =
    R^2 (x - x_0)^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the apex. Defaults to 0.
    y0 : float, optional
        y-coordinate of the apex. Defaults to 0.
    z0 : float, optional
        z-coordinate of the apex. Defaults to 0.
    R2 : float, optional
        Parameter related to the aperature. Defaults to 1.
    name : str, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., y0=0., z0=0., R2=1., name=''):
        super(XCone, self).__init__(surface_id, boundary_type, x0, y0,
                                    z0, R2, name=name)

        self._type = 'x-cone'


class YCone(Cone):
    """A cone parallel to the y-axis of the form :math:`(x - x_0)^2 + (z - z_0)^2 =
    R^2 (y - y_0)^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the apex. Defaults to 0.
    y0 : float, optional
        y-coordinate of the apex. Defaults to 0.
    z0 : float, optional
        z-coordinate of the apex. Defaults to 0.
    R2 : float, optional
        Parameter related to the aperature. Defaults to 1.
    name : str, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., y0=0., z0=0., R2=1., name=''):
        super(YCone, self).__init__(surface_id, boundary_type, x0, y0, z0,
                                    R2, name=name)

        self._type = 'y-cone'


class ZCone(Cone):
    """A cone parallel to the x-axis of the form :math:`(x - x_0)^2 + (y - y_0)^2 =
    R^2 (z - z_0)^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the apex. Defaults to 0.
    y0 : float, optional
        y-coordinate of the apex. Defaults to 0.
    z0 : float, optional
        z-coordinate of the apex. Defaults to 0.
    R2 : float, optional
        Parameter related to the aperature. Defaults to 1.
    name : str, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 x0=0., y0=0., z0=0., R2=1., name=''):
        super(ZCone, self).__init__(surface_id, boundary_type, x0, y0, z0,
                                    R2, name=name)

        self._type = 'z-cone'


class Quadric(Surface):
    """A surface of the form :math:`Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy +
    Jz + K = 0`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    a, b, c, d, e, f, g, h, j, k : float, optional
        coefficients for the surface. All default to 0.
    name : str, optional
        Name of the sphere. If not specified, the name will be the empty string.

    Attributes
    ----------
    a, b, c, d, e, f, g, h, j, k : float
        coefficients for the surface
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic'}
        Boundary condition that defines the behavior for particles hitting the
        surface.
    coefficients : dict
        Dictionary of surface coefficients
    id : int
        Unique identifier for the surface
    name : str
        Name of the surface
    type : str
        Type of the surface

    """

    def __init__(self, surface_id=None, boundary_type='transmission',
                 a=0., b=0., c=0., d=0., e=0., f=0., g=0.,
                 h=0., j=0., k=0., name=''):
        super(Quadric, self).__init__(surface_id, boundary_type, name=name)

        self._type = 'quadric'
        self._coeff_keys = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'j', 'k']
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.g = g
        self.h = h
        self.j = j
        self.k = k

    @property
    def a(self):
        return self.coefficients['a']

    @property
    def b(self):
        return self.coefficients['b']

    @property
    def c(self):
        return self.coefficients['c']

    @property
    def d(self):
        return self.coefficients['d']

    @property
    def e(self):
        return self.coefficients['e']

    @property
    def f(self):
        return self.coefficients['f']

    @property
    def g(self):
        return self.coefficients['g']

    @property
    def h(self):
        return self.coefficients['h']

    @property
    def j(self):
        return self.coefficients['j']

    @property
    def k(self):
        return self.coefficients['k']

    @a.setter
    def a(self, a):
        check_type('a coefficient', a, Real)
        self._coefficients['a'] = a

    @b.setter
    def b(self, b):
        check_type('b coefficient', b, Real)
        self._coefficients['b'] = b

    @c.setter
    def c(self, c):
        check_type('c coefficient', c, Real)
        self._coefficients['c'] = c

    @d.setter
    def d(self, d):
        check_type('d coefficient', d, Real)
        self._coefficients['d'] = d

    @e.setter
    def e(self, e):
        check_type('e coefficient', e, Real)
        self._coefficients['e'] = e

    @f.setter
    def f(self, f):
        check_type('f coefficient', f, Real)
        self._coefficients['f'] = f

    @g.setter
    def g(self, g):
        check_type('g coefficient', g, Real)
        self._coefficients['g'] = g

    @h.setter
    def h(self, h):
        check_type('h coefficient', h, Real)
        self._coefficients['h'] = h

    @j.setter
    def j(self, j):
        check_type('j coefficient', j, Real)
        self._coefficients['j'] = j

    @k.setter
    def k(self, k):
        check_type('k coefficient', k, Real)
        self._coefficients['k'] = k


class Halfspace(Region):
    """A positive or negative half-space region.

    A half-space is either of the two parts into which a two-dimension surface
    divides the three-dimensional Euclidean space. If the equation of the
    surface is :math:`f(x,y,z) = 0`, the region for which :math:`f(x,y,z) < 0`
    is referred to as the negative half-space and the region for which
    :math:`f(x,y,z) > 0` is referred to as the positive half-space.

    Instances of Halfspace are generally not instantiated directly. Rather, they
    can be created from an existing Surface through the __neg__ and __pos__
    operators, as the following example demonstrates:

    >>> sphere = openmc.Sphere(surface_id=1, R=10.0)
    >>> inside_sphere = -sphere
    >>> outside_sphere = +sphere
    >>> type(inside_sphere)
    <class 'openmc.surface.Halfspace'>

    Parameters
    ----------
    surface : openmc.Surface
        Surface which divides Euclidean space.
    side : {'+', '-'}
        Indicates whether the positive or negative half-space is used.

    Attributes
    ----------
    surface : openmc.Surface
        Surface which divides Euclidean space.
    side : {'+', '-'}
        Indicates whether the positive or negative half-space is used.
    bounding_box : tuple of numpy.ndarray
        Lower-left and upper-right coordinates of an axis-aligned bounding box

    """

    def __init__(self, surface, side):
        self.surface = surface
        self.side = side

    def __invert__(self):
        return -self.surface if self.side == '+' else +self.surface

    @property
    def surface(self):
        return self._surface

    @surface.setter
    def surface(self, surface):
        check_type('surface', surface, Surface)
        self._surface = surface

    @property
    def side(self):
        return self._side

    @side.setter
    def side(self, side):
        check_value('side', side, ('+', '-'))
        self._side = side

    @property
    def bounding_box(self):
        return self.surface.bounding_box(self.side)

    def __str__(self):
        return '-' + str(self.surface.id) if self.side == '-' \
            else str(self.surface.id)
