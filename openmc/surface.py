from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from numbers import Real, Integral
from xml.etree import ElementTree as ET
from warnings import warn

import numpy as np

from openmc.checkvalue import check_type, check_value
from openmc.region import Region, Intersection, Union
from openmc.mixin import IDManagerMixin


_BOUNDARY_TYPES = ['transmission', 'vacuum', 'reflective', 'periodic', 'white']

_WARNING_UPPER = """\
"{}(...) accepts an argument named '{}', not '{}'. Future versions of OpenMC \
will not accept the capitalized version.\
"""


class Surface(IDManagerMixin, metaclass=ABCMeta):
    """An implicit surface with an associated boundary condition.

    An implicit surface is defined as the set of zeros of a function of the
    three Cartesian coordinates. Surfaces in OpenMC are limited to a set of
    algebraic surfaces, i.e., surfaces that are polynomial in x, y, and z.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}
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

    next_id = 1
    used_ids = set()

    def __init__(self, surface_id=None, boundary_type='transmission', name=''):
        self.id = surface_id
        self.name = name
        self.boundary_type = boundary_type

        # A dictionary of the quadratic surface coefficients
        # Key      - coefficient name
        # Value    - coefficient value
        self._coefficients = {}

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

    @name.setter
    def name(self, name):
        if name is not None:
            check_type('surface name', name, str)
            self._name = name
        else:
            self._name = ''

    @boundary_type.setter
    def boundary_type(self, boundary_type):
        check_type('boundary type', boundary_type, str)
        check_value('boundary type', boundary_type, _BOUNDARY_TYPES)
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

    def clone(self, memo=None):
        """Create a copy of this surface with a new unique ID.

        Parameters
        ----------
        memo : dict or None
            A nested dictionary of previously cloned objects. This parameter
            is used internally and should not be specified by the user.

        Returns
        -------
        clone : openmc.Surface
            The clone of this surface

        """

        if memo is None:
            memo = {}

        # If no nemoize'd clone exists, instantiate one
        if self not in memo:
            clone = deepcopy(self)
            clone.id = None

            # Memoize the clone
            memo[self] = clone

        return memo[self]

    @abstractmethod
    def evaluate(self, point):
        pass

    @abstractmethod
    def translate(self, vector):
        pass

    def to_xml_element(self):
        """Return XML representation of the surface

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing source data

        """
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

    @staticmethod
    def from_xml_element(elem):
        """Generate surface from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.Surface
            Instance of a surface subclass

        """

        # Determine appropriate class
        surf_type = elem.get('type')
        surface_classes = {
            'plane': Plane,
            'x-plane': XPlane,
            'y-plane': YPlane,
            'z-plane': ZPlane,
            'x-cylinder': XCylinder,
            'y-cylinder': YCylinder,
            'z-cylinder': ZCylinder,
            'sphere': Sphere,
            'x-cone': XCone,
            'y-cone': YCone,
            'z-cone': ZCone,
            'quadric': Quadric,
        }
        cls = surface_classes[surf_type]

        # Determine ID, boundary type, coefficients
        kwargs = {}
        kwargs['surface_id'] = int(elem.get('id'))
        kwargs['boundary_type'] = elem.get('boundary', 'transmission')
        coeffs = [float(x) for x in elem.get('coeffs').split()]
        kwargs.update(dict(zip(cls._coeff_keys, coeffs)))

        return cls(**kwargs)

    @staticmethod
    def from_hdf5(group):
        """Create surface from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.Surface
            Instance of surface subclass

        """
        surface_id = int(group.name.split('/')[-1].lstrip('surface '))
        name = group['name'][()].decode() if 'name' in group else ''
        surf_type = group['type'][()].decode()
        bc = group['boundary_type'][()].decode()
        coeffs = group['coefficients'][...]

        # Create the Surface based on its type
        if surf_type == 'x-plane':
            x0 = coeffs[0]
            surface = XPlane(x0, bc, name, surface_id)

        elif surf_type == 'y-plane':
            y0 = coeffs[0]
            surface = YPlane(y0, bc, name, surface_id)

        elif surf_type == 'z-plane':
            z0 = coeffs[0]
            surface = ZPlane(z0, bc, name, surface_id)

        elif surf_type == 'plane':
            A, B, C, D = coeffs
            surface = Plane(A, B, C, D, bc, name, surface_id)

        elif surf_type == 'x-cylinder':
            y0, z0, r = coeffs
            surface = XCylinder(y0, z0, r, bc, name, surface_id)

        elif surf_type == 'y-cylinder':
            x0, z0, r = coeffs
            surface = YCylinder(x0, z0, r, bc, name, surface_id)

        elif surf_type == 'z-cylinder':
            x0, y0, r = coeffs
            surface = ZCylinder(x0, y0, r, bc, name, surface_id)

        elif surf_type == 'sphere':
            x0, y0, z0, r = coeffs
            surface = Sphere(x0, y0, z0, r, bc, name, surface_id)

        elif surf_type in ['x-cone', 'y-cone', 'z-cone']:
            x0, y0, z0, r2 = coeffs
            if surf_type == 'x-cone':
                surface = XCone(x0, y0, z0, r2, bc, name, surface_id)
            elif surf_type == 'y-cone':
                surface = YCone(x0, y0, z0, r2, bc, name, surface_id)
            elif surf_type == 'z-cone':
                surface = ZCone(x0, y0, z0, r2, bc, name, surface_id)

        elif surf_type == 'quadric':
            a, b, c, d, e, f, g, h, j, k = coeffs
            surface = Quadric(a, b, c, d, e, f, g, h, j, k, bc, name, surface_id)

        return surface


class Plane(Surface):
    """An arbitrary plane of the form :math:`Ax + By + Cz = D`.

    Parameters
    ----------
    a : float, optional
        The 'A' parameter for the plane. Defaults to 1.
    b : float, optional
        The 'B' parameter for the plane. Defaults to 0.
    c : float, optional
        The 'C' parameter for the plane. Defaults to 0.
    d : float, optional
        The 'D' parameter for the plane. Defaults to 0.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the plane. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

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
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'plane'
    _coeff_keys = ('a', 'b', 'c', 'd')

    def __init__(self, a=1., b=0., c=0., d=0., boundary_type='transmission',
                 name='', surface_id=None, **kwargs):
        super().__init__(surface_id, boundary_type, name=name)
        self._periodic_surface = None
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        for k, v in kwargs.items():
            if k in 'ABCD':
                warn(_WARNING_UPPER.format(type(self).__name__, k.lower(), k),
                    FutureWarning)
                setattr(self, k.lower(), v)

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
    def periodic_surface(self):
        return self._periodic_surface

    @a.setter
    def a(self, a):
        check_type('A coefficient', a, Real)
        self._coefficients['a'] = a

    @b.setter
    def b(self, b):
        check_type('B coefficient', b, Real)
        self._coefficients['b'] = b

    @c.setter
    def c(self, c):
        check_type('C coefficient', c, Real)
        self._coefficients['c'] = c

    @d.setter
    def d(self, d):
        check_type('D coefficient', d, Real)
        self._coefficients['d'] = d

    @periodic_surface.setter
    def periodic_surface(self, periodic_surface):
        check_type('periodic surface', periodic_surface, Plane)
        self._periodic_surface = periodic_surface
        periodic_surface._periodic_surface = self

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`Ax' + By' + Cz' - D`

        """

        x, y, z = point
        return self.a*x + self.b*y + self.c*z - self.d

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.Plane
            Translated surface

        """
        vx, vy, vz = vector
        d = self.d + self.a*vx + self.b*vy + self.c*vz
        if d == self.d:
            return self
        else:
            return type(self)(a=self.a, b=self.b, c=self.c, d=d)

    def to_xml_element(self):
        """Return XML representation of the surface

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing source data

        """
        element = super().to_xml_element()

        # Add periodic surface pair information
        if self.boundary_type == 'periodic':
            if self.periodic_surface is not None:
                element.set("periodic_surface_id", str(self.periodic_surface.id))
        return element

    @classmethod
    def from_points(cls, p1, p2, p3, **kwargs):
        """Return a plane given three points that pass through it.

        Parameters
        ----------
        p1, p2, p3 : 3-tuples
            Points that pass through the plane
        kwargs : dict
            Keyword arguments passed to the :class:`Plane` constructor

        Returns
        -------
        Plane
            Plane that passes through the three points

        """
        # Convert to numpy arrays
        p1 = np.asarray(p1)
        p2 = np.asarray(p2)
        p3 = np.asarray(p3)

        # Find normal vector to plane by taking cross product of two vectors
        # connecting p1->p2 and p1->p3
        n = np.cross(p2 - p1, p3 - p1)

        # The equation of the plane will by n·(<x,y,z> - p1) = 0. Determine
        # coefficients a, b, c, and d based on that
        a, b, c = n
        d = np.dot(n, p1)
        return cls(a=a, b=b, c=c, d=d, **kwargs)


class XPlane(Plane):
    """A plane perpendicular to the x axis of the form :math:`x - x_0 = 0`

    Parameters
    ----------
    x0 : float, optional
        Location of the plane. Defaults to 0.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface. Only axis-aligned periodicity is
        supported, i.e., x-planes can only be paired with x-planes.
    name : str, optional
        Name of the plane. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    x0 : float
        Location of the plane
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}
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

    _type = 'x-plane'
    _coeff_keys = ('x0',)

    def __init__(self, x0=0., boundary_type='transmission',
                 name='', surface_id=None):
        super().__init__(surface_id=surface_id, boundary_type=boundary_type, name=name)
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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`x' - x_0`

        """
        return point[0] - self.x0

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.XPlane
            Translated surface

        """
        vx = vector[0]
        if vx == 0:
            return self
        else:
            return type(self)(x0=self.x0 + vx)


class YPlane(Plane):
    """A plane perpendicular to the y axis of the form :math:`y - y_0 = 0`

    Parameters
    ----------
    y0 : float, optional
        Location of the plane
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface. Only axis-aligned periodicity is
        supported, i.e., x-planes can only be paired with x-planes.
    name : str, optional
        Name of the plane. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    y0 : float
        Location of the plane
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}
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

    _type = 'y-plane'
    _coeff_keys = ('y0',)

    def __init__(self, y0=0., boundary_type='transmission',
                 name='', surface_id=None):
        super().__init__(surface_id=surface_id, boundary_type=boundary_type, name=name)
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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`y' - y_0`

        """
        return point[1] - self.y0

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.YPlane
            Translated surface

        """
        vy = vector[1]
        if vy == 0.0:
            return self
        else:
            return type(self)(y0=self.y0 + vy)


class ZPlane(Plane):
    """A plane perpendicular to the z axis of the form :math:`z - z_0 = 0`

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}
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

    _type = 'z-plane'
    _coeff_keys = ('z0',)

    def __init__(self, z0=0., boundary_type='transmission',
                 name='', surface_id=None):
        super().__init__(surface_id=surface_id, boundary_type=boundary_type, name=name)
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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`z' - z_0`

        """
        return point[2] - self.z0

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.ZPlane
            Translated surface

        """
        vz = vector[2]
        if vz == 0.0:
            return self
        else:
            return type(self)(z0=self.z0 + vz)


class Cylinder(Surface):
    """A cylinder whose length is parallel to the x-, y-, or z-axis.

    Parameters
    ----------
    r : float, optional
        Radius of the cylinder. Defaults to 1.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the cylinder. If not specified, the name will be the empty
        string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    r : float
        Radius of the cylinder
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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
    def __init__(self, r=1., boundary_type='transmission',
                 name='', surface_id=None):
        super().__init__(surface_id, boundary_type, name=name)
        self.r = r

    @property
    def r(self):
        return self.coefficients['r']

    @r.setter
    def r(self, r):
        check_type('r coefficient', r, Real)
        self._coefficients['r'] = r


class XCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the x-axis of the form
    :math:`(y - y_0)^2 + (z - z_0)^2 = r^2`.

    Parameters
    ----------
    y0 : float, optional
        y-coordinate of the center of the cylinder. Defaults to 0.
    z0 : float, optional
        z-coordinate of the center of the cylinder. Defaults to 0.
    r : float, optional
        Radius of the cylinder. Defaults to 0.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the cylinder. If not specified, the name will be the empty
        string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    y0 : float
        y-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'x-cylinder'
    _coeff_keys = ('y0', 'z0', 'r')

    def __init__(self, y0=0., z0=0., r=1., boundary_type='transmission',
                 name='', surface_id=None, *, R=None):
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'), FutureWarning)
            r = R
        super().__init__(r, boundary_type, name, surface_id)
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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`(y' - y_0)^2 + (z' - z_0)^2 - r^2`

        """
        y = point[1] - self.y0
        z = point[2] - self.z0
        return y**2 + z**2 - self.r**2

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.XCylinder
            Translated surface

        """
        vx, vy, vz = vector
        if vy == 0.0 and vz == 0.0:
            return self
        else:
            y0 = self.y0 + vy
            z0 = self.z0 + vz
            return type(self)(y0=y0, z0=z0, r=self.r)


class YCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the y-axis of the form
    :math:`(x - x_0)^2 + (z - z_0)^2 = r^2`.

    Parameters
    ----------
    x0 : float, optional
        x-coordinate of the center of the cylinder. Defaults to 0.
    z0 : float, optional
        z-coordinate of the center of the cylinder. Defaults to 0.
    r : float, optional
        Radius of the cylinder. Defaults to 1.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the cylinder. If not specified, the name will be the empty
        string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    x0 : float
        x-coordinate of the center of the cylinder
    z0 : float
        z-coordinate of the center of the cylinder
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'y-cylinder'
    _coeff_keys = ('x0', 'z0', 'r')

    def __init__(self, x0=0., z0=0., r=1., boundary_type='transmission',
                 name='', surface_id=None, *, R=None):
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'), FutureWarning)
            r = R
        super().__init__(r, boundary_type, name, surface_id)
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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`(x' - x_0)^2 + (z' - z_0)^2 - r^2`

        """
        x = point[0] - self.x0
        z = point[2] - self.z0
        return x**2 + z**2 - self.r**2

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.YCylinder
            Translated surface

        """
        vx, vy, vz = vector
        if vx == 0.0 and vz == 0.0:
            return self
        else:
            x0 = self.x0 + vx
            z0 = self.z0 + vz
            return type(self)(x0=x0, z0=z0, r=self.r)


class ZCylinder(Cylinder):
    """An infinite cylinder whose length is parallel to the z-axis of the form
    :math:`(x - x_0)^2 + (y - y_0)^2 = r^2`.

    Parameters
    ----------
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    x0 : float, optional
        x-coordinate of the center of the cylinder. Defaults to 0.
    y0 : float, optional
        y-coordinate of the center of the cylinder. Defaults to 0.
    r : float, optional
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
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'z-cylinder'
    _coeff_keys = ('x0', 'y0', 'r')

    def __init__(self, x0=0., y0=0., r=1., boundary_type='transmission',
                 name='', surface_id=None, *, R=None):
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'), FutureWarning)
            r = R
        super().__init__(r, boundary_type, name, surface_id)
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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`(x' - x_0)^2 + (y' - y_0)^2 - r^2`

        """
        x = point[0] - self.x0
        y = point[1] - self.y0
        return x**2 + y**2 - self.r**2

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.ZCylinder
            Translated surface

        """
        vx, vy, vz = vector
        if vx == 0.0 and vy == 0.0:
            return self
        else:
            x0 = self.x0 + vx
            y0 = self.y0 + vy
            return type(self)(x0=x0, y0=y0, r=self.r)


class Sphere(Surface):
    """A sphere of the form :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 = r^2`.

    Parameters
    ----------
    x0 : float, optional
        x-coordinate of the center of the sphere. Defaults to 0.
    y0 : float, optional
        y-coordinate of the center of the sphere. Defaults to 0.
    z0 : float, optional
        z-coordinate of the center of the sphere. Defaults to 0.
    r : float, optional
        Radius of the sphere. Defaults to 1.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the sphere. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    x0 : float
        x-coordinate of the center of the sphere
    y0 : float
        y-coordinate of the center of the sphere
    z0 : float
        z-coordinate of the center of the sphere
    r : float
        Radius of the sphere
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'sphere'
    _coeff_keys = ('x0', 'y0', 'z0', 'r')

    def __init__(self, x0=0., y0=0., z0=0., r=1., boundary_type='transmission',
                 name='', surface_id=None, *, R=None):
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'), FutureWarning)
            r = R
        super().__init__(surface_id, boundary_type, name=name)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r = r

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
        return self.coefficients['r']

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
    def r(self, r):
        check_type('r coefficient', r, Real)
        self._coefficients['r'] = r

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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`(x' - x_0)^2 + (y' - y_0)^2 + (z' - z_0)^2 - r^2`

        """
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return x**2 + y**2 + z**2 - self.r**2

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.Sphere
            Translated surface

        """
        vx, vy, vz = vector
        if vx == 0.0 and vy == 0.0 and vz == 0.0:
            return self
        else:
            x0 = self.x0 + vx
            y0 = self.y0 + vy
            z0 = self.z0 + vz
            return type(self)(x0=x0, y0=y0, z0=z0, r=self.r)


class Cone(Surface):
    """A conical surface parallel to the x-, y-, or z-axis.

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
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
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
    r2 : float
        Parameter related to the aperature
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _coeff_keys = ('x0', 'y0', 'z0', 'r2')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., boundary_type='transmission',
                 name='', surface_id=None, *, R2=None):
        if R2 is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r2', 'R2'), FutureWarning)
            r2 = R2
        super().__init__(surface_id, boundary_type, name=name)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r2 = r2

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
    def r2(self, r2):
        check_type('r^2 coefficient', r2, Real)
        self._coefficients['r2'] = r2

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.Cone
            Translated surface

        """
        vx, vy, vz = vector
        if vx == 0.0 and vy == 0.0 and vz == 0.0:
            return self
        else:
            x0 = self.x0 + vx
            y0 = self.y0 + vy
            z0 = self.z0 + vz
            return type(self)(x0=x0, y0=y0, z0=z0, r2=self.r2)


class XCone(Cone):
    """A cone parallel to the x-axis of the form :math:`(y - y_0)^2 + (z - z_0)^2 =
    r^2 (x - x_0)^2`.

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
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the cone. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    r2 : float
        Parameter related to the aperature
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'x-cone'

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`(y' - y_0)^2 + (z' - z_0)^2 - r^2(x' - x_0)^2`

        """
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return y**2 + z**2 - self.r2*x**2


class YCone(Cone):
    """A cone parallel to the y-axis of the form :math:`(x - x_0)^2 + (z - z_0)^2 =
    r^2 (y - y_0)^2`.

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
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the cone. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    r2 : float
        Parameter related to the aperature
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'y-cone'

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`(x' - x_0)^2 + (z' - z_0)^2 - r^2(y' - y_0)^2`

        """
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return x**2 + z**2 - self.r2*y**2


class ZCone(Cone):
    """A cone parallel to the x-axis of the form :math:`(x - x_0)^2 + (y - y_0)^2 =
    r^2 (z - z_0)^2`.

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
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the cone. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    x0 : float
        x-coordinate of the apex
    y0 : float
        y-coordinate of the apex
    z0 : float
        z-coordinate of the apex
    r2 : float
        Parameter related to the aperature
    boundary_type : {'transmission, 'vacuum', 'reflective', 'white'}
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

    _type = 'z-cone'

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`(x' - x_0)^2 + (y' - y_0)^2 - r^2(z' - z_0)^2`

        """
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return x**2 + y**2 - self.r2*z**2


class Quadric(Surface):
    """A surface of the form :math:`Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy +
    Jz + K = 0`.

    Parameters
    ----------
    a, b, c, d, e, f, g, h, j, k : float, optional
        coefficients for the surface. All default to 0.
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}, optional
        Boundary condition that defines the behavior for particles hitting the
        surface. Defaults to transmissive boundary condition where particles
        freely pass through the surface.
    name : str, optional
        Name of the surface. If not specified, the name will be the empty string.
    surface_id : int, optional
        Unique identifier for the surface. If not specified, an identifier will
        automatically be assigned.

    Attributes
    ----------
    a, b, c, d, e, f, g, h, j, k : float
        coefficients for the surface
    boundary_type : {'transmission, 'vacuum', 'reflective', 'periodic', 'white'}
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

    _type = 'quadric'
    _coeff_keys = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'j', 'k')

    def __init__(self, a=0., b=0., c=0., d=0., e=0., f=0., g=0., h=0., j=0.,
                 k=0., boundary_type='transmission', name='', surface_id=None):
        super().__init__(surface_id, boundary_type, name=name)
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

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            The Cartesian coordinates, :math:`(x',y',z')`, at which the surface
            equation should be evaluated.

        Returns
        -------
        float
            :math:`Ax'^2 + By'^2 + Cz'^2 + Dx'y' + Ey'z' + Fx'z' + Gx' + Hy' +
            Jz' + K = 0`

        """
        x, y, z = point
        return x*(self.a*x + self.d*y + self.g) + \
            y*(self.b*y + self.e*z + self.h) + \
            z*(self.c*z + self.f*x + self.j) + self.k

    def translate(self, vector):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated

        Returns
        -------
        openmc.Quadric
            Translated surface

        """
        vx, vy, vz = vector
        a, b, c, d, e, f, g, h, j, k = (getattr(self, key) for key in
                                        self._coeff_keys)
        k = (k + vx*vx + vy*vy + vz*vz + d*vx*vy + e*vy*vz + f*vx*vz
             - g*vx - h*vy - j*vz)
        g = g - 2*a*vx - d*vy - f*vz
        h = h - 2*b*vy - d*vx - e*vz
        j = j - 2*c*vz - e*vy - f*vx
        return type(self)(a=a, b=b, c=c, d=d, e=e, f=f, g=g, h=h, j=j, k=k)


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

    >>> sphere = openmc.Sphere(surface_id=1, r=10.0)
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

    def __and__(self, other):
        if isinstance(other, Intersection):
            return Intersection([self] + other[:])
        else:
            return Intersection((self, other))

    def __or__(self, other):
        if isinstance(other, Union):
            return Union([self] + other[:])
        else:
            return Union((self, other))

    def __invert__(self):
        return -self.surface if self.side == '+' else +self.surface

    def __contains__(self, point):
        """Check whether a point is contained in the half-space.

        Parameters
        ----------
        point : 3-tuple of float
            Cartesian coordinates, :math:`(x',y',z')`, of the point

        Returns
        -------
        bool
            Whether the point is in the half-space

        """

        val = self.surface.evaluate(point)
        return val >= 0. if self.side == '+' else val < 0.

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

    def get_surfaces(self, surfaces=None):
        """
        Returns the surface that this is a halfspace of.

        Parameters
        ----------
        surfaces: collections.OrderedDict, optional
            Dictionary mapping surface IDs to :class:`openmc.Surface` instances

        Returns
        -------
        surfaces: collections.OrderedDict
            Dictionary mapping surface IDs to :class:`openmc.Surface` instances

        """
        if surfaces is None:
            surfaces = OrderedDict()

        surfaces[self.surface.id] = self.surface
        return surfaces

    def clone(self, memo=None):
        """Create a copy of this halfspace, with a cloned surface with a
        unique ID.

        Parameters
        ----------
        memo : dict or None
            A nested dictionary of previously cloned objects. This parameter
            is used internally and should not be specified by the user.

        Returns
        -------
        clone : openmc.Halfspace
            The clone of this halfspace

        """

        if memo is None:
            memo = dict

        clone = deepcopy(self)
        clone.surface = self.surface.clone(memo)
        return clone

    def translate(self, vector, memo=None):
        """Translate half-space in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which region should be translated
        memo : dict or None
            Dictionary used for memoization

        Returns
        -------
        openmc.Halfspace
            Translated half-space

        """
        if memo is None:
            memo = {}

        # If translated surface not in memo, add it
        key = (self.surface, tuple(vector))
        if key not in memo:
            memo[key] = self.surface.translate(vector)

        # Return translated surface
        return type(self)(memo[key], self.side)
