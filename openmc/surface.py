from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from collections.abc import Iterable
from copy import deepcopy
from numbers import Real, Integral
from xml.etree import ElementTree as ET
from warnings import warn

import numpy as np

from openmc.checkvalue import check_type, check_value, check_length
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

    @classmethod
    def get_subclasses(cls):
        """Recursively find all subclasses of this class"""
        return set(cls.__subclasses__()).union([s for c in cls.__subclasses__()
                                                for s in get_subclasses(c)])

    @classmethod
    def get_subclass_map(cls):
        """Generate mapping of class _type attributes to classes"""
        return {c._type: c for c in cls.get_subclasses()}

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
        cls = _SURFACE_CLASSES[surf_type]

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
        kwargs = {'boundary_type': bc, 'name': name, 'surface_id': surface_id}

        cls = _SURFACE_CLASSES[surf_type]

        return cls(*coeffs, **kwargs)


_SURFACE_CLASSES = Surface.get_subclass_map()


class PlaneMeta(metaclass=ABCMeta):
    """A Plane Meta class for all operations on order 1 surfaces"""
#    def __new__(cls, *args, **kwargs):
#       for key in cls._coeff_keys, kwargs.pop(key) to get arguments for class
#        """Simplify this plane if possible to an XPlane, YPlane, or ZPlane"""
#        pass
#
#        if cls is Plane:
#            if np.all(np.isclose((self.b, self.c), 0., atol=atol)):
#                x0 = self.d / self.a
#                return XPlane(x0=x0, boundary_type=self.boundary_type,
#                              name=self.name, surface_id=self.id)
#
#            elif np.all(np.isclose((self.a, self.c), 0., atol=atol)):
#                y0 = self.d / self.b
#                return YPlane(y0=y0, boundary_type=self.boundary_type,
#                              name=self.name, surface_id=self.id)
#
#            elif np.all(np.isclose((self.a, self.b), 0., atol=atol)):
#                z0 = self.d / self.c
#                return ZPlane(z0=z0, boundary_type=self.boundary_type,
#                              name=self.name, surface_id=self.id)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._periodic_surface = None

    @property
    def periodic_surface(self):
        return self._periodic_surface

    @periodic_surface.setter
    def periodic_surface(self, periodic_surface):
        check_type('periodic surface', periodic_surface, Plane)
        self._periodic_surface = periodic_surface
        periodic_surface._periodic_surface = self

    @abstractmethod
    def _get_base_coeffs(self):
        pass

    @abstractmethod
    def _update_from_base_coeffs(self):
        pass

    def _neg_bounds(self):
        return (np.array([-np.inf, -np.inf, -np.inf]),
                np.array([np.inf, np.inf, np.inf]))

    def _pos_bounds(self):
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
            :math:`Ax' + By' + Cz' - D`

        """

        x, y, z = point
        a, b, c, d = self._get_base_coeffs()
        return a*x + b*y + c*z - d

    def translate(self, vector, clone=False):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated
        clone : boolean
            Whether or not to return a new instance of a Plane or to modify the
            coefficients of this plane.

        Returns
        -------
        openmc.Plane
            Translated surface

        """
        vx, vy, vz = vector
        a, b, c, d = self._get_base_coeffs()
        d = d + a*vx + b*vy + c*vz

        if clone:
            surf = self.clone()
        else:
            surf = self

        surf._update_from_base_coeffs(a, b, c, d)
        return surf

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
            return self._neg_bounds()
        elif side == '+':
            return self._pos_bounds()

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
                element.set("periodic_surface_id",
                            str(self.periodic_surface.id))
        return element


class Plane(PlaneMeta, Surface):
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

    def __init__(self, a=1., b=0., c=0., d=0., **kwargs):
        # work around until capital letter kwargs are deprecated
        oldkwargs = deepcopy(kwargs)
        for k in 'ABCD':
            kwargs.pop(k, None)

        super().__init__(**kwargs)
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        for k, v in oldkwargs.items():
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

    def _get_base_coeffs(self):
        return (self.a, self.b, self.c, self.d)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d = coeffs
        self.a = a
        self.b = b
        self.c = c
        self.d = d

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


class XPlane(PlaneMeta, Surface):
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

    def __init__(self, x0=0., **kwargs):
        super().__init__(**kwargs)
        self.x0 = x0

    @property
    def x0(self):
        return self.coefficients['x0']

    @x0.setter
    def x0(self, x0):
        check_type('x0 coefficient', x0, Real)
        self._coefficients['x0'] = x0

    def _get_base_coeffs(self):
        return (1., 0., 0., self.x0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d = coeffs
        self.x0 = d / a

    def _neg_bounds(self):
        """Return the lower and upper bounds of the negative half space"""
        return (np.array([-np.inf, -np.inf, -np.inf]),
                np.array([self.x0, np.inf, np.inf]))

    def _pos_bounds(self):
        """Return the lower and upper bounds of the positive half space"""
        return (np.array([self.x0, -np.inf, -np.inf]),
                np.array([np.inf, np.inf, np.inf]))


class YPlane(PlaneMeta, Surface):
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

    def __init__(self, y0=0., **kwargs):
        super().__init__(**kwargs)
        self.y0 = y0

    @property
    def y0(self):
        return self.coefficients['y0']

    @y0.setter
    def y0(self, y0):
        check_type('y0 coefficient', y0, Real)
        self._coefficients['y0'] = y0

    def _get_base_coeffs(self):
        """Return generalized coefficients for a plane"""
        return (0., 1., 0., self.y0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d = coeffs
        self.y0 = d / b

    def _neg_bounds(self):
        """Return the lower and upper bounds of the negative half space"""
        return (np.array([-np.inf, -np.inf, -np.inf]),
                np.array([np.inf, self.y0, np.inf]))

    def _pos_bounds(self):
        """Return the lower and upper bounds of the positive half space"""
        return (np.array([-np.inf, self.y0, -np.inf]),
                np.array([np.inf, np.inf, np.inf]))


class ZPlane(PlaneMeta, Surface):
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

    def __init__(self, z0=0., **kwargs):
        super().__init__(**kwargs)
        self.z0 = z0

    @property
    def z0(self):
        return self.coefficients['z0']

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coefficients['z0'] = z0

    def _get_base_coeffs(self):
        """Return generalized coefficients for a plane"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d = coeffs
        self.z0 = d / c

    def _neg_bounds(self):
        """Return the lower and upper bounds of the negative half space"""
        return (np.array([-np.inf, -np.inf, -np.inf]),
                np.array([np.inf, np.inf, self.z0]))

    def _pos_bounds(self):
        """Return the lower and upper bounds of the positive half space"""
        return (np.array([-np.inf, -np.inf, self.z0]),
                np.array([np.inf, np.inf, np.inf]))


class QuadricMeta(metaclass=ABCMeta):
    """A Meta class implementing common functionality for quadric surfaces"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @abstractmethod
    def _get_base_coeffs(self):
        pass

    @abstractmethod
    def _update_from_base_coeffs(self):
        pass

    def _neg_bounds(self):
        return (np.array([-np.inf, -np.inf, -np.inf]),
                np.array([np.inf, np.inf, np.inf]))

    def _pos_bounds(self):
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
            :math:`Ax'^2 + By'^2 + Cz'^2 + Dx'y' + Ey'z' + Fx'z' + Gx' + Hy' +
            Jz' + K = 0`

        """
        x, y, z = point
        a, b, c, d, e, f, g, h, j, k = self._get_base_coeffs()
        return x*(a*x + d*y + g) + y*(b*y + e*z + h) + z*(c*z + f*x + j) + k

    def translate(self, vector, clone=False):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated
        clone : bool
            Whether to return a clone of the Surface or the Surface itself

        Returns
        -------
        openmc.QuadricMeta
            Translated surface

        """
        vx, vy, vz = vector
        a, b, c, d, e, f, g, h, j, k = self._get_base_coeffs()
        k = (k + vx*vx + vy*vy + vz*vz + d*vx*vy + e*vy*vz + f*vx*vz
             - g*vx - h*vy - j*vz)
        g = g - 2*a*vx - d*vy - f*vz
        h = h - 2*b*vy - d*vx - e*vz
        j = j - 2*c*vz - e*vy - f*vx

        if clone:
            surf = self.clone()
        else:
            surf = self

        surf._update_from_base_coeffs((a, b, c, d, e, f, g, h, j, k))

        return surf


class Cylinder(QuadricMeta, Surface):
    """A cylinder with radius r, centered on the point (x0, y0, z0) with an
    axis specified by the line through points (x0, y0, z0) and (x0+u, y0+v,
    z0+w)

    Parameters
    ----------
    x0 : float, optional
        x-coordinate for the origin of the Cylinder. Defaults to 0
    y0 : float, optional
        y-coordinate for the origin of the Cylinder. Defaults to 0
    z0 : float, optional
        z-coordinate for the origin of the Cylinder. Defaults to 0
    r : float, optional
        Radius of the cylinder. Defaults to 1.
    u : float, optional
        x-component of the vector representing the axis of the cylinder.
        Defaults to 0.
    v : float, optional
        y-component of the vector representing the axis of the cylinder.
        Defaults to 0.
    w : float, optional
        z-component of the vector representing the axis of the cylinder.
        Defaults to 1.
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
        x-coordinate for the origin of the Cylinder
    y0 : float
        y-coordinate for the origin of the Cylinder
    z0 : float
        z-coordinate for the origin of the Cylinder
    r : float
        Radius of the cylinder
    u : float
        x-component of the vector representing the axis of the cylinder
    v : float
        y-component of the vector representing the axis of the cylinder
    w : float
        z-component of the vector representing the axis of the cylinder
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
    _type = 'cylinder'
    _coeff_keys = ('x0', 'y0', 'z0', 'r', 'u', 'v','w')
    def __init__(self, x0=0., y0=0., z0=0. r=1., u=0., v=0., w=1., **kwargs):
        super().__init__(**kwargs)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r = r
        self.u = u
        self.v = v
        self.w = w

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

    @property
    def u(self):
        return self.coefficients['u']

    @property
    def v(self):
        return self.coefficients['v']

    @property
    def w(self):
        return self.coefficients['w']

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

    @u.setter
    def u(self, u):
        check_type('u coefficient', u, Real)
        self._coefficients['u'] = u

    @v.setter
    def v(self, v):
        check_type('v coefficient', v, Real)
        self._coefficients['v'] = v

    @w.setter
    def w(self, w):
        check_type('w coefficient', w, Real)
        self._coefficients['w'] = w

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        x0, y0, z0 = self.x0, self.y0, self.z0
        x1, y1, z1 = self.x0 + self.u, self.y0 + self.v, self.z0 + self.w
        dx, dy, dz = x1 - x0, y1 - y0, z1 - z0

        # Set coefficients for Quadric surface that represents a cylinder of
        # radius r whose axis is the line defined by p1 and p2
        a = dy**2 + dz**2
        b = dx**2 + dz**2
        c = dx**2 + dy**2
        d = -2*dx*dy
        e = -2*dy*dz
        f = -2*dx*dz
        g = -2*((z1 - z0)*(x0*z1 - x1*z0) + (y1 - y0)*(x0*y1-x1*y0))
        h = 2*((x1 - x0)*(x0*y1 - x1*y0) - (z1 - z0)*(y0*z1 - y1*z0))
        j = 2*((x1 - x0)*(x0*z1 - x1*z0) + (y1 - y0)*(y0*z1 - y1*z0))
        k = (y0*z1 - y1*z0)**2 + (x0*z1 - x1*z0)**2 + (x0*y1 - x1*y0)**2 \
            - r**2*(dx**2 + dy**2 + dz**2)

        return (a, b, c, d, e, f, g, h, j, k)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs

    @classmethod
    def from_points(cls, p1, p2, r=1., **kwargs):
        """Return a cylinder given points that define the axis and a radius.

        Parameters
        ----------
        p1, p2 : 3-tuples
            Points that pass through the plane
        r : float, optional
            Radius of the cylinder. Defaults to 1.
        kwargs : dict
            Keyword arguments passed to the :class:`Quadric` constructor

        Returns
        -------
        Cylinder
            Cylinder that has an axis through the points p1 and p2

        """
        # Convert to numpy arrays
        p1 = np.asarray(p1)
        p2 = np.asarray(p2)
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        dx, dy, dz = p2 - p1

        # Set coefficients for Quadric surface that represents a cylinder of
        # radius r whose axis is the line defined by p1 and p2
        a = dy**2 + dz**2
        b = dx**2 + dz**2
        c = dx**2 + dy**2
        d = -2*dx*dy
        e = -2*dy*dz
        f = -2*dx*dz
        g = -2*((z2 - z1)*(x1*z2 - x2*z1) + (y2 - y1)*(x1*y2-x2*y1))
        h = 2*((x2 - x1)*(x1*y2 - x2*y1) - (z2 - z1)*(y1*z2 - y2*z1))
        j = 2*((x2 - x1)*(x1*z2 - x2*z1) + (y2 - y1)*(y1*z2 - y2*z1))
        k = (y1*z2 - y2*z1)**2 + (x1*z2 - x2*z1)**2 + (x1*y2 - x2*y1)**2 \
            - r**2*(dx**2 + dy**2 + dz**2)

        cyl = cls(x0=x1, y0=y1, z0=z1, r=r, **kwargs)
        cyl._update_from_base_coeffs((a, b, c, d, e, f, g, h, j, k))
        return cyl


class XCylinder(QuadricMeta, Surface):
    """An infinite cylinder whose length is parallel to the x-axis of the form
    :math:`(y - y_0)^2 + (z - z_0)^2 = r^2`.

    Parameters
    ----------
    y0 : float, optional
        y-coordinate for the origin of the Cylinder. Defaults to 0
    z0 : float, optional
        z-coordinate for the origin of the Cylinder. Defaults to 0
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
    y0 : float
        y-coordinate for the origin of the Cylinder
    z0 : float
        z-coordinate for the origin of the Cylinder
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

    _type = 'x-cylinder'
    _coeff_keys = ('y0', 'z0', 'r')

    def __init__(self, y0=0., z0=0., r=1., **kwargs):
        R = kwargs.pop('R', None)
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'),
                 FutureWarning)
            r = R
        super().__init__(**kwargs)
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

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs

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


class YCylinder(QuadricMeta):
    """An infinite cylinder whose length is parallel to the y-axis of the form
    :math:`(x - x_0)^2 + (z - z_0)^2 = r^2`.

    Parameters
    ----------
    x0 : float, optional
        x-coordinate for the origin of the Cylinder. Defaults to 0
    z0 : float, optional
        z-coordinate for the origin of the Cylinder. Defaults to 0
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
        x-coordinate for the origin of the Cylinder
    z0 : float
        z-coordinate for the origin of the Cylinder
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

    _type = 'y-cylinder'
    _coeff_keys = ('x0', 'z0', 'r')

    def __init__(self, x0=0., z0=0., r=1., **kwargs):
        R = kwargs.pop('R', None)
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'), FutureWarning)
            r = R
        super().__init__(**kwargs)
        self.x0 = x0
        self.z0 = z0
        self.r = r

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

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs

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


class ZCylinder(QuadricMeta, Surface):
    """An infinite cylinder whose length is parallel to the z-axis of the form
    :math:`(x - x_0)^2 + (y - y_0)^2 = r^2`.

    Parameters
    ----------
    x0 : float, optional
        x-coordinate for the origin of the Cylinder. Defaults to 0
    y0 : float, optional
        y-coordinate for the origin of the Cylinder. Defaults to 0
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
        x-coordinate for the origin of the Cylinder
    y0 : float
        y-coordinate for the origin of the Cylinder
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

    _type = 'z-cylinder'
    _coeff_keys = ('x0', 'y0', 'r')

    def __init__(self, x0=0., y0=0., r=1., **kwargs):
        R = kwargs.pop('R', None)
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'), FutureWarning)
            r = R
        super().__init__(**kwargs)
        self.x0 = x0
        self.y0 = y0
        self.r = r

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

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs

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


class Sphere(QuadricMeta, Surface):
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

    def __init__(self, x0=0., y0=0., z0=0., r=1., **kwargs):
        R = kwargs.pop('R', None)
        if R is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'), FutureWarning)
            r = R
        super().__init__(**kwargs)
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

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs

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


class Cone(QuadricMeta, Surface):
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
    u : float, optional
        x-component of the vector representing the axis of the cone.
        Defaults to 0.
    v : float, optional
        y-component of the vector representing the axis of the cone.
        Defaults to 0.
    w : float, optional
        z-component of the vector representing the axis of the cone.
        Defaults to 1.
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
    u : float
        x-component of the vector representing the axis of the cone.
    v : float
        y-component of the vector representing the axis of the cone.
    w : float
        z-component of the vector representing the axis of the cone.
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

    _type = 'cylinder'
    _coeff_keys = ('x0', 'y0', 'z0', 'r2', 'u', 'v', 'w')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., u=0., v=0., w=1., **kwargs):
        R2 = kwargs.pop('R2', None)
        if R2 is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r2', 'R2'),
                 FutureWarning)
            r2 = R2
        super().__init__(**kwargs)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.r2 = r2
        self.u = u
        self.v = v
        self.w = w

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

    @property
    def u(self):
        return self.coefficients['u']

    @property
    def v(self):
        return self.coefficients['v']

    @property
    def w(self):
        return self.coefficients['w']

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

    @u.setter
    def u(self, u):
        check_type('u coefficient', u, Real)
        self._coefficients['u'] = u

    @v.setter
    def v(self, v):
        check_type('v coefficient', v, Real)
        self._coefficients['v'] = v

    @w.setter
    def w(self, w):
        check_type('w coefficient', w, Real)
        self._coefficients['w'] = w

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs


class XCone(QuadricMeta, Surface):
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
    _coeff_keys = ('x0', 'y0', 'z0', 'r2')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., **kwargs):
        R2 = kwargs.pop('R2', None)
        if R2 is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r2', 'R2'),
                 FutureWarning)
            r2 = R2
        super().__init__(**kwargs)
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

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs


class YCone(QuadricMeta, Surface):
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
    _coeff_keys = ('x0', 'y0', 'z0', 'r2')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., **kwargs):
        R2 = kwargs.pop('R2', None)
        if R2 is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r2', 'R2'),
                 FutureWarning)
            r2 = R2
        super().__init__(**kwargs)
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

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs


class ZCone(QuadricMeta, Surface):
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
    _coeff_keys = ('x0', 'y0', 'z0', 'r2')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., **kwargs):
        R2 = kwargs.pop('R2', None)
        if R2 is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r2', 'R2'),
                 FutureWarning)
            r2 = R2
        super().__init__(**kwargs)
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

    def _get_base_coeffs(self):
        """Return generalized coefficients for a cylinder"""
        return (0., 0., 1., self.z0)

    def _update_from_base_coeffs(self, coeffs):
        a, b, c, d, e, f, g, h, j, k = coeffs


class Quadric(QuadricMeta, Surface):
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
                 k=0., **kwargs):
        super().__init__(**kwargs)
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

    def remove_redundant_surfaces(self, redundant_surfaces):
        """Recursively remove all redundant surfaces referenced by this region

        Parameters
        ----------
        redundant_surfaces : dict
            Dictionary mapping redundant surface IDs to surface IDs for the 
            :class:`openmc.Surface` instances that should replace them.

        """

        surf = redundant_surfaces.get(self.surface.id)
        if surf is not None:
            self.surface = surf

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
