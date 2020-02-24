from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from numbers import Real
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

_WARNING_KWARGS = """\
"{}(...) accepts keyword arguments only for '{}'. Future versions of OpenMC \
will not accept positional parameters for superclass arguments.\
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
    _atol = 1.e-12

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

    def normalize(self, coeffs=None):
        """Normalize coefficients by first nonzero value

        Parameters
        ----------
        coeffs : tuple, optional
            Tuple of surface coefficients to normalize. Defaults to None. If no
            coefficients are supplied then the coefficients will be taken from
            the current Surface.

        Returns
        -------
        tuple of normalized coefficients

        """
        if coeffs is None:
            coeffs = self._get_base_coeffs()
        coeffs = np.asarray(coeffs)
        nonzeros = ~np.isclose(coeffs, 0., rtol=0., atol=self._atol)
        norm_factor = np.abs(coeffs[nonzeros][0])
        return tuple([c/norm_factor for c in coeffs])

    def is_equal(self, other):
        """Determine if this Surface is equivalent to another

        Parameters
        ----------
        other : instance of openmc.Surface
            Instance of openmc.Surface that should be compared to the current
            surface

        """
        coeffs1 = self.normalize(self._get_base_coeffs())
        coeffs2 = self.normalize(other._get_base_coeffs())

        return np.allclose(coeffs1, coeffs2, rtol=0., atol=self._atol)

    @abstractmethod
    def _get_base_coeffs(self):
        """Return polynomial coefficients representing the implicit surface
        equation.

        """

    @abstractmethod
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
            Evaluation of the surface polynomial at point :math:`(x',y',z')`

        """

    @abstractmethod
    def translate(self, vector, inplace=False):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated
        inplace : boolean
            Whether or not to return a new instance of this Surface or to
            modify the coefficients of this Surface. Defaults to False

        Returns
        -------
        instance of openmc.Surface
            Translated surface

        """

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


class PlaneMixin(metaclass=ABCMeta):
    """A Plane mixin class for all operations on order 1 surfaces"""
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

    def translate(self, vector, inplace=False):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated
        inplace : boolean
            Whether or not to return a new instance of a Plane or to modify the
            coefficients of this plane. Defaults to False

        Returns
        -------
        openmc.Plane
            Translated surface

        """
        vx, vy, vz = vector
        a, b, c, d = self._get_base_coeffs()
        d = d + a*vx + b*vy + c*vz

        surf = self if inplace else self.clone()

        setattr(surf, surf._coeff_keys[-1], d)

        return surf

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


class Plane(PlaneMixin, Surface):
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

    def __init__(self, a=1., b=0., c=0., d=0., *args, **kwargs):
        # *args should ultimately be limited to a, b, c, d as specified in
        # __init__, but to preserve the API it is allowed to accept Surface
        # parameters for now, but will raise warnings if this is done.
        argtup = ('a', 'b', 'c', 'd', 'boundary_type', 'name', 'surface_id')
        kwargs.update(dict(zip(argtup, args)))

        # Warn if Surface parameters are passed by position, not by keyword
        superkwargs = {}
        for k in ('boundary_type', 'name', 'surface_id'):
            val = kwargs.get(k, None)
            if val is not None:
                superkwargs[k] = val
                warn(_WARNING_KWARGS.format(type(self), k),
                     FutureWarning)

        super().__init__(**superkwargs)

        for key, val in zip(self._coeff_keys, (a, b, c, d)):
            setattr(self, key, val)

        # Warn if capital letter arguments are passed
        for k in 'ABCD':
            val = kwargs.pop(k, None)
            if val is not None:
                warn(_WARNING_UPPER.format(type(self), k.lower(), k),
                     FutureWarning)
                setattr(self, k.lower(), val)

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


class XPlane(PlaneMixin, Surface):
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

    def __init__(self, x0=0., *args, **kwargs):
        # work around for accepting Surface kwargs as positional parameters
        # until they are deprecated
        argsdict = dict(zip(('boundary_type', 'name', 'surface_id'), args))
        for k in argsdict:
            warn(_WARNING_KWARGS.format(type(self).__name__, k), FutureWarning)
        kwargs.update(argsdict)

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

    def bounding_box(self, side):
        if side == '-':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([self.x0, np.inf, np.inf]))
        elif side == '+':
            return (np.array([self.x0, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))

    def evaluate(self, point):
        return point[0] - self.x0


Plane.register(XPlane)


class YPlane(PlaneMixin, Surface):
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

    def __init__(self, y0=0., *args, **kwargs):
        # work around for accepting Surface kwargs as positional parameters
        # until they are deprecated
        argsdict = dict(zip(('boundary_type', 'name', 'surface_id'), args))
        for k in argsdict:
            warn(_WARNING_KWARGS.format(type(self).__name__, k), FutureWarning)
        kwargs.update(argsdict)

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
        return (0., 1., 0., self.y0)

    def bounding_box(self, side):
        if side == '-':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, self.y0, np.inf]))
        elif side == '+':
            return (np.array([-np.inf, self.y0, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))

    def evaluate(self, point):
        return point[1] - self.y0


Plane.register(YPlane)


class ZPlane(PlaneMixin, Surface):
    """A plane perpendicular to the z axis of the form :math:`z - z_0 = 0`

    Parameters
    ----------
    z0 : float, optional
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

    def __init__(self, z0=0., *args, **kwargs):
        # work around for accepting Surface kwargs as positional parameters
        # until they are deprecated
        argsdict = dict(zip(('boundary_type', 'name', 'surface_id'), args))
        for k in argsdict:
            warn(_WARNING_KWARGS.format(type(self).__name__, k), FutureWarning)
        kwargs.update(argsdict)

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
        return (0., 0., 1., self.z0)

    def bounding_box(self, side):
        if side == '-':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, self.z0]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, self.z0]),
                    np.array([np.inf, np.inf, np.inf]))

    def evaluate(self, point):
        return point[2] - self.z0


Plane.register(ZPlane)


class QuadricMixin(metaclass=ABCMeta):
    """A Mixin class implementing common functionality for quadric surfaces"""

    def get_Abc(self, coeffs=None):
        """Compute matrix, vector, and scalar coefficients for this surface or
        for a specified set of coefficients.

        Parameters
        ----------
        coeffs : tuple, optional
            Tuple of coefficients from which to compute the quadric elements.
            If none are supplied the coefficients of this surface will be used.
        """
        if coeffs is None:
            a, b, c, d, e, f, g, h, j, k = self._get_base_coeffs()
        else:
            a, b, c, d, e, f, g, h, j, k = coeffs

        A = np.array([[a, d/2, f/2], [d/2, b, e/2], [f/2, e/2, c]])
        bvec = np.array([g, h, j])

        return A, bvec, k

    def eigh(self, coeffs=None):
        """Wrapper method for returning eigenvalues and eigenvectors of this
        quadric surface which is used for transformations.

        Parameters
        ----------
        coeffs : tuple, optional
            Tuple of coefficients from which to compute the quadric elements.
            If none are supplied the coefficients of this surface will be used.

        Returns
        -------
        w, v : tuple of numpy arrays with shapes (3,) and (3,3) respectively
            Returns the eigenvalues and eigenvectors of the quadric matrix A
            that represents the supplied coefficients. The vector w contains
            the eigenvalues in ascending order and the matrix v contains the
            eigenvectors such that v[:,i] is the eigenvector corresponding to
            the eigenvalue w[i].

        """
        return np.linalg.eigh(self.get_Abc(coeffs=coeffs)[0])

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
        x = np.asarray(point)
        A, b, c = self.get_Abc()
        return np.matmul(x.T, np.matmul(A, x)) + np.matmul(b.T, x) + c

    def translate(self, vector, inplace=False):
        """Translate surface in given direction

        Parameters
        ----------
        vector : iterable of float
            Direction in which surface should be translated
        inplace : boolean
            Whether to return a clone of the Surface or the Surface itself.
            Defaults to False

        Returns
        -------
        openmc.Surface
            Translated surface

        """
        vector = np.asarray(vector)

        surf = self if inplace else self.clone()

        if set(('x0', 'y0', 'z0')).intersection(set(surf._coeff_keys)):
            for vi, xi in zip(vector, ('x0', 'y0', 'z0')):
                val = getattr(surf, xi, None)
                if val is not None:
                    setattr(surf, xi, val + vi)
        else:
            A, bvec, cnst = self.get_Abc()

            g, h, j = bvec - 2*np.matmul(vector.T, A)
            k = cnst + np.matmul(vector.T, np.matmul(A, vector)) \
                - np.matmul(bvec.T, vector)

            for key, val in zip(('g', 'h', 'j', 'k'), (g, h, j, k)):
                setattr(surf, key, val)

        return surf


class Cylinder(QuadricMixin, Surface):
    """A cylinder with radius r, centered on the point (x0, y0, z0) with an
    axis specified by the line through points (x0, y0, z0) and (x0+dx, y0+dy,
    z0+dz)

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
    dx : float, optional
        x-component of the vector representing the axis of the cylinder.
        Defaults to 0.
    dy : float, optional
        y-component of the vector representing the axis of the cylinder.
        Defaults to 0.
    dz : float, optional
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
    dx : float
        x-component of the vector representing the axis of the cylinder
    dy : float
        y-component of the vector representing the axis of the cylinder
    dz : float
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
    _coeff_keys = ('x0', 'y0', 'z0', 'r', 'dx', 'dy', 'dz')

    def __init__(self, x0=0., y0=0., z0=0., r=1., dx=0., dy=0., dz=1., **kwargs):
        raise NotImplementedError('There is no C++ implementation for general '
                                  'Cylinders yet, please use '
                                  'openmc.model.funcs.cylinder_from_points to '
                                  'return a Quadric instance instead for now')

        super().__init__(**kwargs)

        for key, val in zip(self._coeff_keys, (x0, y0, z0, r, dx, dy, dz)):
            setattr(self, key, val)

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
    def dx(self):
        return self.coefficients['dx']

    @property
    def dy(self):
        return self.coefficients['dy']

    @property
    def dz(self):
        return self.coefficients['dz']

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

    @dx.setter
    def dx(self, dx):
        check_type('dx coefficient', dx, Real)
        self._coefficients['dx'] = dx

    @dy.setter
    def dy(self, dy):
        check_type('dy coefficient', dy, Real)
        self._coefficients['dy'] = dy

    @dz.setter
    def dz(self, dz):
        check_type('dz coefficient', dz, Real)
        self._coefficients['dz'] = dz

    def _get_base_coeffs(self):
        # Get x, y, z coordinates of two points
        x1, y1, z1 = self.x0, self.y0, self.z0
        x2, y2, z2 = x1 + self.dx, y1 + self.dy, z1 + self.dz
        r = self.r 

        # Define intermediate terms
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        cx = y1*z2 - y2*z1
        cy = x2*z1 - x1*z2
        cz = x1*y2 - x2*y1

        # Given p=(x,y,z), p1=(x1, y1, z1), p2=(x2, y2, z2), the equation
        # for the cylinder can be derived as
        # r = |(p - p1) ⨯ (p - p2)| / |p2 - p1|.
        # Expanding out all terms and grouping according to what Quadric
        # expects gives the following coefficients.
        a = dy*dy + dz*dz
        b = dx*dx + dz*dz
        c = dx*dx + dy*dy
        d = -2*dx*dy
        e = -2*dy*dz
        f = -2*dx*dz
        g = 2*(cy*dz - cz*dy)
        h = 2*(cz*dx - cx*dz)
        j = 2*(cx*dy - cy*dx)
        k = cx*cx + cy*cy + cz*cz - (dx*dx + dy*dy + dz*dz)*r*r

        return (a, b, c, d, e, f, g, h, j, k)

    @classmethod
    def from_points(cls, p1, p2, r=1., **kwargs):
        """Return a cylinder given points that define the axis and a radius.

        Parameters
        ----------
        p1, p2 : 3-tuples
            Points that pass through the plane, p1 will be used as (x0, y0, z0)
        r : float, optional
            Radius of the cylinder. Defaults to 1.
        kwargs : dict
            Keyword arguments passed to the :class:`Cylinder` constructor

        Returns
        -------
        Cylinder
            Cylinder that has an axis through the points p1 and p2, and a
            radius r.

        """
        raise NotImplementedError('There is no C++ implementation for general '
                                  'Cylinders yet, please use '
                                  'openmc.model.funcs.cylinder_from_points to '
                                  'return a Quadric instance instead for now')
        # Convert to numpy arrays
        p1 = np.asarray(p1)
        p2 = np.asarray(p2)
        x0, y0, z0 = p1
        dx, dy, dz = p2 - p1

        return cls(x0=x0, y0=y0, z0=z0, r=r, dx=dx, dy=dy, dz=dz, **kwargs)


class XCylinder(QuadricMixin, Surface):
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

        for key, val in zip(self._coeff_keys, (y0, z0, r)):
            setattr(self, key, val)

    @property
    def y0(self):
        return self.coefficients['y0']

    @property
    def z0(self):
        return self.coefficients['z0']

    @property
    def r(self):
        return self.coefficients['r']

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
        y0, z0, r = self.y0, self.z0, self.r

        a = d = e = f = g = 0.
        b = c = 1.
        h, j, k = -2*y0, -2*z0, y0*y0 + z0*z0 - r*r

        return (a, b, c, d, e, f, g, h, j, k)

    def bounding_box(self, side):
        if side == '-':
            return (np.array([-np.inf, self.y0 - self.r, self.z0 - self.r]),
                    np.array([np.inf, self.y0 + self.r, self.z0 + self.r]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))

    def evaluate(self, point):
        y = point[1] - self.y0
        z = point[2] - self.z0
        return y*y + z*z - self.r**2


Cylinder.register(XCylinder)


class YCylinder(QuadricMixin, Surface):
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
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'),
                 FutureWarning)
            r = R
        super().__init__(**kwargs)

        for key, val in zip(self._coeff_keys, (x0, z0, r)):
            setattr(self, key, val)

    @property
    def x0(self):
        return self.coefficients['x0']

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

    @z0.setter
    def z0(self, z0):
        check_type('z0 coefficient', z0, Real)
        self._coefficients['z0'] = z0

    @r.setter
    def r(self, r):
        check_type('r coefficient', r, Real)
        self._coefficients['r'] = r

    def _get_base_coeffs(self):
        x0, z0, r = self.x0, self.z0, self.r

        b = d = e = f = h = 0.
        a = c = 1.
        g, j, k = -2*x0, -2*z0, x0*x0 + z0*z0 - r*r

        return (a, b, c, d, e, f, g, h, j, k)

    def bounding_box(self, side):
        if side == '-':
            return (np.array([self.x0 - self.r, -np.inf, self.z0 - self.r]),
                    np.array([self.x0 + self.r, np.inf, self.z0 + self.r]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))

    def evaluate(self, point):
        x = point[0] - self.x0
        z = point[2] - self.z0
        return x*x + z*z - self.r**2


Cylinder.register(YCylinder)


class ZCylinder(QuadricMixin, Surface):
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
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'),
                 FutureWarning)
            r = R
        super().__init__(**kwargs)

        for key, val in zip(self._coeff_keys, (x0, y0, r)):
            setattr(self, key, val)

    @property
    def x0(self):
        return self.coefficients['x0']

    @property
    def y0(self):
        return self.coefficients['y0']

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

    @r.setter
    def r(self, r):
        check_type('r coefficient', r, Real)
        self._coefficients['r'] = r

    def _get_base_coeffs(self):
        x0, y0, r = self.x0, self.y0, self.r

        c = d = e = f = j = 0.
        a = b = 1.
        g, h, k = -2*x0, -2*y0, x0*x0 + y0*y0 - r*r

        return (a, b, c, d, e, f, g, h, j, k)

    def bounding_box(self, side):
        if side == '-':
            return (np.array([self.x0 - self.r, self.y0 - self.r, -np.inf]),
                    np.array([self.x0 + self.r, self.y0 + self.r, np.inf]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))

    def evaluate(self, point):
        x = point[0] - self.x0
        y = point[1] - self.y0
        return x*x + y*y - self.r**2


Cylinder.register(ZCylinder)


class Sphere(QuadricMixin, Surface):
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
            warn(_WARNING_UPPER.format(type(self).__name__, 'r', 'R'),
                 FutureWarning)
            r = R
        super().__init__(**kwargs)

        for key, val in zip(self._coeff_keys, (x0, y0, z0, r)):
            setattr(self, key, val)

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
        x0, y0, z0, r = self.x0, self.y0, self.z0, self.r
        a = b = c = 1.
        d = e = f = 0.
        g, h, j = -2*x0, -2*y0, -2*z0
        k = x0*x0 + y0*y0 + z0*z0 - r*r

        return (a, b, c, d, e, f, g, h, j, k)

    def bounding_box(self, side):
        if side == '-':
            return (np.array([self.x0 - self.r, self.y0 - self.r,
                              self.z0 - self.r]),
                    np.array([self.x0 + self.r, self.y0 + self.r,
                              self.z0 + self.r]))
        elif side == '+':
            return (np.array([-np.inf, -np.inf, -np.inf]),
                    np.array([np.inf, np.inf, np.inf]))

    def evaluate(self, point):
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return x*x + y*y + z*z - self.r**2


class Cone(QuadricMixin, Surface):
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
    dx : float, optional
        x-component of the vector representing the axis of the cone.
        Defaults to 0.
    dy : float, optional
        y-component of the vector representing the axis of the cone.
        Defaults to 0.
    dz : float, optional
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
    dx : float
        x-component of the vector representing the axis of the cone.
    dy : float
        y-component of the vector representing the axis of the cone.
    dz : float
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

    _type = 'cone'
    _coeff_keys = ('x0', 'y0', 'z0', 'r2', 'dx', 'dy', 'dz')

    def __init__(self, x0=0., y0=0., z0=0., r2=1., dx=0., dy=0., dz=1., **kwargs):
        raise NotImplementedError('There is no C++ implementation for general '
                                  'Cones yet, this functionality should be '
                                  'added soon.')
        R2 = kwargs.pop('R2', None)
        if R2 is not None:
            warn(_WARNING_UPPER.format(type(self).__name__, 'r2', 'R2'),
                 FutureWarning)
            r2 = R2
        super().__init__(**kwargs)

        for key, val in zip(self._coeff_keys, (x0, y0, z0, r2, dx, dy, dz)):
            setattr(self, key, val)

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
    def dx(self):
        return self.coefficients['dx']

    @property
    def dy(self):
        return self.coefficients['dy']

    @property
    def dz(self):
        return self.coefficients['dz']

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

    @dx.setter
    def dx(self, dx):
        check_type('dx coefficient', dx, Real)
        self._coefficients['dx'] = dx

    @dy.setter
    def dy(self, dy):
        check_type('dy coefficient', dy, Real)
        self._coefficients['dy'] = dy

    @dz.setter
    def dz(self, dz):
        check_type('dz coefficient', dz, Real)
        self._coefficients['dz'] = dz

    def _get_base_coeffs(self):
        # The equation for a general cone with vertex at point p = (x0, y0, z0)
        # and axis specified by the unit vector d = (dx, dy, dz) and opening
        # half angle theta can be described by the equation
        #
        # (d*(r - p))^2 - (r - p)*(r - p)cos^2(theta) = 0
        #
        # where * is the dot product and the vector r is the evaulation point
        # r = (x, y, z)
        #
        # The argument r2 for cones is actually tan^2(theta) so that
        # cos^2(theta) = 1 / (1 + r2)

        x0, y0, z0, r2 = self.x0, self.y0, self.z0, self.r2
        dx, dy, dz = self.dx, self.dy, self.dz
        dnorm = dx*dx + dy*dy + dz*dz
        dx /= dnorm
        dy /= dnorm
        dz /= dnorm
        cos2 = 1 / (1 + r2)

        a = dx*dx - cos2
        b = dy*dy - cos2
        c = dz*dz - cos2
        d = 2*dx*dy
        e = 2*dy*dz
        f = 2*dx*dz
        g = -2*(dx*dx*x0 + dx*dy*y0 + dx*dz*z0 - cos2)
        h = -2*(dy*dy*y0 + dx*dy*x0 + dy*dz*z0 - cos2)
        j = -2*(dz*dz*y0 + dx*dz*x0 + dy*dz*y0 - cos2)
        k = (dx*x0 + dy*y0 + dz*z0)**2 - cos2*(x0*x0 + y0*y0 + z0*z0)

        return (a, b, c, d, e, f, g, h, j, k)


class XCone(QuadricMixin, Surface):
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

        for key, val in zip(self._coeff_keys, (x0, y0, z0, r2)):
            setattr(self, key, val)

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
        x0, y0, z0, r2 = self.x0, self.y0, self.z0, self.r2

        a = -r2
        b = c = 1.
        d = e = f = 0.
        g, h, j = 2*x0*r2, -2*y0, -2*z0
        k = y0*y0 + z0*z0 - r2*x0*x0

        return (a, b, c, d, e, f, g, h, j, k)

    def evaluate(self, point):
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return y*y + z*z - self.r2*x*x


Cone.register(XCone)


class YCone(QuadricMixin, Surface):
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

        for key, val in zip(self._coeff_keys, (x0, y0, z0, r2)):
            setattr(self, key, val)

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
        x0, y0, z0, r2 = self.x0, self.y0, self.z0, self.r2

        b = -r2
        a = c = 1.
        d = e = f = 0.
        g, h, j = -2*x0, 2*y0*r2, -2*z0
        k = x0*x0 + z0*z0 - r2*y0*y0

        return (a, b, c, d, e, f, g, h, j, k)

    def evaluate(self, point):
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return x*x + z*z - self.r2*y*y


Cone.register(YCone)


class ZCone(QuadricMixin, Surface):
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

        for key, val in zip(self._coeff_keys, (x0, y0, z0, r2)):
            setattr(self, key, val)

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
        x0, y0, z0, r2 = self.x0, self.y0, self.z0, self.r2

        c = -r2
        a = b = 1.
        d = e = f = 0.
        g, h, j = -2*x0, -2*y0, 2*z0*r2
        k = x0*x0 + y0*y0 - r2*z0*z0

        return (a, b, c, d, e, f, g, h, j, k)

    def evaluate(self, point):
        x = point[0] - self.x0
        y = point[1] - self.y0
        z = point[2] - self.z0
        return x*x + y*y - self.r2*z*z


Cone.register(ZCone)


class Quadric(QuadricMixin, Surface):
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

        for key, val in zip(self._coeff_keys, (a, b, c, d, e, f, g, h, j, k)):
            setattr(self, key, val)

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

    def _get_base_coeffs(self):
        return tuple(getattr(self, c) for c in self._coeff_keys)


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


_SURFACE_CLASSES = {cls._type: cls for cls in Surface.__subclasses__()}
