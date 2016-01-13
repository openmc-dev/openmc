from abc import ABCMeta, abstractmethod
from collections import Iterable
from math import pi
from numbers import Real
import sys
from xml.etree import ElementTree as ET

import numpy as np
from numpy.linalg import norm

import openmc.checkvalue as cv
from openmc.stats.univariate import Univariate, Uniform

if sys.version_info[0] >= 3:
    basestring = str


class UnitSphere(object):
    """Distribution of points on the unit sphere.

    This abstract class is used for angular distributions, since a direction is
    represented as a unit vector (i.e., vector on the unit sphere).

    Parameters
    ----------
    name : str
        Name of the distribution
    reference_uvw : Iterable of Real
        Direction from which polar angle is measured

    Attributes
    ----------
    name : str
        Name of the distribution
    reference_uvw : Iterable of Real
        Direction from which polar angle is measured

    """

    __metaclass__ = ABCMeta

    def __init__(self, name, reference_uvw=None):
        self.name = name
        self._reference_uvw = None
        if reference_uvw is not None:
            self.reference_uvw = reference_uvw

    @property
    def name(self):
        return self._name

    @property
    def reference_uvw(self):
        return self._reference_uvw

    @name.setter
    def name(self, name):
        cv.check_type('name', name, basestring)
        self._name = name

    @reference_uvw.setter
    def reference_uvw(self, uvw):
        cv.check_type('reference direction', uvw, Iterable, Real)
        uvw = np.asarray(uvw)
        self._reference_uvw = uvw/norm(uvw)

    @abstractmethod
    def to_xml(self):
        return ''


class PolarAzimuthal(UnitSphere):
    """Angular distribution represented by polar and azimuthal angles

    This distribution allows one to specify the distribution of the cosine of
    the polar angle and the azimuthal angle independently of once another.

    Parameters
    ----------
    mu : openmc.stats.Univariate
        Distribution of the cosine of the polar angle
    phi : openmc.stats.Univariate
        Distribution of the azimuthal angle
    name : str, optional
        Name of the distribution. Defaults to 'angle'.
    reference_uvw : Iterable of Real
        Direction from which polar angle is measured. Defaults to the positive
        z-direction.

    Attributes
    ----------
    mu : openmc.stats.Univariate
        Distribution of the cosine of the polar angle
    phi : openmc.stats.Univariate
        Distribution of the azimuthal angle

    """

    def __init__(self, mu=None, phi=None, name='angle',reference_uvw=[0., 0., 1.]):
        super(PolarAzimuthal, self).__init__(name, reference_uvw)
        if mu is not None:
            self.mu = mu
        else:
            self.mu = Uniform('mu', -1., 1.)

        if phi is not None:
            self.phi = phi
        else:
            self.phi = Uniform('phi', 0., 2*pi)

    @property
    def mu(self):
        return self._mu

    @property
    def phi(self):
        return self._phi

    @mu.setter
    def mu(self, mu):
        cv.check_type('cosine of polar angle', mu, Univariate)
        self._mu = mu

    @phi.setter
    def phi(self, phi):
        cv.check_type('azimuthal angle', phi, Univariate)
        self._phi = phi

    def to_xml(self):
        element = ET.Element(self.name)
        element.set("type", "mu-phi")
        if self.reference_uvw is not None:
            element.set("reference_uvw", ' '.join(map(str, self.reference_uvw)))
        element.append(self.mu.to_xml())
        element.append(self.phi.to_xml())
        return element


class Isotropic(UnitSphere):
    """Isotropic angular distribution.

    Parameters
    ----------
    name : str, optional
        Name of the distribution. Defaults to 'angle'.

    """

    def __init__(self, name='angle'):
        super(Isotropic, self).__init__(name)

    def to_xml(self):
        element = ET.Element(self.name)
        element.set("type", "isotropic")
        return element


class Monodirectional(UnitSphere):
    """Monodirectional angular distribution.

    A monodirectional angular distribution is one for which the polar and
    azimuthal angles are always the same. It is completely specified by the
    reference direction vector.

    Parameters
    ----------
    name : str, optional
        Name of the distribution. Defaults to 'angle'.
    reference_uvw : Iterable of Real
        Direction from which polar angle is measured. Defaults to the positive
        x-direction.

    """


    def __init__(self, name='angle', reference_uvw=[1., 0., 0.]):
        super(Monodirectional, self).__init__(name, reference_uvw)

    def to_xml(self):
        element = ET.Element(self.name)
        element.set("type", "monodirectional")
        if self.reference_uvw is not None:
            element.set("reference_uvw", ' '.join(map(str, self.reference_uvw)))
        return element


class Spatial(object):
    """Distribution of locations in three-dimensional Euclidean space.

    Classes derived from this abstract class can be used for spatial
    distributions of source sites.

    Parameters
    ----------
    name : str
        Name of the distribution

    Attributes
    ----------
    name : str
        Name of the distribution

    """

    __metaclass__ = ABCMeta

    def __init__(self, name):
        self.name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        cv.check_type('name', name, basestring)
        self._name = name

    @abstractmethod
    def to_xml(self):
        return ''


class SpatialIndependent(Spatial):
    """Spatial distribution with independent x, y, and z distributions.

    This distribution allows one to specify a coordinates whose x-, y-, and z-
    components are sampled independently from one another.

    Parameters
    ----------
    x : openmc.stats.Univariate
        Distribution of x-coordinates
    y : openmc.stats.Univariate
        Distribution of y-coordinates
    z : openmc.stats.Univariate
        Distribution of z-coordinates
    name : str
        Name of the distribution

    Attributes
    ----------
    x : openmc.stats.Univariate
        Distribution of x-coordinates
    y : openmc.stats.Univariate
        Distribution of y-coordinates
    z : openmc.stats.Univariate
        Distribution of z-coordinates

    """


    def __init__(self, x, y, z, name='space'):
        super(SpatialIndependent, self).__init__(name)
        self.x = x
        self.y = y
        self.z = z

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @x.setter
    def x(self, x):
        cv.check_type('x coordinate', x, Univariate)
        self._x = x

    @y.setter
    def y(self, y):
        cv.check_type('y coordinate', y, Univariate)
        self._y = y

    @z.setter
    def z(self, z):
        cv.check_type('z coordinate', z, Univariate)
        self._z = z

    def to_xml(self):
        element = ET.Element(self.name)
        element.set("type", "independent")
        element.append(self.x.to_xml())
        element.append(self.y.to_xml())
        element.append(self.z.to_xml())
        return element


class SpatialBox(Spatial):
    """Uniform distribution of coordinates in a rectangular cuboid.

    Parameters
    ----------
    lower_left : Iterable of Real
        Lower-left coordinates of cuboid
    upper_right : Iterable of Real
        Upper-right coordinates of cuboid
    name : str, optional
        Name of the distribution
    only_fissionable : bool, optional
        Whether spatial sites should only be accepted if they occur in
        fissionable materials

    Attributes
    ----------
    lower_left : Iterable of Real
        Lower-left coordinates of cuboid
    upper_right : Iterable of Real
        Upper-right coordinates of cuboid
    only_fissionable : bool, optional
        Whether spatial sites should only be accepted if they occur in
        fissionable materials

    """


    def __init__(self, lower_left, upper_right, name='space', only_fissionable=False):
        super(SpatialBox, self).__init__(name)
        self.lower_left = lower_left
        self.upper_right = upper_right
        self.only_fissionable = only_fissionable

    @property
    def lower_left(self):
        return self._lower_left

    @property
    def upper_right(self):
        return self._upper_right

    @property
    def only_fissionable(self):
        return self._only_fissionable

    @lower_left.setter
    def lower_left(self, lower_left):
        cv.check_type('lower left coordinate', lower_left, Iterable, Real)
        cv.check_length('lower left coordinate', lower_left, 3)
        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):
        cv.check_type('upper right coordinate', upper_right, Iterable, Real)
        cv.check_length('upper right coordinate', upper_right, 3)
        self._upper_right = upper_right

    @only_fissionable.setter
    def only_fissionable(self, only_fissionable):
        cv.check_type('only fissionable', only_fissionable, bool)
        self._only_fissionable = only_fissionable

    def to_xml(self):
        element = ET.Element(self.name)
        if self.only_fissionable:
            element.set("type", "fission")
        else:
            element.set("type", "box")
        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.lower_left)) + ' ' + \
                      ' '.join(map(str, self.upper_right))
        return element


class SpatialPoint(Spatial):
    """Delta function in three dimensions.

    This spatial distribution can be used for a point source where sites are
    emitted at a specific location given by its Cartesian coordinates.

    Parameters
    ----------
    xyz : Iterable of Real
        Cartesian coordinates of location
    name : str, optional
        Name of the distribution

    Attributes
    ----------
    xyz : Iterable of Real
        Cartesian coordinates of location

    """

    def __init__(self, xyz, name='space'):
        super(SpatialPoint, self).__init__(name)
        self.xyz = xyz

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        cv.check_type('coordinate', xyz, Iterable, Real)
        cv.check_length('coordinate', xyz, 3)
        self._xyz = xyz

    def to_xml(self):
        element = ET.Element(self.name)
        element.set("type", "point")
        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.xyz))
        return element
