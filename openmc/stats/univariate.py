from abc import ABCMeta, abstractmethod
from collections import Iterable
from numbers import Real
import sys
from xml.etree import ElementTree as ET

import openmc.checkvalue as cv

if sys.version_info[0] >= 3:
    basestring = str


class Univariate(object):
    """Probability distribution of a single random variable.

    The Univariate class is an abstract class that can be derived to implement a
    specific probability distribution.

    """

    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def to_xml(self):
        return ''


class Discrete(Univariate):
    """Distribution characterized by a probability mass function.

    The Discrete distribution assigns probability values to discrete values of a
    random variable, rather than expressing the distribution as a continuous
    random variable.

    Parameters
    ----------
    x : Iterable of float
        Values of the random variable
    p : Iterable of float
        Discrete probability for each value

    Attributes
    ----------
    x : Iterable of float
        Values of the random variable
    p : Iterable of float
        Discrete probability for each value

    """

    def __init__(self, x, p):
        super(Discrete, self).__init__()
        self.x = x
        self.p = p

    @property
    def x(self):
        return self._x

    @property
    def p(self):
        return self._p

    @x.setter
    def x(self, x):
        if cv._isinstance(x, Real):
            x = [x]
        cv.check_type('discrete values', x, Iterable, Real)
        self._x = x

    @p.setter
    def p(self, p):
        if cv._isinstance(p, Real):
            p = [p]
        cv.check_type('discrete probabilities', p, Iterable, Real)
        for pk in p:
            cv.check_greater_than('discrete probability', pk, 0.0, True)
        self._p = p

    def to_xml(self, element_name):
        element = ET.Element(element_name)
        element.set("type", "discrete")

        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.x)) + ' ' + ' '.join(map(str, self.p))

        return element


class Uniform(Univariate):
    """Distribution with constant probability over a finite interval [a,b]

    Parameters
    ----------
    a : float, optional
        Lower bound of the sampling interval. Defaults to zero.
    b : float, optional
        Upper bound of the sampling interval. Defaults to unity.

    Attributes
    ----------
    a : float
        Lower bound of the sampling interval
    b : float
        Upper bound of the sampling interval

    """

    def __init__(self, a=0.0, b=1.0):
        super(Uniform, self).__init__()
        self.a = a
        self.b = b

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @a.setter
    def a(self, a):
        cv.check_type('Uniform a', a, Real)
        self._a = a

    @b.setter
    def b(self, b):
        cv.check_type('Uniform b', b, Real)
        self._b = b

    def to_xml(self, element_name):
        element = ET.Element(element_name)
        element.set("type", "uniform")
        element.set("parameters", '{} {}'.format(self.a, self.b))
        return element


class Maxwell(Univariate):
    """Maxwellian distribution in energy.

    The Maxwellian distribution in energy is characterized by a single parameter
    :math:`\theta` and has a density function :math:`p(E) dE = c E e^{-E/\theta}
    dE`.

    Parameters
    ----------
    theta : float
        Effective temperature for distribution

    Attributes
    ----------
    theta : float
        Effective temperature for distribution

    """

    def __init__(self, theta):
        super(Maxwell, self).__init__()
        self.theta = theta

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, theta):
        cv.check_type('Maxwell temperature', theta, Real)
        cv.check_greater_than('Maxwell temperature', theta, 0.0)
        self._theta = theta

    def to_xml(self, element_name):
        element = ET.Element(element_name)
        element.set("type", "maxwell")
        element.set("parameters", str(self.theta))
        return element


class Watt(Univariate):
    """Watt fission energy spectrum.

    The Watt fission energy spectrum is characterized by two parameters
    :math:`a` and :math:`b` and has density function :math:`p(E) dE = c e^{-E/a}
    \sinh \sqrt{b \, E} dE`.

    Parameters
    ----------
    a : float
        First parameter of distribution
    b : float
        Second parameter of distribution

    Attributes
    ----------
    a : float
        First parameter of distribution
    b : float
        Second parameter of distribution

    """

    def __init__(self, a=0.988, b=2.249):
        super(Watt, self).__init__()
        self.a = a
        self.b = b

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @a.setter
    def a(self, a):
        cv.check_type('Watt a', a, Real)
        cv.check_greater_than('Watt a', a, 0.0)
        self._a = a

    @b.setter
    def b(self, b):
        cv.check_type('Watt b', b, Real)
        cv.check_greater_than('Watt b', b, 0.0)
        self._b = b

    def to_xml(self, element_name):
        element = ET.Element(element_name)
        element.set("type", "watt")
        element.set("parameters", '{} {}'.format(self.a, self.b))
        return element


class Tabular(Univariate):
    """Piecewise continuous probability distribution.

    This class is used to represent a probability distribution whose density
    function is tabulated at specific values and is either histogram or linearly
    interpolated between points.

    Parameters
    ----------
    x : Iterable of float
        Tabulated values of the random variable
    p : Iterable of float
        Tabulated probabilities
    interpolation : {'histogram', 'linear-linear'}, optional
        Indicate whether the density function is constant between tabulated
        points or linearly-interpolated.

    Attributes
    ----------
    x : Iterable of float
        Tabulated values of the random variable
    p : Iterable of float
        Tabulated probabilities
    interpolation : {'histogram', 'linear-linear'}, optional
        Indicate whether the density function is constant between tabulated
        points or linearly-interpolated.

    """

    def __init__(self, x, p, interpolation='linear-linear'):
        super(Tabular, self).__init__()
        self.x = x
        self.p = p
        self.interpolation = interpolation

    @property
    def x(self):
        return self._x

    @property
    def p(self):
        return self._p

    @property
    def interpolation(self):
        return self._interpolation

    @x.setter
    def x(self, x):
        cv.check_type('tabulated values', x, Iterable, Real)
        self._x = x

    @p.setter
    def p(self, p):
        cv.check_type('tabulated probabilities', p, Iterable, Real)
        for pk in p:
            cv.check_greater_than('tabulated probability', pk, 0.0, True)
        self._p = p

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_value('interpolation', interpolation,
                       ['linear-linear', 'histogram'])
        self._interpolation = interpolation

    def to_xml(self, element_name):
        element = ET.Element(element_name)
        element.set("type", "tabular")
        element.set("interpolation", self.interpolation)

        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.x)) + ' ' + ' '.join(map(str, self.p))

        return element
