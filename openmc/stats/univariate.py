from abc import ABCMeta, abstractmethod
from collections.abc import Iterable
from numbers import Real
import sys
from xml.etree import ElementTree as ET

import numpy as np

import openmc.checkvalue as cv
from openmc._xml import get_text
from openmc.mixin import EqualityMixin


_INTERPOLATION_SCHEMES = ['histogram', 'linear-linear', 'linear-log',
                          'log-linear', 'log-log']


class Univariate(EqualityMixin, metaclass=ABCMeta):
    """Probability distribution of a single random variable.

    The Univariate class is an abstract class that can be derived to implement a
    specific probability distribution.

    """
    def __init__(self):
        pass

    @abstractmethod
    def to_xml_element(self, element_name):
        return ''

    @abstractmethod
    def __len__(self):
        return 0

    @classmethod
    @abstractmethod
    def from_xml_element(cls, elem):
        distribution = get_text(elem, 'type')
        if distribution == 'discrete':
            return Discrete.from_xml_element(elem)
        elif distribution == 'uniform':
            return Uniform.from_xml_element(elem)
        elif distribution == 'maxwell':
            return Maxwell.from_xml_element(elem)
        elif distribution == 'watt':
            return Watt.from_xml_element(elem)
        elif distribution == 'normal':
            return Normal.from_xml_element(elem)
        elif distribution == 'muir':
            return Muir.from_xml_element(elem)
        elif distribution == 'tabular':
            return Tabular.from_xml_element(elem)
        elif distribution == 'legendre':
            return Legendre.from_xml_element(elem)
        elif distribution == 'mixture':
            return Mixture.from_xml_element(elem)


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
        super().__init__()
        self.x = x
        self.p = p

    def __len__(self):
        return len(self.x)

    @property
    def x(self):
        return self._x

    @property
    def p(self):
        return self._p

    @x.setter
    def x(self, x):
        if isinstance(x, Real):
            x = [x]
        cv.check_type('discrete values', x, Iterable, Real)
        self._x = x

    @p.setter
    def p(self, p):
        if isinstance(p, Real):
            p = [p]
        cv.check_type('discrete probabilities', p, Iterable, Real)
        for pk in p:
            cv.check_greater_than('discrete probability', pk, 0.0, True)
        self._p = p

    def to_xml_element(self, element_name):
        """Return XML representation of the discrete distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing discrete distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "discrete")

        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.x)) + ' ' + ' '.join(map(str, self.p))

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate discrete distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Discrete
            Discrete distribution generated from XML element

        """
        params = [float(x) for x in get_text(elem, 'parameters').split()]
        x = params[:len(params)//2]
        p = params[len(params)//2:]
        return cls(x, p)


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
        super().__init__()
        self.a = a
        self.b = b

    def __len__(self):
        return 2

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

    def to_tabular(self):
        prob = 1./(self.b - self.a)
        t = Tabular([self.a, self.b], [prob, prob], 'histogram')
        t.c = [0., 1.]
        return t

    def to_xml_element(self, element_name):
        """Return XML representation of the uniform distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing uniform distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "uniform")
        element.set("parameters", '{} {}'.format(self.a, self.b))
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate uniform distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Uniform
            Uniform distribution generated from XML element

        """
        params = get_text(elem, 'parameters').split()
        return cls(*map(float, params))


class Maxwell(Univariate):
    """Maxwellian distribution in energy.

    The Maxwellian distribution in energy is characterized by a single parameter
    :math:`\theta` and has a density function :math:`p(E) dE = c E e^{-E/\theta}
    dE`.

    Parameters
    ----------
    theta : float
        Effective temperature for distribution in eV

    Attributes
    ----------
    theta : float
        Effective temperature for distribution in eV

    """

    def __init__(self, theta):
        super().__init__()
        self.theta = theta

    def __len__(self):
        return 1

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, theta):
        cv.check_type('Maxwell temperature', theta, Real)
        cv.check_greater_than('Maxwell temperature', theta, 0.0)
        self._theta = theta

    def to_xml_element(self, element_name):
        """Return XML representation of the Maxwellian distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing Maxwellian distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "maxwell")
        element.set("parameters", str(self.theta))
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate Maxwellian distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Maxwell
            Maxwellian distribution generated from XML element

        """
        theta = float(get_text(elem, 'parameters'))
        return cls(theta)


class Watt(Univariate):
    r"""Watt fission energy spectrum.

    The Watt fission energy spectrum is characterized by two parameters
    :math:`a` and :math:`b` and has density function :math:`p(E) dE = c e^{-E/a}
    \sinh \sqrt{b \, E} dE`.

    Parameters
    ----------
    a : float
        First parameter of distribution in units of eV
    b : float
        Second parameter of distribution in units of 1/eV

    Attributes
    ----------
    a : float
        First parameter of distribution in units of eV
    b : float
        Second parameter of distribution in units of 1/eV

    """

    def __init__(self, a=0.988e6, b=2.249e-6):
        super().__init__()
        self.a = a
        self.b = b

    def __len__(self):
        return 2

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

    def to_xml_element(self, element_name):
        """Return XML representation of the Watt distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing Watt distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "watt")
        element.set("parameters", '{} {}'.format(self.a, self.b))
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate Watt distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Watt
            Watt distribution generated from XML element

        """
        params = get_text(elem, 'parameters').split()
        return cls(*map(float, params))


class Normal(Univariate):
    r"""Normally distributed sampling.

    The Normal Distribution is characterized by two parameters
    :math:`\mu` and :math:`\sigma` and has density function
    :math:`p(X) dX = 1/(\sqrt{2\pi}\sigma) e^{-(X-\mu)^2/(2\sigma^2)}`

    Parameters
    ----------
    mean_value : float
        Mean value of the  distribution
    std_dev : float
        Standard deviation of the Normal distribution

    Attributes
    ----------
    mean_value : float
        Mean of the Normal distribution
    std_dev : float
        Standard deviation of the Normal distribution
    """

    def __init__(self, mean_value, std_dev):
        super().__init__()
        self.mean_value = mean_value
        self.std_dev = std_dev

    def __len__(self):
        return 2

    @property
    def mean_value(self):
        return self._mean_value

    @property
    def std_dev(self):
        return self._std_dev

    @mean_value.setter
    def mean_value(self, mean_value):
        cv.check_type('Normal mean_value', mean_value, Real)
        self._mean_value = mean_value

    @std_dev.setter
    def std_dev(self, std_dev):
        cv.check_type('Normal std_dev', std_dev, Real)
        cv.check_greater_than('Normal std_dev', std_dev, 0.0)
        self._std_dev = std_dev

    def to_xml_element(self, element_name):
        """Return XML representation of the Normal distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing Watt distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "normal")
        element.set("parameters", '{} {}'.format(self.mean_value, self.std_dev))
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate Normal distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Normal
            Normal distribution generated from XML element

        """
        params = get_text(elem, 'parameters').split()
        return cls(*map(float, params))


class Muir(Univariate):
    """Muir energy spectrum.

    The Muir energy spectrum is a Gaussian spectrum, but for
    convenience reasons allows the user 3 parameters to define
    the distribution, e0 the mean energy of particles, the mass
    of reactants m_rat, and the ion temperature kt.

    Parameters
    ----------
    e0 : float
        Mean of the Muir distribution in units of eV
    m_rat : float
        Ratio of the sum of the masses of the reaction inputs to an
        AMU
    kt : float
         Ion temperature for the Muir distribution in units of eV

    Attributes
    ----------
    e0 : float
        Mean of the Muir distribution in units of eV
    m_rat : float
        Ratio of the sum of the masses of the reaction inputs to an
        AMU
    kt : float
         Ion temperature for the Muir distribution in units of eV

    """

    def __init__(self, e0=14.08e6, m_rat = 5., kt = 20000.):
        super().__init__()
        self.e0 = e0
        self.m_rat = m_rat
        self.kt = kt

    def __len__(self):
        return 3

    @property
    def e0(self):
        return self._e0

    @property
    def m_rat(self):
        return self._m_rat

    @property
    def kt(self):
        return self._kt

    @e0.setter
    def e0(self, e0):
        cv.check_type('Muir e0', e0, Real)
        cv.check_greater_than('Muir e0', e0, 0.0)
        self._e0 = e0

    @m_rat.setter
    def m_rat(self, m_rat):
        cv.check_type('Muir m_rat', m_rat, Real)
        cv.check_greater_than('Muir m_rat', m_rat, 0.0)
        self._m_rat = m_rat

    @kt.setter
    def kt(self, kt):
        cv.check_type('Muir kt', kt, Real)
        cv.check_greater_than('Muir kt', kt, 0.0)
        self._kt = kt

    def to_xml_element(self, element_name):
        """Return XML representation of the Watt distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing Watt distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "muir")
        element.set("parameters", '{} {} {}'.format(self._e0, self._m_rat, self._kt))
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate Muir distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Muir
            Muir distribution generated from XML element

        """
        params = get_text(elem, 'parameters').split()
        return cls(*map(float, params))


class Tabular(Univariate):
    """Piecewise continuous probability distribution.

    This class is used to represent a probability distribution whose density
    function is tabulated at specific values with a specified interpolation
    scheme.

    Parameters
    ----------
    x : Iterable of float
        Tabulated values of the random variable
    p : Iterable of float
        Tabulated probabilities
    interpolation : {'histogram', 'linear-linear', 'linear-log', 'log-linear', 'log-log'}, optional
        Indicate whether the density function is constant between tabulated
        points or linearly-interpolated. Defaults to 'linear-linear'.
    ignore_negative : bool
        Ignore negative probabilities

    Attributes
    ----------
    x : Iterable of float
        Tabulated values of the random variable
    p : Iterable of float
        Tabulated probabilities
    interpolation : {'histogram', 'linear-linear', 'linear-log', 'log-linear', 'log-log'}, optional
        Indicate whether the density function is constant between tabulated
        points or linearly-interpolated.

    """

    def __init__(self, x, p, interpolation='linear-linear',
                 ignore_negative=False):
        super().__init__()
        self._ignore_negative = ignore_negative
        self.x = x
        self.p = p
        self.interpolation = interpolation

    def __len__(self):
        return len(self.x)

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
        if not self._ignore_negative:
            for pk in p:
                cv.check_greater_than('tabulated probability', pk, 0.0, True)
        self._p = p

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_value('interpolation', interpolation, _INTERPOLATION_SCHEMES)
        self._interpolation = interpolation

    def to_xml_element(self, element_name):
        """Return XML representation of the tabular distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing tabular distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "tabular")
        element.set("interpolation", self.interpolation)

        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.x)) + ' ' + ' '.join(map(str, self.p))

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate tabular distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Tabular
            Tabular distribution generated from XML element

        """
        interpolation = get_text(elem, 'interpolation')
        params = [float(x) for x in get_text(elem, 'parameters').split()]
        x = params[:len(params)//2]
        p = params[len(params)//2:]
        return cls(x, p, interpolation)


class Legendre(Univariate):
    r"""Probability density given by a Legendre polynomial expansion
    :math:`\sum\limits_{\ell=0}^N \frac{2\ell + 1}{2} a_\ell P_\ell(\mu)`.

    Parameters
    ----------
    coefficients : Iterable of Real
        Expansion coefficients :math:`a_\ell`. Note that the :math:`(2\ell +
        1)/2` factor should not be included.

    Attributes
    ----------
    coefficients : Iterable of Real
        Expansion coefficients :math:`a_\ell`. Note that the :math:`(2\ell +
        1)/2` factor should not be included.

    """

    def __init__(self, coefficients):
        self.coefficients = coefficients
        self._legendre_poly = None

    def __call__(self, x):
        # Create Legendre polynomial if we haven't yet
        if self._legendre_poly is None:
            l = np.arange(len(self._coefficients))
            coeffs = (2.*l + 1.)/2. * self._coefficients
            self._legendre_poly = np.polynomial.Legendre(coeffs)

        return self._legendre_poly(x)

    def __len__(self):
        return len(self._coefficients)

    @property
    def coefficients(self):
        return self._coefficients

    @coefficients.setter
    def coefficients(self, coefficients):
        self._coefficients = np.asarray(coefficients)

    def to_xml_element(self, element_name):
        raise NotImplementedError

    @classmethod
    def from_xml_element(cls, elem):
        raise NotImplementedError


class Mixture(Univariate):
    """Probability distribution characterized by a mixture of random variables.

    Parameters
    ----------
    probability : Iterable of Real
        Probability of selecting a particular distribution
    distribution : Iterable of Univariate
        List of distributions with corresponding probabilities

    Attributes
    ----------
    probability : Iterable of Real
        Probability of selecting a particular distribution
    distribution : Iterable of Univariate
        List of distributions with corresponding probabilities

    """

    def __init__(self, probability, distribution):
        super().__init__()
        self.probability = probability
        self.distribution = distribution

    def __len__(self):
        return sum(len(d) for d in self.distribution)

    @property
    def probability(self):
        return self._probability

    @property
    def distribution(self):
        return self._distribution

    @probability.setter
    def probability(self, probability):
        cv.check_type('mixture distribution probabilities', probability,
                      Iterable, Real)
        for p in probability:
            cv.check_greater_than('mixture distribution probabilities',
                                  p, 0.0, True)
        self._probability = probability

    @distribution.setter
    def distribution(self, distribution):
        cv.check_type('mixture distribution components', distribution,
                      Iterable, Univariate)
        self._distribution = distribution

    def to_xml_element(self, element_name):
        raise NotImplementedError

    @classmethod
    def from_xml_element(cls, elem):
        raise NotImplementedError
