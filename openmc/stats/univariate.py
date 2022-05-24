from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Iterable
from numbers import Real
from xml.etree import ElementTree as ET

import numpy as np

import openmc.checkvalue as cv
from .._xml import get_text
from ..mixin import EqualityMixin


_INTERPOLATION_SCHEMES = [
    'histogram',
    'linear-linear',
    'linear-log',
    'log-linear',
    'log-log'
]


class Univariate(EqualityMixin, ABC):
    """Probability distribution of a single random variable.

    The Univariate class is an abstract class that can be derived to implement a
    specific probability distribution.

    """
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
        elif distribution == 'powerlaw':
            return PowerLaw.from_xml_element(elem)
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

    @abstractmethod
    def sample(n_samples=1, seed=None):
        """Sample the univariate distribution

        Parameters
        ----------
        n_samples : int
            Number of sampled values to generate
        seed : int or None
            Initial random number seed.

        Returns
        -------
        numpy.ndarray
            A 1-D array of sampled values
        """
        pass


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

    def cdf(self):
        return np.insert(np.cumsum(self.p), 0, 0.0)

    def sample(self, n_samples=1, seed=None):
        np.random.seed(seed)
        return np.random.choice(self.x, n_samples, p=self.p)

    def normalize(self):
        """Normalize the probabilities stored on the distribution"""
        norm = sum(self.p)
        self.p = [val / norm for val in self.p]

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

    @classmethod
    def merge(cls, dists, probs):
        """Merge multiple discrete distributions into a single distribution

        Parameters
        ----------
        dists : iterable of openmc.stats.Discrete
            Discrete distributions to combine
        probs : iterable of float
            Probability of each distribution

        Returns
        -------
        openmc.stats.Discrete
            Combined discrete distribution

        """
        if len(dists) != len(probs):
            raise ValueError("Number of distributions and probabilities must match.")

        # Combine distributions accounting for duplicate x values
        x_merged = set()
        p_merged = defaultdict(float)
        for dist, p_dist in zip(dists, probs):
            for x, p in zip(dist.x, dist.p):
                x_merged.add(x)
                p_merged[x] += p*p_dist

        # Create values and probabilities as arrays
        x_arr = np.array(sorted(x_merged))
        p_arr = np.array([p_merged[x] for x in x_arr])
        return cls(x_arr, p_arr)

    def integral(self):
        """Return integral of distribution

        Returns
        -------
        float
            Integral of discrete distribution
        """
        return np.sum(self.p)

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

    def sample(self, n_samples=1, seed=None):
        np.random.seed(seed)
        return np.random.uniform(self.a, self.b, n_samples)

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


class PowerLaw(Univariate):
    """Distribution with power law probability over a finite interval [a,b]

    The power law distribution has density function :math:`p(x) dx = c x^n dx`.

    .. versionadded:: 0.13.0

    Parameters
    ----------
    a : float, optional
        Lower bound of the sampling interval. Defaults to zero.
    b : float, optional
        Upper bound of the sampling interval. Defaults to unity.
    n : float, optional
        Power law exponent. Defaults to zero, which is equivalent to a uniform
        distribution.

    Attributes
    ----------
    a : float
        Lower bound of the sampling interval
    b : float
        Upper bound of the sampling interval
    n : float
        Power law exponent

    """

    def __init__(self, a=0.0, b=1.0, n=0):
        self.a = a
        self.b = b
        self.n = n

    def __len__(self):
        return 3

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def n(self):
        return self._n

    @a.setter
    def a(self, a):
        cv.check_type('interval lower bound', a, Real)
        self._a = a

    @b.setter
    def b(self, b):
        cv.check_type('interval upper bound', b, Real)
        self._b = b

    @n.setter
    def n(self, n):
        cv.check_type('power law exponent', n, Real)
        self._n = n

    def sample(self, n_samples=1, seed=None):
        np.random.seed(seed)
        xi = np.random.rand(n_samples)
        pwr = self.n + 1
        offset = self.a**pwr
        span = self.b**pwr - offset
        return np.power(offset + xi * span, 1/pwr)

    def to_xml_element(self, element_name):
        """Return XML representation of the power law distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "powerlaw")
        element.set("parameters", f'{self.a} {self.b} {self.n}')
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate power law distribution from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.PowerLaw
            Distribution generated from XML element

        """
        params = get_text(elem, 'parameters').split()
        return cls(*map(float, params))


class Maxwell(Univariate):
    r"""Maxwellian distribution in energy.

    The Maxwellian distribution in energy is characterized by a single parameter
    :math:`\theta` and has a density function :math:`p(E) dE = c \sqrt{E}
    e^{-E/\theta} dE`.

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

    def sample(self, n_samples=1, seed=None):
        np.random.seed(seed)
        return self.sample_maxwell(self.theta, n_samples)

    @staticmethod
    def sample_maxwell(t, n_samples):
        r1, r2, r3 = np.random.rand(3, n_samples)
        c = np.cos(0.5 * np.pi * r3)
        return -t * (np.log(r1) + np.log(r2) * c * c)

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

    def sample(self, n_samples=1, seed=None):
        np.random.seed(seed)
        w = Maxwell.sample_maxwell(self.a, n_samples)
        u = np.random.uniform(-1., 1., n_samples)
        aab = self.a * self.a * self.b
        return w + 0.25*aab + u*np.sqrt(aab*w)

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

    def sample(self, n_samples=1, seed=None):
        np.random.seed(seed)
        return np.random.normal(self.mean_value, self.std_dev, n_samples)

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

    @property
    def std_dev(self):
        return np.sqrt(4.*self.e0*self.kt/self.m_rat)

    def sample(self, n_samples=1, seed=None):
        # Based on LANL report LA-05411-MS
        np.random.seed(seed)
        return np.random.normal(self.e0, self.std_dev, n_samples)

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

    def cdf(self):
        if not self.interpolation in ('histogram', 'linear-linear'):
            raise NotImplementedError('Can only generate CDFs for tabular '
                                      'distributions using histogram or '
                                      'linear-linear interpolation')
        c = np.zeros_like(self.x)
        x = np.asarray(self.x)
        p = np.asarray(self.p)

        if self.interpolation == 'histogram':
            c[1:] = p[:-1] * np.diff(x)
        elif self.interpolation == 'linear-linear':
            c[1:] = 0.5 * (p[:-1] + p[1:]) * np.diff(x)

        return np.cumsum(c)

    def mean(self):
        """Compute the mean of the tabular distribution"""
        if not self.interpolation in ('histogram', 'linear-linear'):
            raise NotImplementedError('Can only compute mean for tabular '
                                      'distributions using histogram '
                                      'or linear-linear interpolation.')
        if self.interpolation == 'linear-linear':
            mean = 0.0
            self.normalize()
            for i in range(1, len(self.x)):
                y_min = self.p[i-1]
                y_max = self.p[i]
                x_min = self.x[i-1]
                x_max = self.x[i]

                m = (y_max - y_min) / (x_max - x_min)

                exp_val = (1./3.) * m * (x_max**3 - x_min**3)
                exp_val += 0.5 * m * x_min * (x_min**2 - x_max**2)
                exp_val += 0.5 * y_min * (x_max**2 - x_min**2)
                mean += exp_val

        elif self.interpolation == 'histogram':
            mean = 0.5 * (self.x[:-1] + self.x[1:])
            mean *= np.diff(self.cdf())
            mean = sum(mean)

        return mean

    def normalize(self):
        """Normalize the probabilities stored on the distribution"""
        self.p = np.asarray(self.p) / self.cdf().max()

    def sample(self, n_samples=1, seed=None):
        if not self.interpolation in ('histogram', 'linear-linear'):
            raise NotImplementedError('Can only sample tabular distributions '
                                      'using histogram or '
                                      'linear-linear interpolation')
        np.random.seed(seed)
        xi = np.random.rand(n_samples)
        cdf = self.cdf()
        cdf /= cdf.max()
        # always use normalized probabilities when sampling
        p = self.p / cdf.max()

        # get CDF bins that are above the
        # sampled values
        c_i = np.full(n_samples, cdf[0])
        cdf_idx = np.zeros(n_samples, dtype=int)
        for i, val in enumerate(cdf[:-1]):
            mask = xi > val
            c_i[mask] = val
            cdf_idx[mask] = i

        # get table values at each index where
        # the random number is less than the next cdf
        # entry
        x_i = self.x[cdf_idx]
        p_i = self.p[cdf_idx]

        if self.interpolation == 'histogram':
            # mask where probability is greater than zero
            pos_mask = p_i > 0.0
            # probabilities greater than zero are set proportional to the
            # position of the random numebers in relation to the cdf value
            p_i[pos_mask] = x_i[pos_mask] + (xi[pos_mask] - c_i[pos_mask]) \
                           / p_i[pos_mask]
            # probabilities smaller than zero are set to the random number value
            p_i[~pos_mask] = x_i[~pos_mask]

            samples_out = p_i

        elif self.interpolation == 'linear-linear':
            # get variable and probability values for the
            # next entry
            x_i1 = self.x[cdf_idx + 1]
            p_i1 = self.p[cdf_idx + 1]
            # compute slope between entries
            m = (p_i1 - p_i) / (x_i1 - x_i)
            # set values for zero slope
            zero = m == 0.0
            m[zero] = x_i[zero] + (xi[zero] - c_i[zero]) / p_i[zero]
            # set values for non-zero slope
            non_zero = ~zero
            quad = np.power(p_i[non_zero], 2) + 2.0 * m[non_zero] * (xi[non_zero] - c_i[non_zero])
            quad[quad < 0.0] = 0.0
            m[non_zero] = x_i[non_zero] + (np.sqrt(quad) - p_i[non_zero]) / m[non_zero]
            samples_out = m

        assert all(samples_out < self.x[-1])
        return samples_out

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

    def integral(self):
        """Return integral of distribution

        Returns
        -------
        float
            Integral of tabular distrbution
        """
        if self.interpolation == 'histogram':
            return np.sum(np.diff(self.x) * self.p[:-1])
        elif self.interpolation == 'linear-linear':
            return np.trapz(self.p, self.x)
        else:
            raise NotImplementedError(
                f'integral() not supported for {self.inteprolation} interpolation')


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

    def sample(self, n_samples=1, seed=None):
        raise NotImplementedError

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

    def cdf(self):
        return np.insert(np.cumsum(self.probability), 0, 0.0)

    def sample(self, n_samples=1, seed=None):
        np.random.seed(seed)
        idx = np.random.choice(self.distribution, n_samples, p=self.probability)

        out = np.zeros_like(idx)
        for i in np.unique(idx):
            n_dist_samples = np.count_nonzero(idx == i)
            samples = self.distribution[i].sample(n_dist_samples)
            out[idx == i] = samples
        return out

    def normalize(self):
        """Normalize the probabilities stored on the distribution"""
        norm = sum(self.probability)
        self.probability = [val / norm for val in self.probability]

    def to_xml_element(self, element_name):
        """Return XML representation of the mixture distribution

        .. versionadded:: 0.13.0

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing mixture distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "mixture")

        for p, d in zip(self.probability, self.distribution):
          data = ET.SubElement(element, "pair")
          data.set("probability", str(p))
          data.append(d.to_xml_element("dist"))

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate mixture distribution from an XML element

        .. versionadded:: 0.13.0

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

        Returns
        -------
        openmc.stats.Mixture
            Mixture distribution generated from XML element

        """
        probability = []
        distribution = []
        for pair in elem.findall('pair'):
            probability.append(float(get_text(pair, 'probability')))
            distribution.append(Univariate.from_xml_element(pair.find("dist")))

        return cls(probability, distribution)

    def integral(self):
        """Return integral of the distribution

        Returns
        -------
        float
            Integral of the distribution
        """
        return sum([
            p*dist.integral()
            for p, dist in zip(self.probability, self.distribution)
        ])
