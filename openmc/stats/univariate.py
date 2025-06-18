from __future__ import annotations
import math
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Iterable, Sequence
from copy import deepcopy
from numbers import Real
from warnings import warn

import lxml.etree as ET
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import exprel

import openmc.checkvalue as cv
from .._xml import get_text
from ..mixin import EqualityMixin

_INTERPOLATION_SCHEMES = {
    'histogram',
    'linear-linear',
    'linear-log',
    'log-linear',
    'log-log'
}

def log1prel(x):
    """Evaluate log(1+x)/x without loss of precision near 0"""
    return np.where(np.abs(x) < 1e-16, 1.0, np.log1p(x) / x)

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
            # Support older files where Muir had its own class
            params = [float(x) for x in get_text(elem, 'parameters').split()]
            return muir(*params)
        elif distribution == 'tabular':
            return Tabular.from_xml_element(elem)
        elif distribution == 'legendre':
            return Legendre.from_xml_element(elem)
        elif distribution == 'mixture':
            return Mixture.from_xml_element(elem)

    @abstractmethod
    def sample(n_samples: int = 1, seed: int | None = None):
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

    def integral(self):
        """Return integral of distribution

        .. versionadded:: 0.13.1

        Returns
        -------
        float
            Integral of distribution
        """
        return 1.0


def _intensity_clip(intensity: Sequence[float], tolerance: float = 1e-6) -> np.ndarray:
    """Clip low-importance points from an array of intensities.

    Given an array of intensities, this function returns an array of indices for
    points that contribute non-negligibly to the total sum of intensities.

    Parameters
    ----------
    intensity : sequence of float
        Intensities in arbitrary units.
    tolerance : float
        Maximum fraction of intensities that will be discarded.

    Returns
    -------
    Array of indices

    """
    # Get indices of intensities from largest to smallest
    index_sort = np.argsort(intensity)[::-1]

    # Get intensities from largest to smallest
    sorted_intensity = np.asarray(intensity)[index_sort]

    # Determine cumulative sum of probabilities
    cumsum = np.cumsum(sorted_intensity)
    cumsum /= cumsum[-1]

    # Find index that satisfies cutoff
    index_cutoff = np.searchsorted(cumsum, 1.0 - tolerance)

    # Now get indices up to cutoff
    new_indices = index_sort[:index_cutoff + 1]

    # Put back in the order of the original array and return
    new_indices.sort()
    return new_indices


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
    x : numpy.ndarray
        Values of the random variable
    p : numpy.ndarray
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

    @x.setter
    def x(self, x):
        if isinstance(x, Real):
            x = [x]
        cv.check_type('discrete values', x, Iterable, Real)
        self._x = np.array(x, dtype=float)

    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, p):
        if isinstance(p, Real):
            p = [p]
        cv.check_type('discrete probabilities', p, Iterable, Real)
        for pk in p:
            cv.check_greater_than('discrete probability', pk, 0.0, True)
        self._p = np.array(p, dtype=float)

    def cdf(self):
        return np.insert(np.cumsum(self.p), 0, 0.0)

    def sample(self, n_samples=1, seed=None):
        rng = np.random.RandomState(seed)
        p = self.p / self.p.sum()
        return rng.choice(self.x, n_samples, p=p)

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
        element : lxml.etree._Element
            XML element containing discrete distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "discrete")

        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.x)) + ' ' + ' '.join(map(str, self.p))

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate discrete distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
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
    def merge(
        cls,
        dists: Sequence[Discrete],
        probs: Sequence[int]
    ):
        """Merge multiple discrete distributions into a single distribution

        .. versionadded:: 0.13.1

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

        .. versionadded:: 0.13.1

        Returns
        -------
        float
            Integral of discrete distribution
        """
        return np.sum(self.p)

    def clip(self, tolerance: float = 1e-6, inplace: bool = False) -> Discrete:
        r"""Remove low-importance points from discrete distribution.

        Given a probability mass function :math:`p(x)` with :math:`\{x_1, x_2,
        x_3, \dots\}` the possible values of the random variable with
        corresponding probabilities :math:`\{p_1, p_2, p_3, \dots\}`, this
        function will remove any low-importance points such that :math:`\sum_i
        x_i p_i` is preserved to within some threshold.

        .. versionadded:: 0.14.0

        Parameters
        ----------
        tolerance : float
            Maximum fraction of :math:`\sum_i x_i p_i` that will be discarded.
        inplace : bool
            Whether to modify the current object in-place or return a new one.

        Returns
        -------
        Discrete distribution with low-importance points removed

        """
        cv.check_less_than("tolerance", tolerance, 1.0, equality=True)
        cv.check_greater_than("tolerance", tolerance, 0.0, equality=True)

        # Compute intensities
        intensity = self.p * self.x

        # Get indices for intensities above threshold
        indices = _intensity_clip(intensity, tolerance=tolerance)

        # Create new discrete distribution
        if inplace:
            self.x = self.x[indices]
            self.p = self.p[indices]
            return self
        else:
            new_x = self.x[indices]
            new_p = self.p[indices]
            return type(self)(new_x, new_p)


def delta_function(value: float, intensity: float = 1.0) -> Discrete:
    """Return a discrete distribution with a single point.

    .. versionadded:: 0.15.1

    Parameters
    ----------
    value : float
        Value of the random variable.
    intensity : float, optional
        When used for an energy distribution, this can be used to assign an
        intensity.

    Returns
    -------
    Discrete distribution with a single point

    """
    return Discrete([value], [intensity])


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

    def __init__(self, a: float = 0.0, b: float = 1.0):
        self.a = a
        self.b = b

    def __len__(self):
        return 2

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, a):
        cv.check_type('Uniform a', a, Real)
        self._a = a

    @property
    def b(self):
        return self._b

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
        rng = np.random.RandomState(seed)
        return rng.uniform(self.a, self.b, n_samples)

    def to_xml_element(self, element_name: str):
        """Return XML representation of the uniform distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : lxml.etree._Element
            XML element containing uniform distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "uniform")
        element.set("parameters", f'{self.a} {self.b}')
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate uniform distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
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

    def __init__(self, a: float = 0.0, b: float = 1.0, n: float = 0.):
        self.a = a
        self.b = b
        self.n = n

    def __len__(self):
        return 3

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, a):
        cv.check_type('interval lower bound', a, Real)
        self._a = a

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, b):
        cv.check_type('interval upper bound', b, Real)
        self._b = b

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, n):
        cv.check_type('power law exponent', n, Real)
        self._n = n

    def sample(self, n_samples=1, seed=None):
        rng = np.random.RandomState(seed)
        xi = rng.random(n_samples)
        f = np.log(self.b/self.a)*exprel((self.n+1)*np.log(self.b/self.a))
        return self.a*np.exp(f*xi*log1prel(((self.n+1)*f)*xi))

    def to_xml_element(self, element_name: str):
        """Return XML representation of the power law distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : lxml.etree._Element
            XML element containing distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "powerlaw")
        element.set("parameters", f'{self.a} {self.b} {self.n}')
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate power law distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
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
        rng = np.random.RandomState(seed)
        return self.sample_maxwell(self.theta, n_samples, rng=rng)

    @staticmethod
    def sample_maxwell(t, n_samples: int, rng=None):
        if rng is None:
            rng = np.random.default_rng()
        r1, r2, r3 = rng.random((3, n_samples))
        c = np.cos(0.5 * np.pi * r3)
        return -t * (np.log(r1) + np.log(r2) * c * c)

    def to_xml_element(self, element_name: str):
        """Return XML representation of the Maxwellian distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : lxml.etree._Element
            XML element containing Maxwellian distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "maxwell")
        element.set("parameters", str(self.theta))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate Maxwellian distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
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

    @a.setter
    def a(self, a):
        cv.check_type('Watt a', a, Real)
        cv.check_greater_than('Watt a', a, 0.0)
        self._a = a

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, b):
        cv.check_type('Watt b', b, Real)
        cv.check_greater_than('Watt b', b, 0.0)
        self._b = b

    def sample(self, n_samples=1, seed=None):
        rng = np.random.RandomState(seed)
        w = Maxwell.sample_maxwell(self.a, n_samples, rng=rng)
        u = rng.uniform(-1., 1., n_samples)
        aab = self.a * self.a * self.b
        return w + 0.25*aab + u*np.sqrt(aab*w)

    def to_xml_element(self, element_name: str):
        """Return XML representation of the Watt distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : lxml.etree._Element
            XML element containing Watt distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "watt")
        element.set("parameters", f'{self.a} {self.b}')
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate Watt distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
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

    @mean_value.setter
    def mean_value(self, mean_value):
        cv.check_type('Normal mean_value', mean_value, Real)
        self._mean_value = mean_value

    @property
    def std_dev(self):
        return self._std_dev

    @std_dev.setter
    def std_dev(self, std_dev):
        cv.check_type('Normal std_dev', std_dev, Real)
        cv.check_greater_than('Normal std_dev', std_dev, 0.0)
        self._std_dev = std_dev

    def sample(self, n_samples=1, seed=None):
        rng = np.random.RandomState(seed)
        return rng.normal(self.mean_value, self.std_dev, n_samples)

    def to_xml_element(self, element_name: str):
        """Return XML representation of the Normal distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : lxml.etree._Element
            XML element containing Watt distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "normal")
        element.set("parameters", f'{self.mean_value} {self.std_dev}')
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate Normal distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.Normal
            Normal distribution generated from XML element

        """
        params = get_text(elem, 'parameters').split()
        return cls(*map(float, params))


def muir(e0: float, m_rat: float, kt: float):
    """Generate a Muir energy spectrum

    The Muir energy spectrum is a normal distribution, but for convenience
    reasons allows the user to specify three parameters to define the
    distribution: the mean energy of particles ``e0``, the mass of reactants
    ``m_rat``, and the ion temperature ``kt``.

    .. versionadded:: 0.13.2

    Parameters
    ----------
    e0 : float
        Mean of the Muir distribution in [eV]
    m_rat : float
        Ratio of the sum of the masses of the reaction inputs to 1 amu
    kt : float
         Ion temperature for the Muir distribution in [eV]

    Returns
    -------
    openmc.stats.Normal
        Corresponding normal distribution

    """
    # https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-05411-MS
    std_dev = math.sqrt(2 * e0 * kt / m_rat)
    return Normal(e0, std_dev)


# Retain deprecated name for the time being
def Muir(*args, **kwargs):
    # warn of name change
    warn(
        "The Muir(...) class has been replaced by the muir(...) function and "
        "will be removed in a future version of OpenMC. Use muir(...) instead.",
        FutureWarning
    )
    return muir(*args, **kwargs)


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
        Tabulated probabilities. For histogram interpolation, if the length of
        `p` is the same as `x`, the last value is ignored. Probabilities `p` are
        given per unit of `x`.
    interpolation : {'histogram', 'linear-linear', 'linear-log', 'log-linear', 'log-log'}, optional
        Indicates how the density function is interpolated between tabulated
        points. Defaults to 'linear-linear'.
    ignore_negative : bool
        Ignore negative probabilities

    Attributes
    ----------
    x : numpy.ndarray
        Tabulated values of the random variable
    p : numpy.ndarray
        Tabulated probabilities
    interpolation : {'histogram', 'linear-linear', 'linear-log', 'log-linear', 'log-log'}
        Indicates how the density function is interpolated between tabulated
        points. Defaults to 'linear-linear'.

    Notes
    -----
    The probabilities `p` are interpreted per unit of the corresponding
    independent variable `x`. This follows the definition of a probability
    density function (PDF) in probability theory, where the PDF represents the
    relative likelihood of the random variable taking on a particular value per
    unit of the variable. For example, if `x` represents energy in eV, then `p`
    should represent probabilities per eV.

    """

    def __init__(
            self,
            x: Sequence[float],
            p: Sequence[float],
            interpolation: str = 'linear-linear',
            ignore_negative: bool = False
        ):
        self.interpolation = interpolation

        cv.check_type('tabulated values', x, Iterable, Real)
        cv.check_type('tabulated probabilities', p, Iterable, Real)

        x = np.array(x, dtype=float)
        p = np.array(p, dtype=float)

        if p.size > x.size:
            raise ValueError('Number of probabilities exceeds number of table values.')
        if self.interpolation != 'histogram' and x.size != p.size:
            raise ValueError(f'Tabulated values ({x.size}) and probabilities '
                             f'({p.size}) should have the same length')

        if not ignore_negative:
            for pk in p:
                cv.check_greater_than('tabulated probability', pk, 0.0, True)

        self._x = x
        self._p = p

    def __len__(self):
        return self.p.size

    @property
    def x(self):
        return self._x

    @property
    def p(self):
        return self._p

    @property
    def interpolation(self):
        return self._interpolation

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_value('interpolation', interpolation, _INTERPOLATION_SCHEMES)
        self._interpolation = interpolation

    def cdf(self):
        c = np.zeros_like(self.x)
        x = self.x
        p = self.p

        if self.interpolation == 'histogram':
            c[1:] = p[:x.size-1] * np.diff(x)
        elif self.interpolation == 'linear-linear':
            c[1:] = 0.5 * (p[:-1] + p[1:]) * np.diff(x)
        else:
            raise NotImplementedError('Can only generate CDFs for tabular '
                                      'distributions using histogram or '
                                      'linear-linear interpolation')


        return np.cumsum(c)

    def mean(self):
        """Compute the mean of the tabular distribution"""
        if self.interpolation == 'linear-linear':
            mean = 0.0
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
            x_l = self.x[:-1]
            x_r = self.x[1:]
            p_l = self.p[:self.x.size-1]
            mean = (0.5 * (x_l + x_r) * (x_r - x_l) * p_l).sum()
        else:
            raise NotImplementedError('Can only compute mean for tabular '
                                      'distributions using histogram '
                                      'or linear-linear interpolation.')

        # Normalize for when integral of distribution is not 1
        mean /= self.integral()

        return mean

    def normalize(self):
        """Normalize the probabilities stored on the distribution"""
        self._p /= self.cdf().max()

    def sample(self, n_samples: int = 1, seed: int | None = None):
        rng = np.random.RandomState(seed)
        xi = rng.random(n_samples)

        # always use normalized probabilities when sampling
        cdf = self.cdf()
        p = self.p / cdf.max()
        cdf /= cdf.max()

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
        p_i = p[cdf_idx]

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
            p_i1 = p[cdf_idx + 1]
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

        else:
            raise NotImplementedError('Can only sample tabular distributions '
                                      'using histogram or '
                                      'linear-linear interpolation')

        assert all(samples_out < self.x[-1])
        return samples_out

    def to_xml_element(self, element_name: str):
        """Return XML representation of the tabular distribution

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : lxml.etree._Element
            XML element containing tabular distribution data

        """
        element = ET.Element(element_name)
        element.set("type", "tabular")
        element.set("interpolation", self.interpolation)

        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.x)) + ' ' + ' '.join(map(str, self.p))

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate tabular distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.Tabular
            Tabular distribution generated from XML element

        """
        interpolation = get_text(elem, 'interpolation')
        params = [float(x) for x in get_text(elem, 'parameters').split()]
        m = (len(params) + 1)//2  # +1 for when len(params) is odd
        x = params[:m]
        p = params[m:]
        return cls(x, p, interpolation)

    def integral(self):
        """Return integral of distribution

        .. versionadded:: 0.13.1

        Returns
        -------
        float
            Integral of tabular distrbution
        """
        if self.interpolation == 'histogram':
            return np.sum(np.diff(self.x) * self.p[:self.x.size-1])
        elif self.interpolation == 'linear-linear':
            return trapezoid(self.p, self.x)
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

    def __init__(self, coefficients: Sequence[float]):
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

    def __init__(
        self,
        probability: Sequence[float],
        distribution: Sequence[Univariate]
    ):
        self.probability = probability
        self.distribution = distribution

    def __len__(self):
        return sum(len(d) for d in self.distribution)

    @property
    def probability(self):
        return self._probability

    @probability.setter
    def probability(self, probability):
        cv.check_type('mixture distribution probabilities', probability,
                      Iterable, Real)
        for p in probability:
            cv.check_greater_than('mixture distribution probabilities',
                                  p, 0.0, True)
        self._probability = np.array(probability, dtype=float)

    @property
    def distribution(self):
        return self._distribution

    @distribution.setter
    def distribution(self, distribution):
        cv.check_type('mixture distribution components', distribution,
                      Iterable, Univariate)
        self._distribution = distribution

    def cdf(self):
        return np.insert(np.cumsum(self.probability), 0, 0.0)

    def sample(self, n_samples=1, seed=None):
        rng = np.random.RandomState(seed)

        # Get probability of each distribution accounting for its intensity
        p = np.array([prob*dist.integral() for prob, dist in
                      zip(self.probability, self.distribution)])
        p /= p.sum()

        # Sample from the distributions
        idx = rng.choice(range(len(self.distribution)), n_samples, p=p)

        # Draw samples from the distributions sampled above
        out = np.empty_like(idx, dtype=float)
        for i in np.unique(idx):
            n_dist_samples = np.count_nonzero(idx == i)
            samples = self.distribution[i].sample(n_dist_samples)
            out[idx == i] = samples
        return out

    def normalize(self):
        """Normalize the probabilities stored on the distribution"""
        norm = sum(self.probability)
        self.probability = [val / norm for val in self.probability]

    def to_xml_element(self, element_name: str):
        """Return XML representation of the mixture distribution

        .. versionadded:: 0.13.0

        Parameters
        ----------
        element_name : str
            XML element name

        Returns
        -------
        element : lxml.etree._Element
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
    def from_xml_element(cls, elem: ET.Element):
        """Generate mixture distribution from an XML element

        .. versionadded:: 0.13.0

        Parameters
        ----------
        elem : lxml.etree._Element
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

        .. versionadded:: 0.13.1

        Returns
        -------
        float
            Integral of the distribution
        """
        return sum([
            p*dist.integral()
            for p, dist in zip(self.probability, self.distribution)
        ])

    def clip(self, tolerance: float = 1e-6, inplace: bool = False) -> Mixture:
        r"""Remove low-importance points / distributions

        Like :meth:`Discrete.clip`, this method will remove low-importance
        points from discrete distributions contained within the mixture but it
        will also clip any distributions that have negligible contributions to
        the overall intensity.

        .. versionadded:: 0.14.0

        Parameters
        ----------
        tolerance : float
            Maximum fraction of intensities that will be discarded.
        inplace : bool
            Whether to modify the current object in-place or return a new one.

        Returns
        -------
        Distribution with low-importance points / distributions removed

        """
        # Determine integral of original distribution to compare later
        original_integral = self.integral()

        # Determine indices for any distributions that contribute non-negligibly
        # to overall intensity
        intensities = [prob*dist.integral() for prob, dist in
                       zip(self.probability, self.distribution)]
        indices = _intensity_clip(intensities, tolerance=tolerance)

        # Clip mixture of distributions
        probability = self.probability[indices]
        distribution = [self.distribution[i] for i in indices]

        # Clip points from Discrete distributions
        distribution = [
            dist.clip(tolerance, inplace) if isinstance(dist, Discrete) else dist
            for dist in distribution
        ]

        if inplace:
            # Set attributes of current object and return
            self.probability = probability
            self.distribution = distribution
            new_dist = self
        else:
            # Create new distribution
            new_dist = type(self)(probability, distribution)

        # Show warning if integral of new distribution is not within
        # tolerance of original
        diff = (original_integral - new_dist.integral())/original_integral
        if diff > tolerance:
            warn("Clipping mixture distribution resulted in an integral that is "
                    f"lower by a fraction of {diff} when tolerance={tolerance}.")

        return new_dist


def combine_distributions(
    dists: Sequence[Univariate],
    probs: Sequence[float]
):
    """Combine distributions with specified probabilities

    This function can be used to combine multiple instances of
    :class:`~openmc.stats.Discrete` and `~openmc.stats.Tabular`. Multiple
    discrete distributions are merged into a single distribution and the
    remainder of the distributions are put into a :class:`~openmc.stats.Mixture`
    distribution.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    dists : iterable of openmc.stats.Univariate
        Distributions to combine
    probs : iterable of float
        Probability (or intensity) of each distribution

    """
    # Get copy of distribution list so as not to modify the argument
    dist_list = deepcopy(dists)

    # Get list of discrete/continuous distribution indices
    discrete_index = [i for i, d in enumerate(dist_list) if isinstance(d, Discrete)]
    cont_index = [i for i, d in enumerate(dist_list) if isinstance(d, Tabular)]

    # Apply probabilites to continuous distributions
    for i in cont_index:
        dist = dist_list[i]
        dist._p *= probs[i]

    if discrete_index:
        # Create combined discrete distribution
        dist_discrete = [dist_list[i] for i in discrete_index]
        discrete_probs = [probs[i] for i in discrete_index]
        combined_dist = Discrete.merge(dist_discrete, discrete_probs)

        # Replace multiple discrete distributions with merged
        for idx in reversed(discrete_index):
            dist_list.pop(idx)
        dist_list.append(combined_dist)

    # Combine discrete and continuous if present
    if len(dist_list) > 1:
        probs = [1.0]*len(dist_list)
        dist_list[:] = [Mixture(probs, dist_list.copy())]

    return dist_list[0]
