import math
import numpy as np
import openmc
from collections.abc import Iterable


def legendre_from_expcoef(coef, domain=(-1, 1)):
    """Return a Legendre series object based on expansion coefficients.

    Given a list of coefficients from FET tally and a array of down, return
    the numpy Legendre object.

    Parameters
    ----------
    coef : Iterable of float
        A list of coefficients of each term in Legendre polynomials
    domain : (2,) List of float
        Domain of the Legendre polynomial

    Returns
    -------
    numpy.polynomial.Legendre
        A numpy Legendre series class

    """

    n = np.arange(len(coef))
    c = (2*n + 1) * np.asarray(coef) / (domain[1] - domain[0])
    return np.polynomial.Legendre(c, domain)


class Polynomial:
    """Abstract Polynomial Class for creating polynomials.
    """
    def __init__(self, coef):
        self.coef = np.asarray(coef)


class ZernikeRadial(Polynomial):
    """Create radial only Zernike polynomials given coefficients and domain.

    The radial only Zernike polynomials are defined as in
    :class:`ZernikeRadialFilter`.

    Parameters
    ----------
    coef : Iterable of float
        A list of coefficients of each term in radial only Zernike polynomials
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.

    Attributes
    ----------
    order : int
        The maximum (even) order of Zernike polynomials.
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.
    norm_coef : iterable of float
        The list of coefficients of each term in the polynomials after
        normailization.

    """
    def __init__(self, coef, radius=1):
        super().__init__(coef)
        self._order = 2 * (len(self.coef) - 1)
        self.radius = radius
        norm_vec = (2 * np.arange(len(self.coef)) + 1) / (math.pi * radius**2)
        self._norm_coef = norm_vec * self.coef

    @property
    def order(self):
        return self._order

    def __call__(self, r):
        import openmc.lib as lib
        if isinstance(r, Iterable):
            return [np.sum(self._norm_coef * lib.calc_zn_rad(self.order, r_i / self.radius))
                    for r_i in r]
        else:
            return np.sum(self._norm_coef * lib.calc_zn_rad(self.order, r / self.radius))


class Zernike(Polynomial):
    r"""Create Zernike polynomials given coefficients and domain.

    The azimuthal Zernike polynomials are defined as in :class:`ZernikeFilter`.

    Parameters
    ----------
    coef : Iterable of float
        A list of coefficients of each term in radial only Zernike polynomials
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.

    Attributes
    ----------
    order : int
        The maximum (even) order of Zernike polynomials.
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.
    theta : float
        Azimuthal of Zernike polynomial to be applied on. Default is 0.
    norm_coef : iterable of float
        The list of coefficients of each term in the polynomials after
        normailization.
    """
    def __init__(self, coef, radius=1):
        super().__init__(coef)
        # Solve order from number of coefficients
        # N = (order + 1)(order + 2) / 2
        self._order = int((math.sqrt(8 * len(self.coef) + 1) - 3) / 2)
        self.radius = radius
        norm_vec = np.ones(len(self.coef))
        for n in range(self._order + 1):
            for m in range(-n, n + 1, 2):
                j = int((n*(n + 2) + m)/2)
                if m == 0:
                    norm_vec[j] = n + 1
                else:
                    norm_vec[j] = 2*n + 2
        norm_vec /= (math.pi * radius**2)
        self._norm_coef = norm_vec * self.coef

    @property
    def order(self):
        return self._order

    def __call__(self, r, theta=0.0):
        import openmc.lib as lib
        if isinstance(r, Iterable) and isinstance(theta, Iterable):
            return [[np.sum(self._norm_coef * lib.calc_zn(self.order, r_i / self.radius, theta_i))
                    for r_i in r] for theta_i in theta]
        elif isinstance(r, Iterable) and not isinstance(theta, Iterable):
            return [np.sum(self._norm_coef * lib.calc_zn(self.order, r_i / self.radius, theta))
                    for r_i in r]
        elif not isinstance(r, Iterable) and isinstance(theta, Iterable):
            return [np.sum(self._norm_coef * lib.calc_zn(self.order, r / self.radius, theta_i))
                    for theta_i in theta]
        else:
            return np.sum(self._norm_coef * lib.calc_zn(self.order, r / self.radius, theta))
