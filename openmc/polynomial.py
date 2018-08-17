import numpy as np
import openmc
import openmc.capi as capi


def legendre_from_expcoef(coef, domain):
    """Download file from a URL

    Parameters
    ----------
    coef : Iterable of float
        A list of coefficients of each term in Legendre polynomials

    Returns
    -------
    Legendre class
        A numpy Legendre series object

    """
    n = np.arange(len(coef) + 1)
    c = (2*n + 1) * np.asarray(coef) / (domain[1] - domain[0])
    return np.polynomial.Legendre(c, domain)


class Polynomial(object):
    """Abstract Polynomial Class for creating polynomials.
    """
    def __init__(self, coef):
        self.coef = np.asarray(coef)


class ZernikeRadial(Polynomial):
    """Create radial only Zernike polynomials given coefficients and domian.

    The radial only Zernike polynomials are defined as in
    :class:`ZernikeRadialFilter`.

    Parameters
    ----------
    coef : Iterable of float
        A list of coefficients of each term in radial only Zernike polynomials
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.
    r : float
        Position to be evaluated, normalized on radius [0,1]

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
        self._order = 2 * (self.n_coef - 1)
        self.radius = radius
        norm_vec = (2 * np.arange(self.n_coef) + 1) / (np.pi * radius**2)
        self._norm_coef = norm_vec * self.coef

    @property
    def order(self):
        return self._order

    def __call__(self, r):
        zn_rad = capi.calc_zn_rad(self.order, r)
        return np.sum(self._norm_coef * zn_rad)
