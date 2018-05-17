from ctypes import (c_int, c_double, POINTER)

import numpy as np
from numpy.ctypeslib import ndpointer

from . import _dll


_dll.t_percentile.restype = c_double
_dll.t_percentile.argtypes = [POINTER(c_double), POINTER(c_int)]

_dll.calc_pn.restype = None
_dll.calc_pn.argtypes = [POINTER(c_int), POINTER(c_double),
                         ndpointer(c_double)]

_dll.evaluate_legendre.restype = c_double
_dll.evaluate_legendre.argtypes = [POINTER(c_int), POINTER(c_double),
                                   POINTER(c_double)]

_dll.calc_rn.restype = None
_dll.calc_rn.argtypes = [POINTER(c_int), ndpointer(c_double),
                         ndpointer(c_double)]

_dll.calc_zn.restype = None
_dll.calc_zn.argtypes = [POINTER(c_int), POINTER(c_double), POINTER(c_double),
                         ndpointer(c_double)]

_dll.rotate_angle.restype = None
_dll.rotate_angle.argtypes = [ndpointer(c_double), POINTER(c_double),
                              ndpointer(c_double), POINTER(c_double)]
_dll.maxwell_spectrum.restype = c_double
_dll.maxwell_spectrum.argtypes = [POINTER(c_double)]

_dll.watt_spectrum.restype = c_double
_dll.watt_spectrum.argtypes = [POINTER(c_double), POINTER(c_double)]

_dll.broaden_wmp_polynomials.restype = None
_dll.broaden_wmp_polynomials.argtypes = [POINTER(c_double), POINTER(c_double),
                                         POINTER(c_int), ndpointer(c_double)]


def t_percentile(p, df):
    """ Calculate the percentile of the Student's t distribution with a
    specified probability level and number of degrees of freedom

    Parameters
    ----------
    p : float
        Probability level
    df : int
        Degrees of freedom

    Returns
    -------
    float
        Corresponding t-value

    """

    return _dll.t_percentile(c_double(p), c_int(df))


def calc_pn(n, x):
    """ Calculate the n-th order Legendre polynomial at the value of x.

    Parameters
    ----------
    n : int
        Legendre order
    x : float
        Independent variable to evaluate the Legendre at

    Returns
    -------
    float
        Corresponding Legendre polynomial result

    """

    pnx = np.empty(n + 1, dtype=np.float64)
    _dll.calc_pn(c_int(n), c_double(x), pnx)
    return pnx


def evaluate_legendre(data, x):
    """ Finds the value of f(x) given a set of Legendre coefficients
    and the value of x.

    Parameters
    ----------
    data : iterable of float
        Legendre coefficients
    x : float
        Independent variable to evaluate the Legendre at

    Returns
    -------
    float
        Corresponding Legendre expansion result

    """

    data_arr = np.array(data, dtype=np.float64)
    return _dll.evaluate_legendre(c_int(len(data)),
                                  data_arr.ctypes.data_as(POINTER(c_double)),
                                  c_double(x))


def calc_rn(n, uvw):
    """ Calculate the n-th order real Spherical Harmonics for a given angle;
    all Rn,m values are provided for all n (where -n <= m <= n).

    Parameters
    ----------
    n : int
        Harmonics order
    uvw : iterable of float
        Independent variable to evaluate the Legendre at

    Returns
    -------
    numpy.ndarray
        Corresponding real harmonics value

    """

    num_nm = (n + 1) * (n + 1)
    rn = np.empty(num_nm, dtype=np.float64)
    uvw_arr = np.array(uvw, dtype=np.float64)
    _dll.calc_rn(c_int(n), uvw_arr, rn)
    return rn


def calc_zn(n, rho, phi):
    """ Calculate the n-th order modified Zernike polynomial moment for a
    given angle (rho, theta) location in the unit disk. The normalization of
    the polynomials is such that the integral of Z_pq*Z_pq over the unit disk
    is exactly pi

    Parameters
    ----------
    n : int
        Maximum order
    rho : float
        Radial location in the unit disk
    phi : float
        Theta (radians) location in the unit disk

    Returns
    -------
    numpy.ndarray
        Corresponding resulting list of coefficients

    """

    num_bins = ((n + 1) * (n + 2)) // 2
    zn = np.zeros(num_bins, dtype=np.float64)
    _dll.calc_zn(c_int(n), c_double(rho), c_double(phi), zn)
    return zn


def rotate_angle(uvw0, mu, phi=None):
    """ Rotates direction cosines through a polar angle whose cosine is
    mu and through an azimuthal angle sampled uniformly.

    Parameters
    ----------
    uvw0 : iterable of float
        Original direction cosine
    mu : float
        Polar angle cosine to rotate
    phi : float, optional
        Azimuthal angle; if None, one will be sampled uniformly

    Returns
    -------
    numpy.ndarray
        Rotated direction cosine

    """

    uvw = np.zeros(3, dtype=np.float64)
    uvw0_arr = np.array(uvw0, dtype=np.float64)

    if phi is None:
        _dll.rotate_angle(uvw0_arr, c_double(mu), uvw, None)
    else:
        _dll.rotate_angle(uvw0_arr, c_double(mu), uvw, c_double(phi))

    return uvw


def maxwell_spectrum(T):
    """ Samples an energy from the Maxwell fission distribution based
    on a direct sampling scheme.

    Parameters
    ----------
    T : float
        Spectrum parameter

    Returns
    -------
    float
        Sampled outgoing energy

    """

    return _dll.maxwell_spectrum(c_double(T))


def watt_spectrum(a, b):
    """ Samples an energy from the Watt energy-dependent fission spectrum.

    Parameters
    ----------
    a : float
        Spectrum parameter a
    b : float
        Spectrum parameter b

    Returns
    -------
    float
        Sampled outgoing energy

    """

    return _dll.watt_spectrum(c_double(a), c_double(b))


def broaden_wmp_polynomials(E, dopp, n):
    """ Doppler broadens the windowed multipole curvefit.  The curvefit is a
    polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E) ...

    Parameters
    ----------
    E : float
        Energy to evaluate at
    dopp : float
        sqrt(atomic weight ratio / kT), with kT given in eV
    n : int
        Number of components to the polynomial

    Returns
    -------
    numpy.ndarray
        Resultant leading coefficients

    """

    factors = np.zeros(n, dtype=np.float64)
    _dll.broaden_wmp_polynomials(c_double(E), c_double(dopp), c_int(n),
                                 factors)
    return factors
