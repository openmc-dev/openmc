from ctypes import (c_int, c_double, POINTER)

import numpy as np
from numpy.ctypeslib import ndpointer

from . import _dll


_dll.t_percentile_c.restype = c_double
_dll.t_percentile_c.argtypes = [c_double, c_int]

_dll.calc_pn_c.restype = None
_dll.calc_pn_c.argtypes = [c_int, c_double, ndpointer(c_double)]

_dll.evaluate_legendre_c.restype = c_double
_dll.evaluate_legendre_c.argtypes = [c_int, POINTER(c_double), c_double]

_dll.calc_rn_c.restype = None
_dll.calc_rn_c.argtypes = [c_int, ndpointer(c_double), ndpointer(c_double)]

_dll.calc_zn_c.restype = None
_dll.calc_zn_c.argtypes = [c_int, c_double, c_double, ndpointer(c_double)]

_dll.calc_zn_rad_c.restype = None
_dll.calc_zn_rad_c.argtypes = [c_int, c_double, ndpointer(c_double)]

_dll.rotate_angle_c.restype = None
_dll.rotate_angle_c.argtypes = [ndpointer(c_double), c_double,
                                POINTER(c_double)]
_dll.maxwell_spectrum_c.restype = c_double
_dll.maxwell_spectrum_c.argtypes = [c_double]

_dll.watt_spectrum_c.restype = c_double
_dll.watt_spectrum_c.argtypes = [c_double, c_double]

_dll.broaden_wmp_polynomials_c.restype = None
_dll.broaden_wmp_polynomials_c.argtypes = [c_double, c_double, c_int,
                                           ndpointer(c_double)]


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

    return _dll.t_percentile_c(p, df)


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
    _dll.calc_pn_c(n, x, pnx)
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
    return _dll.evaluate_legendre_c(len(data),
                                    data_arr.ctypes.data_as(POINTER(c_double)),
                                    x)


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
    _dll.calc_rn_c(n, uvw_arr, rn)
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
    _dll.calc_zn_c(n, rho, phi, zn)
    return zn


def calc_zn_rad(n, rho):
    """ Calculate the even orders in n-th order modified Zernike polynomial
    moment with no azimuthal dependency (m=0) for a given radial location in
    the unit disk. The normalization of the polynomials is such that the
    integral of Z_pq*Z_pq over the unit disk is exactly pi.

    Parameters
    ----------
    n : int
        Maximum order
    rho : float
        Radial location in the unit disk

    Returns
    -------
    numpy.ndarray
        Corresponding resulting list of coefficients

    """

    num_bins = n // 2 + 1
    zn_rad = np.zeros(num_bins, dtype=np.float64)
    _dll.calc_zn_rad_c(n, rho, zn_rad)
    return zn_rad
    

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

    uvw0_arr = np.array(uvw0, dtype=np.float64)

    if phi is None:
        _dll.rotate_angle_c(uvw0_arr, mu, None)
    else:
        _dll.rotate_angle_c(uvw0_arr, mu, c_double(phi))
    uvw = uvw0_arr

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

    return _dll.maxwell_spectrum_c(T)


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

    return _dll.watt_spectrum_c(a, b)


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
    _dll.broaden_wmp_polynomials_c(E, dopp, n, factors)
    return factors
