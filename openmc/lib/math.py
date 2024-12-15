from ctypes import c_int, c_double

import numpy as np
from numpy.ctypeslib import ndpointer

from . import _dll


_dll.calc_zn.restype = None
_dll.calc_zn.argtypes = [c_int, c_double, c_double, ndpointer(c_double)]

_dll.calc_zn_rad.restype = None
_dll.calc_zn_rad.argtypes = [c_int, c_double, ndpointer(c_double)]


def calc_zn(n, rho, phi):
    """Calculate the n-th order modified Zernike polynomial moment for a
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
    _dll.calc_zn(n, rho, phi, zn)
    return zn


def calc_zn_rad(n, rho):
    """Calculate the even orders in n-th order modified Zernike polynomial
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
    _dll.calc_zn_rad(n, rho, zn_rad)
    return zn_rad
