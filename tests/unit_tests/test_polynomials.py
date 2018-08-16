import numpy as np
import scipy as sp

import openmc
import openmc.capi

def test_legendre():
    coeff = np.asarray([2.5, 1.2, -0.2, 0.1])
    pn = openmc.Legendre(coeff)
    assert pn.order == 3
    assert pn.domain_min == -1
    assert pn.domain_max == 1

    norm_vec = (2 * np.arange(4) + 1) / 2
    assert np.allclose(pn.norm_vec,norm_vec)
    norm_coeff = np.multiply(norm_vec,coeff)
    assert np.allclose(pn.norm_coeff, norm_coeff)

    coeff = np.asarray([2.5, 1.2, -0.2])
    pn = openmc.Legendre(coeff, -10, 10)
    assert pn.order == 2
    assert pn.domain_min == -10
    assert pn.domain_max == 10
    norm_vec = (2 * np.arange(3) + 1) / 20
    assert np.allclose(pn.norm_vec,norm_vec)
    norm_coeff = np.multiply(norm_vec,coeff)
    assert np.allclose(pn.norm_coeff, norm_coeff)

    test_xs = np.linspace(-1., 1., num=5, endpoint=True)
    ref_vals = np.polynomial.legendre.legval(test_xs, norm_coeff)

    test_vals = np.array([pn(x) for x in test_xs])
    assert np.allclose(ref_vals, test_vals)


def test_zernike_radial():
    coeff = np.asarray([1.3, -3.0, 9e-1, -6e-1, 0.11])
    zn_rad = openmc.ZernikeRadial(coeff)
    assert zn_rad.order == 8
    assert zn_rad.radius == 1
    norm_vec = (2 * np.arange(5) + 1) / np.pi
    assert np.allclose(zn_rad.norm_vec,norm_vec)
    norm_coeff = np.multiply(norm_vec,coeff)
    assert np.allclose(zn_rad.norm_coeff,norm_coeff)

    coeff = np.asarray([1.3, -3.0, 9e-1, -6e-1, 0.11, 0.222])
    zn_rad = openmc.ZernikeRadial(coeff, 0.392)
    assert zn_rad.order == 10
    assert zn_rad.radius == 0.392
    norm_vec = (2 * np.arange(6) + 1) / (np.pi * 0.392 ** 2)
    assert np.allclose(zn_rad.norm_vec,norm_vec)
    norm_coeff = np.multiply(norm_vec,coeff)
    assert np.allclose(zn_rad.norm_coeff,norm_coeff)


    rho = 0.5
    # Reference solution from running the Fortran implementation   
    raw_zn = np.array([
        1.00000000e+00, -5.00000000e-01, -1.25000000e-01,
        4.37500000e-01, -2.89062500e-01,-8.98437500e-02])

    ref_vals = np.sum(np.multiply(norm_coeff,raw_zn))

    test_vals = zn_rad(rho)

    assert ref_vals == test_vals

