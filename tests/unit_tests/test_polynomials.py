import numpy as np

import openmc


def test_zernike_radial():
    coeff = np.asarray([1.3, -3.0, 9e-1, -6e-1, 0.11])
    zn_rad = openmc.ZernikeRadial(coeff)
    assert zn_rad.order == 8
    assert zn_rad.radius == 1

    coeff = np.asarray([1.3, -3.0, 9e-1, -6e-1, 0.11, 0.222])
    zn_rad = openmc.ZernikeRadial(coeff, 0.392)
    assert zn_rad.order == 10
    assert zn_rad.radius == 0.392
    norm_vec = (2 * np.arange(6) + 1) / (np.pi * 0.392 ** 2)
    norm_coeff = norm_vec * coeff

    rho = 0.5
    # Reference solution from running the Fortran implementation
    raw_zn = np.array([
        1.00000000e+00, -5.00000000e-01, -1.25000000e-01,
        4.37500000e-01, -2.89062500e-01, -8.98437500e-02])

    ref_vals = np.sum(norm_coeff * raw_zn)

    test_vals = zn_rad(rho)

    assert ref_vals == test_vals
