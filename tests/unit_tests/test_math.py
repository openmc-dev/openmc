import numpy as np
import scipy as sp

import openmc
import openmc.capi


# def test_normal_percentile():
#     # normal_percentile has three branches to consider:
#     # p < 0.02425; 0.02425 <= p <= 0.97575; and p > 0.97575
#     test_ps = [0.02, 0.4, 0.5, 0.6, 0.98]

#     # The reference solutions come from Scipy
#     ref_zs = [sp.stats.norm.ppf(p) for p in test_ps]

#     test_zs = [openmc.capi.math.normal_percentile(p) for p in test_ps]

#     assert np.allclose(ref_zs, test_zs)


def test_t_percentile():
    # Permutations include 1 DoF, 2 DoF, and > 2 DoF
    # We will test 5 p-values at 3-DoF values
    test_ps = [0.02, 0.4, 0.5, 0.6, 0.98]
    test_dfs = [1, 2, 5]

    # The reference solutions come from Scipy
    ref_ts = [[sp.stats.t.ppf(p, df) for p in test_ps] for df in test_dfs]

    test_ts = [[openmc.capi.math.t_percentile(p, df) for p in test_ps]
               for df in test_dfs]

    # The 5 DoF approximation in openmc.capi.math.t_percentile is off by up to
    # 8e-3 from the scipy solution, so test that one separately with looser
    # tolerance
    assert np.allclose(ref_ts[:-1], test_ts[:-1])
    assert np.allclose(ref_ts[-1], test_ts[-1], atol=1e-2)


def test_calc_pn():
    max_order = 10
    test_ns = np.array([i for i in range(0, max_order + 1)])
    test_xs = np.linspace(-1., 1., num=5, endpoint=True)

    # Reference solutions from scipy
    ref_vals = [sp.special.eval_legendre(n, test_xs) for n in test_ns]

    test_vals = [[openmc.capi.math.calc_pn(n, x) for x in test_xs]
                 for n in test_ns]

    assert np.allclose(ref_vals, test_vals)


def test_evaluate_legendre():
    max_order = 10
    # Coefficients are set to 1, but will incorporate the (2l+1)/2 norm factor
    # for the reference solution
    test_coeffs = [0.5 * (2. * l + 1.) for l in range(max_order + 1)]
    test_xs = np.linspace(-1., 1., num=5, endpoint=True)

    ref_vals = np.polynomial.legendre.legval(test_xs, test_coeffs)

    # Set the coefficients back to 1s for the test values since
    # evaluate legendre includes the (2l+1)/2 term
    test_coeffs = [1. for l in range(max_order + 1)]

    test_vals = np.array([openmc.capi.math.evaluate_legendre(test_coeffs, x)
                          for x in test_xs])

    assert np.allclose(ref_vals, test_vals)


def test_calc_rn():
    max_order = 10
    test_ns = np.array([i for i in range(0, max_order + 1)])
    azi = 0.1  # Longitude
    pol = 0.2  # Latitude
    test_uvw = np.array([np.sin(pol) * np.cos(azi),
                         np.sin(pol) * np.sin(azi),
                         np.cos(pol)])

    # Reference solutions from the equations
    ref_vals = []

    def coeff(n, m):
        return np.sqrt((2. * n + 1) * sp.special.factorial(n - m) /
                       (sp.special.factorial(n + m)))

    def pnm_bar(n, m, mu):
        val = coeff(n, m)
        if m != 0:
            val *= np.sqrt(2.)
        val *= sp.special.lpmv([m], [n], [mu])
        return val[0]

    ref_vals = []
    for n in test_ns:
        for m in range(-n, n + 1):
            if m < 0:
                ylm = pnm_bar(n, np.abs(m), np.cos(pol)) * \
                    np.sin(np.abs(m) * azi)
            else:
                ylm = pnm_bar(n, m, np.cos(pol)) * np.cos(m * azi)

            # Un-normalize for comparison
            ylm /= np.sqrt(2. * n + 1.)
            ref_vals.append(ylm)

    test_vals = []
    for n in test_ns:
        ylms = openmc.capi.math.calc_rn(n, test_uvw)
        test_vals.extend(ylms.tolist())

    assert np.allclose(ref_vals, test_vals)


def test_calc_zn():
    n = 10
    rho = 0.5
    phi = 0.5

    # Reference solution from running the Fortran implementation
    ref_vals = np.array(
        [1.00000000e+00, 4.79425539e-01, 8.77582562e-01,
         5.15293637e-01, -8.66025404e-01, 3.30866239e-01,
         3.52667735e-01, -8.47512624e-01, -1.55136145e+00,
         2.50093775e-02, 1.79715684e-01, -1.33048245e+00,
         -2.79508497e-01, -8.54292956e-01, -8.22482403e-02,
         6.47865100e-02, -1.18780200e+00, 5.18993370e-01,
         9.50010991e-01, -8.42327938e-02, -8.67263404e-02,
         8.25035501e-03, -7.44248626e-01, 1.52505281e+00,
         1.15751620e+00, 9.79225149e-01, 3.40611006e-01,
         -5.78783240e-02, -1.09619759e-02, -3.17938327e-01,
         1.90147482e+00, 2.84658914e-01, 5.21064646e-01,
         1.34842791e-01, 4.25607546e-01, -2.92642715e-02,
         -1.25423479e-02, -4.67751162e-02, 1.50696182e+00,
         -6.13603897e-01, -8.67187500e-01, -3.93990531e-01,
         -6.89672461e-01, 3.28139254e-01, -1.08327149e-02,
         -8.53837419e-03, 7.04712042e-02, 7.73660979e-01,
         -1.63799891e+00, -8.12396290e-01, -1.48708143e+00,
         -1.16158437e-01, -1.03565982e+00, 1.88131088e-01,
         -1.84122553e-03, -4.39233743e-03, 9.01295675e-02,
         1.32511582e-01, -1.83260987e+00, -7.93994967e-01,
         -2.97978009e-01, -5.09818305e-01, 8.38707753e-01,
         -9.29602211e-01, 7.78441102e-02, 1.29931014e-03])

    test_vals = openmc.capi.math.calc_zn(n, rho, phi)

    assert np.allclose(ref_vals, test_vals)


def test_rotate_angle():
    uvw0 = np.array([1., 0., 0.])
    phi = 0.
    mu = 0.

    # reference: mu of 0 pulls the vector the bottom, so:
    ref_uvw = np.array([0., 0., -1.])

    test_uvw = openmc.capi.math.rotate_angle(uvw0, mu, phi)

    assert np.array_equal(ref_uvw, test_uvw)

    # Repeat for mu = 1 (no change)
    mu = 1.
    ref_uvw = np.array([1., 0., 0.])

    test_uvw = openmc.capi.math.rotate_angle(uvw0, mu, phi)

    assert np.array_equal(ref_uvw, test_uvw)

    # Need to test phi=None somehow...


def test_maxwell_spectrum():
    settings = openmc.capi.settings
    settings.seed = 1
    T = 0.5
    ref_val = 0.6129982175261098
    test_val = openmc.capi.math.maxwell_spectrum(T)

    assert ref_val == test_val


def test_watt_spectrum():
    settings = openmc.capi.settings
    settings.seed = 1
    a = 0.5
    b = 0.75
    ref_val = 0.6247242713640233
    test_val = openmc.capi.math.watt_spectrum(a, b)

    assert ref_val == test_val


def test_broaden_wmp_polynomials():
    # Two branches of the code to worry about, beta > 6 and otherwise
    # beta = sqrtE * dopp
    # First lets do beta > 6
    test_E = 0.5
    test_dopp = 100.  # approximately U235 at room temperature
    n = 4
    ref_val = [2., 1.41421356, 1.0001, 0.70731891]
    test_val = openmc.capi.math.broaden_wmp_polynomials(test_E, test_dopp, n)

    assert np.allclose(ref_val, test_val)

    # now beta < 6
    test_dopp = 5.
    ref_val = [1.99999885, 1.41421356, 1.04, 0.79195959]
    test_val = openmc.capi.math.broaden_wmp_polynomials(test_E, test_dopp, n)

    assert np.allclose(ref_val, test_val)

test_calc_zn()