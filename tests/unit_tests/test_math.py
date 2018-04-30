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
    pass


def test_evaluate_legendre():
    max_order = 10
    # Coefficients are set to 1, but will incorporate the (2l+1)/2 norm factor
    # for the reference solution
    test_coeffs = [0.5 * (2. * l + 1.) for l in range(max_order + 1)]
    test_xs = np.linspace(-1., 1., num=5, endpoint=True)

    ref_vals = np.polynomial.legendre.legval(test_xs, test_coeffs)

    # Set the coefficients back to 1s for the test values
    test_coeffs = [1. for l in range(max_order + 1)]
    test_vals = np.array([openmc.capi.math.evaluate_legendre(test_coeffs, x)
                          for x in test_xs])

    assert np.allclose(ref_vals, test_vals)


def test_rotate_angle():
    uvw0 = np.array([1., 0., 0.])
    phi = 0.
    mu = 0.

    # reference: mu of 0 pulls the vector the bottom, so:
    ref_uvw = np.array([0., 0., -1.])

    test_uvw = openmc.capi.math.rotate_angle(uvw0, mu, phi)

    assert np.allclose(ref_uvw, test_uvw)

    # Repeat for mu = 1 (no change)
    mu = 1.
    ref_uvw = np.array([1., 0., 0.])

    test_uvw = openmc.capi.math.rotate_angle(uvw0, mu, phi)

    assert np.allclose(ref_uvw, test_uvw)

    # Need to test phi=None somehow...


def test_maxwell_spectrum():
    settings = openmc.capi.settings
    settings.seed = 1
    T = 0.5
    ref_val = 0.6129982175261098
    test_val = openmc.capi.math.maxwell_spectrum(T)
    print(test_val)
    assert np.isclose(ref_val, test_val)


def test_watt_spectrum():
    settings = openmc.capi.settings
    settings.seed = 1
    a = 0.5
    b = 0.75
    ref_val = 0.6247242713640233
    test_val = openmc.capi.math.watt_spectrum(a, b)
    print(test_val)
    assert np.isclose(ref_val, test_val)


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
