import numpy as np
import scipy as sp

import openmc
import openmc.lib

import pytest

def test_t_percentile():
    # Permutations include 1 DoF, 2 DoF, and > 2 DoF
    # We will test 5 p-values at 3-DoF values
    test_ps = [0.02, 0.4, 0.5, 0.6, 0.98]
    test_dfs = [1, 2, 5]

    # The reference solutions come from Scipy
    ref_ts = [[sp.stats.t.ppf(p, df) for p in test_ps] for df in test_dfs]

    test_ts = [[openmc.lib.math.t_percentile(p, df) for p in test_ps]
               for df in test_dfs]

    # The 5 DoF approximation in openmc.lib.math.t_percentile is off by up to
    # 8e-3 from the scipy solution, so test that one separately with looser
    # tolerance
    assert np.allclose(ref_ts[:-1], test_ts[:-1])
    assert np.allclose(ref_ts[-1], test_ts[-1], atol=1e-2)


def test_calc_pn():
    max_order = 10
    test_xs = np.linspace(-1., 1., num=5, endpoint=True)

    # Reference solutions from scipy
    ref_vals = np.array([sp.special.eval_legendre(n, test_xs)
                         for n in range(0, max_order + 1)])

    test_vals = []
    for x in test_xs:
        test_vals.append(openmc.lib.math.calc_pn(max_order, x).tolist())

    test_vals = np.swapaxes(np.array(test_vals), 0, 1)

    assert np.allclose(ref_vals, test_vals)


def test_evaluate_legendre():
    max_order = 10
    # Coefficients are set to 1, but will incorporate the (2l+1)/2 norm factor
    # for the reference solution
    test_coeffs = [0.5 * (2. * l + 1.) for l in range(max_order + 1)]
    test_xs = np.linspace(-1., 1., num=5, endpoint=True)

    ref_vals = np.polynomial.legendre.legval(test_xs, test_coeffs)

    # Set the coefficients back to 1s for the test values since
    # evaluate legendre incorporates the (2l+1)/2 term on its own
    test_coeffs = [1. for l in range(max_order + 1)]

    test_vals = np.array([openmc.lib.math.evaluate_legendre(test_coeffs, x)
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
    test_vals = openmc.lib.math.calc_rn(max_order, test_uvw)

    assert np.allclose(ref_vals, test_vals)


def test_calc_zn():
    n = 10
    rho = 0.5
    phi = 0.5

    # Reference solution from running the C++ implementation
    ref_vals = np.array([
        1.00000000e+00, 2.39712769e-01, 4.38791281e-01,
        2.10367746e-01, -5.00000000e-01, 1.35075576e-01,
        1.24686873e-01, -2.99640962e-01, -5.48489101e-01,
        8.84215021e-03, 5.68310892e-02, -4.20735492e-01,
        -1.25000000e-01, -2.70151153e-01, -2.60091773e-02,
        1.87022545e-02, -3.42888902e-01, 1.49820481e-01,
        2.74244551e-01, -2.43159131e-02, -2.50357380e-02,
        2.20500013e-03, -1.98908812e-01, 4.07587508e-01,
        4.37500000e-01, 2.61708929e-01, 9.10321205e-02,
        -1.54686328e-02, -2.74049397e-03, -7.94845816e-02,
        4.75368705e-01, 7.11647284e-02, 1.30266162e-01,
        3.37106977e-02, 1.06401886e-01, -7.31606787e-03,
        -2.95625975e-03, -1.10250006e-02, 3.55194307e-01,
        -1.44627826e-01, -2.89062500e-01, -9.28644588e-02,
        -1.62557358e-01, 7.73431638e-02, -2.55329539e-03,
        -1.90923851e-03, 1.57578403e-02, 1.72995854e-01,
        -3.66267690e-01, -1.81657333e-01, -3.32521518e-01,
        -2.59738162e-02, -2.31580576e-01, 4.20673902e-02,
        -4.11710546e-04, -9.36449487e-04, 1.92156884e-02,
        2.82515641e-02, -3.90713738e-01, -1.69280296e-01,
        -8.98437500e-02, -1.08693628e-01, 1.78813094e-01,
        -1.98191857e-01, 1.65964201e-02, 2.77013853e-04])

    test_vals = openmc.lib.math.calc_zn(n, rho, phi)

    assert np.allclose(ref_vals, test_vals)


def test_calc_zn_rad():
    n = 10
    rho = 0.5

    # Reference solution from running the C++ implementation
    ref_vals = np.array([
        1.00000000e+00, -5.00000000e-01, -1.25000000e-01,
        4.37500000e-01, -2.89062500e-01,-8.98437500e-02])

    test_vals = openmc.lib.math.calc_zn_rad(n, rho)

    assert np.allclose(ref_vals, test_vals)


def test_rotate_angle():
    uvw0 = np.array([1., 0., 0.])
    phi = 0.
    mu = 0.

    # reference: mu of 0 pulls the vector the bottom, so:
    ref_uvw = np.array([0., 0., -1.])

    test_uvw = openmc.lib.math.rotate_angle(uvw0, mu, phi)

    assert np.array_equal(ref_uvw, test_uvw)

    # Repeat for mu = 1 (no change)
    mu = 1.
    ref_uvw = np.array([1., 0., 0.])

    test_uvw = openmc.lib.math.rotate_angle(uvw0, mu, phi)

    assert np.array_equal(ref_uvw, test_uvw)

    # Now to test phi is None
    mu = 0.9
    phi = None
    prn_seed = 1

    # When seed = 1, phi will be sampled as 1.9116495709698769
    # The resultant reference is from hand-calculations given the above
    ref_uvw = [0.9, 0.410813051297112, 0.1457142302040]
    test_uvw = openmc.lib.math.rotate_angle(uvw0, mu, phi, prn_seed)

    assert np.allclose(ref_uvw, test_uvw)


def test_maxwell_spectrum():
    prn_seed = 1
    T = 0.5
    ref_val = 0.6129982175261098
    test_val = openmc.lib.math.maxwell_spectrum(T, prn_seed)

    assert ref_val == test_val


def test_watt_spectrum():
    prn_seed = 1
    a = 0.5
    b = 0.75
    ref_val = 0.6247242713640233
    test_val = openmc.lib.math.watt_spectrum(a, b, prn_seed)

    assert ref_val == test_val


def test_normal_dist():
    prn_seed = 1
    a = 14.08
    b = 0.0
    ref_val = 14.08
    test_val = openmc.lib.math.normal_variate(a, b, prn_seed)

    assert ref_val == pytest.approx(test_val)

    prn_seed = 1
    a = 14.08
    b = 1.0
    ref_val = 16.436645416691427
    test_val = openmc.lib.math.normal_variate(a, b, prn_seed)

    assert ref_val == pytest.approx(test_val)


def test_broaden_wmp_polynomials():
    # Two branches of the code to worry about, beta > 6 and otherwise
    # beta = sqrtE * dopp
    # First lets do beta > 6
    test_E = 0.5
    test_dopp = 100.  # approximately U235 at room temperature
    n = 6

    ref_val = [2., 1.41421356, 1.0001, 0.70731891, 0.50030001, 0.353907]
    test_val = openmc.lib.math.broaden_wmp_polynomials(test_E, test_dopp, n)

    assert np.allclose(ref_val, test_val)

    # now beta < 6
    test_dopp = 5.
    ref_val = [1.99999885, 1.41421356, 1.04, 0.79195959, 0.6224, 0.50346003]
    test_val = openmc.lib.math.broaden_wmp_polynomials(test_E, test_dopp, n)

    assert np.allclose(ref_val, test_val)
