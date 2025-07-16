"""
Initially from Jingang Liang: https://github.com/mit-crpg/vectfit.git
"""
import vectfit as m
import numpy as np
import pytest

@pytest.fixture
def ref_poles():
    """Reference poles for real-pole test."""
    return np.array([
             9.709261771920490e+02 + 0.0j,
            -1.120960794075339e+03 + 0.0j,
             1.923889557426567e+00 + 7.543700246109742e+01j,
             1.923889557426567e+00 - 7.543700246109742e+01j,
             1.159741300380281e+02 + 3.595650922556496e-02j,
             1.159741300380281e+02 - 3.595650922556496e-02j,
             1.546932165729394e+02 + 8.728391144940301e-02j,
             1.546932165729394e+02 - 8.728391144940301e-02j,
             2.280349190818197e+02 + 2.814037559718684e-01j,
             2.280349190818197e+02 - 2.814037559718684e-01j,
             2.313004772627853e+02 + 3.004628477692201e-01j,
             2.313004772627853e+02 - 3.004628477692201e-01j,
             2.787470098364861e+02 + 3.414179169920170e-01j,
             2.787470098364861e+02 - 3.414179169920170e-01j,
             3.570711338764254e+02 + 4.485587371149193e-01j,
             3.570711338764254e+02 - 4.485587371149193e-01j,
             4.701059001346060e+02 + 6.598089307174224e-01j,
             4.701059001346060e+02 - 6.598089307174224e-01j,
             7.275819506342254e+02 + 1.189678974845038e+03j,
             7.275819506342254e+02 - 1.189678974845038e+03j
        ])

@pytest.fixture
def ref_residues():
    """Reference residues for real-pole test."""
    return np.array([[
            -3.269879776751686e+07 + 0.0j,
             1.131087935798761e+09 + 0.0j,
             1.634151281869857e+04 + 2.251103589277891e+05j,
             1.634151281869857e+04 - 2.251103589277891e+05j,
             3.281792303833561e+03 - 1.756079516325274e+04j,
             3.281792303833561e+03 + 1.756079516325274e+04j,
             1.110800880243503e+04 - 4.324813594540043e+04j,
             1.110800880243503e+04 + 4.324813594540043e+04j,
             8.812700704117636e+04 - 2.256520243571103e+05j,
             8.812700704117636e+04 + 2.256520243571103e+05j,
             5.842090495551535e+04 - 1.442159380741478e+05j,
             5.842090495551535e+04 + 1.442159380741478e+05j,
             1.339410514130921e+05 - 2.640767909713812e+05j,
             1.339410514130921e+05 + 2.640767909713812e+05j,
             2.211245633333130e+05 - 3.222447758311512e+05j,
             2.211245633333130e+05 + 3.222447758311512e+05j,
             4.124430059785149e+05 - 4.076023108323907e+05j,
             4.124430059785149e+05 + 4.076023108323907e+05j,
             1.607378314999252e+09 - 1.401163320110452e+08j,
             1.607378314999252e+09 + 1.401163320110452e+08j
        ]])

@pytest.fixture
def vector_test_data():
    """Simple 2-signal test with known poles and residues."""
    Ns = 101
    s = np.linspace(3., 7., Ns)
    poles = [5.0 + 0.1j, 5.0 - 0.1j]
    residues = [[0.5 - 11.0j, 0.5 + 11.0j],
                [1.5 - 20.0j, 1.5 + 20.0j]]
    f = np.zeros((2, Ns))
    for i in range(2):
        f[i, :] = np.real(residues[i][0]/(s - poles[0]) + residues[i][1]/(s - poles[1]))
    weight = 1.0 / f
    init_poles = [3.5 + 0.035j, 3.5 - 0.035j]
    return s, poles, residues, f, weight, init_poles

@pytest.fixture
def poly_test_data():
    """Test data with rational function plus polynomial terms."""
    Ns = 201
    s = np.linspace(0., 5., Ns)
    poles = [-20.0 + 30.0j, -20.0 - 30.0j]
    residues = [[5.0 + 10.0j, 5.0 - 10.0j]]
    polys = [[1.0, 2.0, 0.3]]
    f = m.evaluate(s, poles, residues, polys)
    weight = 1.0 / f
    init_poles = [2.5 + 0.025j, 2.5 - 0.025j]
    return s, poles, residues, polys, f, weight, init_poles

@pytest.fixture
def real_poles_data(ref_poles, ref_residues):
    """Large-scale signal using complex and real poles."""
    Ns = 5000
    s = np.linspace(1.0e-2, 5.e3, Ns)
    f = np.zeros((1, Ns))
    for p, r in zip(ref_poles, ref_residues[0]):
        f[0] += (r / (s - p)).real
    weight = 1.0 / f
    poles = np.linspace(1.1e-2, 4.8e3, 10)
    poles = poles + poles * 0.01j
    poles = np.sort(np.append(poles, np.conj(poles)))
    return s, f, weight, pole


@pytest.fixture
def large_test_data():
    """Stress test data with thousands of poles and samples."""
    Ns = 20000
    N = 1000
    s = np.linspace(1.0e-2, 5.e3, Ns)
    poles = np.linspace(1.1e-2, 4.8e3, N // 2) + 0.01j * np.linspace(1.1e-2, 4.8e3, N // 2)
    poles = np.sort(np.append(poles, np.conj(poles)))
    residues = np.linspace(1e2, 1e6, N // 2) + 0.5j * np.linspace(1e2, 1e6, N // 2)
    residues = np.sort(np.append(residues, np.conj(residues))).reshape((1, N))
    f = np.zeros((1, Ns))
    for p, r in zip(poles, residues[0]):
        f[0] += (r / (s - p)).real
    weight = 1.0 / f
    init_poles = np.linspace(1.2e-2, 4.7e3, N // 2) + 0.01j * np.linspace(1.2e-2, 4.7e3, N // 2)
    init_poles = np.sort(np.append(init_poles, np.conj(init_poles)))
    return s, f, weight, init_poles

@pytest.fixture
def eval_test_data():
    """Reference data for evaluating rational + polynomial models."""
    Ns = 101
    s = np.linspace(-5., 5., Ns)
    poles = [-2.0 + 30.0j, -2.0 - 30.0j]
    residues = [5.0 + 10.0j, 5.0 - 10.0j]
    polys = [1.0, 2.0, 0.3]
    return s, poles, residues, polys

def test_vector(vector_test_data):
    """Test vectfit with vector samples and simple poles.
       It is expected to get exact results with one iteration.
    """
    s, poles, residues, f, weight, init_poles = vector_test_data
    poles, residues, cf, fit, rms = m.vectfit(f, test_s, init_poles, weight)
    assert np.allclose(test_poles, poles, rtol=1e-7)
    assert np.allclose(test_residues, residues, rtol=1e-7)
    assert np.allclose(f, fit, rtol=1e-5)


def test_poly(poly_test_data):
    """Test vectfit with polynomials."""
    s, poles, residues, polys, f, weight, init_poles = poly_test_data
    poles, residues, cf, fit, rms = m.vectfit(f, test_s, init_poles, weight, n_polys=3)
    poles, residues, cf, fit, rms = m.vectfit(f, test_s, poles, weight, n_polys=3)
    assert np.allclose(test_poles, poles, rtol=1e-5)
    assert np.allclose(test_residues, residues, rtol=1e-5)
    assert np.allclose(test_polys, cf, rtol=1e-5)
    assert np.allclose(f, fit, rtol=1e-4)


def test_real_poles(real_poles_data, ref_poles, ref_residues):
    """Test vectfit with more poles including real poles"""
    s, f, weight, poles = real_poles_data
    for _ in range(10):
        poles, residues, _, fit, _ = m.vectfit(f, s, poles, weight)
    assert np.allclose(np.sort(ref_poles), np.sort(poles))
    assert np.allclose(np.sort(ref_residues), np.sort(residues))
    assert np.allclose(f, fit, rtol=1e-3)


def test_large(large_test_data):
    """Test vectfit with a large set of poles and samples"""
    s, f, weight, init_poles = large_test_data
    poles_fit, residues_fit, cf, f_fit, rms = m.vectfit(f, s, poles_init, weight)
    assert np.allclose(f, f_fit, rtol=1e-3)


def test_evaluate(eval_test_data):
    """Test evaluate function"""
    # Single signal, no polynomial
    s, poles, residues, polys = eval_test_data
    f_ref[0, :] = np.real(residues[0]/(s - poles[0]) +
                          residues[1]/(s - poles[1]))
    f = m.evaluate(s, poles, residues)
    assert np.allclose(f_ref, f)

    # Single signal, with polynomial
    for n, c in enumerate(polys):
        f_ref[0, :] += c*np.power(s, n)
    f = m.evaluate(s, poles, residues, polys)
    assert np.allclose(f_ref, f)

    # Multi-signal, multi-residue, multi-poly
    poles = [5.0+0.1j, 5.0-0.1j]
    residues = [[0.5-11.0j, 0.5+11.0j],
                [1.5-20.0j, 1.5+20.0j]]
    polys = [[1.0, 2.0, 0.3], [4.0, -2.0, -10.0]]
    f_ref = np.zeros([2, Ns])
    f_ref[0, :] = np.real(residues[0][0]/(s - poles[0]) +
                          residues[0][1]/(s - poles[1]))
    f_ref[1, :] = np.real(residues[1][0]/(s - poles[0]) +
                          residues[1][1]/(s - poles[1]))
    for n, c in enumerate(polys[0]):
        f_ref[0, :] += c*np.power(s, n)
    for n, c in enumerate(polys[1]):
        f_ref[1, :] += c*np.power(s, n)
    f = m.evaluate(s, poles, residues, polys)
    assert np.allclose(f_ref, f)
