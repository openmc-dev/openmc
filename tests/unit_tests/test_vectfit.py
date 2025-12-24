"""
Initially from Jingang Liang: https://github.com/mit-crpg/vectfit.git
"""

import numpy as np
import pytest
from openmc.data.vectfit import evaluate, vectfit


@pytest.fixture
def ref_poles():
    """Reference poles for real-pole test."""
    return np.array(
        [
            9.709261771920490e02 + 0.0j,
            -1.120960794075339e03 + 0.0j,
            1.923889557426567e00 + 7.543700246109742e01j,
            1.923889557426567e00 - 7.543700246109742e01j,
            1.159741300380281e02 + 3.595650922556496e-02j,
            1.159741300380281e02 - 3.595650922556496e-02j,
            1.546932165729394e02 + 8.728391144940301e-02j,
            1.546932165729394e02 - 8.728391144940301e-02j,
            2.280349190818197e02 + 2.814037559718684e-01j,
            2.280349190818197e02 - 2.814037559718684e-01j,
            2.313004772627853e02 + 3.004628477692201e-01j,
            2.313004772627853e02 - 3.004628477692201e-01j,
            2.787470098364861e02 + 3.414179169920170e-01j,
            2.787470098364861e02 - 3.414179169920170e-01j,
            3.570711338764254e02 + 4.485587371149193e-01j,
            3.570711338764254e02 - 4.485587371149193e-01j,
            4.701059001346060e02 + 6.598089307174224e-01j,
            4.701059001346060e02 - 6.598089307174224e-01j,
            7.275819506342254e02 + 1.189678974845038e03j,
            7.275819506342254e02 - 1.189678974845038e03j,
        ]
    )


@pytest.fixture
def ref_residues():
    """Reference residues for real-pole test."""
    return np.array(
        [
            [
                -3.269879776751686e07 + 0.0j,
                1.131087935798761e09 + 0.0j,
                1.634151281869857e04 + 2.251103589277891e05j,
                1.634151281869857e04 - 2.251103589277891e05j,
                3.281792303833561e03 - 1.756079516325274e04j,
                3.281792303833561e03 + 1.756079516325274e04j,
                1.110800880243503e04 - 4.324813594540043e04j,
                1.110800880243503e04 + 4.324813594540043e04j,
                8.812700704117636e04 - 2.256520243571103e05j,
                8.812700704117636e04 + 2.256520243571103e05j,
                5.842090495551535e04 - 1.442159380741478e05j,
                5.842090495551535e04 + 1.442159380741478e05j,
                1.339410514130921e05 - 2.640767909713812e05j,
                1.339410514130921e05 + 2.640767909713812e05j,
                2.211245633333130e05 - 3.222447758311512e05j,
                2.211245633333130e05 + 3.222447758311512e05j,
                4.124430059785149e05 - 4.076023108323907e05j,
                4.124430059785149e05 + 4.076023108323907e05j,
                1.607378314999252e09 - 1.401163320110452e08j,
                1.607378314999252e09 + 1.401163320110452e08j,
            ]
        ]
    )


@pytest.fixture
def vector_test_data():
    """Simple 2-signal test with known poles and residues."""
    Ns = 101
    s = np.linspace(3.0, 7.0, Ns)
    poles = [5.0 + 0.1j, 5.0 - 0.1j]
    residues = [[0.5 - 11.0j, 0.5 + 11.0j], [1.5 - 20.0j, 1.5 + 20.0j]]
    f = np.zeros((2, Ns))
    for i in range(2):
        f[i, :] = np.real(
            residues[i][0] / (s - poles[0]) + residues[i][1] / (s - poles[1])
        )
    weight = 1.0 / f
    init_poles = [3.5 + 0.035j, 3.5 - 0.035j]
    return s, poles, residues, f, weight, init_poles


@pytest.fixture
def poly_test_data():
    """Test data with rational function plus polynomial terms."""
    Ns = 201
    s = np.linspace(0.0, 5.0, Ns)
    poles = [-20.0 + 30.0j, -20.0 - 30.0j]
    residues = [[5.0 + 10.0j, 5.0 - 10.0j]]
    polys = [[1.0, 2.0, 0.3]]
    f = evaluate(s, poles, residues, polys)
    weight = 1.0 / f
    init_poles = [2.5 + 0.025j, 2.5 - 0.025j]
    return s, poles, residues, polys, f, weight, init_poles


@pytest.fixture
def real_poles_data(ref_poles, ref_residues):
    """Large-scale signal using complex and real poles."""
    Ns = 5000
    s = np.linspace(1.0e-2, 5.0e3, Ns)
    f = np.zeros((1, Ns))
    for p, r in zip(ref_poles, ref_residues[0]):
        f[0] += (r / (s - p)).real
    weight = 1.0 / f
    poles = np.linspace(1.1e-2, 4.8e3, 10)
    poles = poles + poles * 0.01j
    poles = np.sort(np.append(poles, np.conj(poles)))
    return s, f, weight, poles


@pytest.fixture
def large_test_data():
    """Stress test data with thousands of poles and samples."""
    Ns = 20000
    N = 1000
    s = np.linspace(1.0e-2, 5.0e3, Ns)
    poles = np.linspace(1.1e-2, 4.8e3, N // 2) + 0.01j * np.linspace(
        1.1e-2, 4.8e3, N // 2
    )
    poles = np.sort(np.append(poles, np.conj(poles)))
    residues = np.linspace(1e2, 1e6, N // 2) + 0.5j * np.linspace(1e2, 1e6, N // 2)
    residues = np.sort(np.append(residues, np.conj(residues))).reshape((1, N))
    f = np.zeros((1, Ns))
    for p, r in zip(poles, residues[0]):
        f[0] += (r / (s - p)).real
    weight = 1.0 / f
    init_poles = np.linspace(1.2e-2, 4.7e3, N // 2) + 0.01j * np.linspace(
        1.2e-2, 4.7e3, N // 2
    )
    init_poles = np.sort(np.append(init_poles, np.conj(init_poles)))
    return s, f, weight, init_poles


@pytest.fixture
def eval_test_data():
    """Reference data for evaluating rational + polynomial models."""
    Ns = 101
    s = np.linspace(-5.0, 5.0, Ns)
    poles = [-2.0 + 30.0j, -2.0 - 30.0j]
    residues = [5.0 + 10.0j, 5.0 - 10.0j]
    polys = [1.0, 2.0, 0.3]
    return s, poles, residues, polys


def test_vector(vector_test_data):
    """Test vectfit with vector samples and simple poles.
    It is expected to get exact results with one iteration.
    """
    s, expected_poles, expected_residues, f, weight, init_poles = vector_test_data
    poles, residues, _, fit, _ = vectfit(f, s, init_poles, weight)
    assert np.allclose(poles, expected_poles, rtol=1e-7)
    assert np.allclose(residues, expected_residues, rtol=1e-7)
    assert np.allclose(f, fit, rtol=1e-5)


def test_poly(poly_test_data):
    """Test vectfit with polynomials."""
    s, expected_poles, expected_residues, expected_polys, f, weight, init_poles = (
        poly_test_data
    )
    poles, residues, cf, fit, _ = vectfit(f, s, init_poles, weight, n_polys=3)
    poles, residues, cf, fit, _ = vectfit(f, s, poles, weight, n_polys=3)
    assert np.allclose(poles, expected_poles, rtol=1e-5)
    assert np.allclose(residues, expected_residues, rtol=1e-5)
    assert np.allclose(cf, expected_polys, rtol=1e-5)
    assert np.allclose(f, fit, rtol=1e-4)


def test_real_poles(real_poles_data, ref_poles, ref_residues):
    """Test vectfit with more poles including real poles"""
    s, f, weight, poles = real_poles_data
    for _ in range(10):
        poles, residues, _, fit, _ = vectfit(f, s, poles, weight)
    assert np.allclose(np.sort(poles), np.sort(ref_poles))
    assert np.allclose(np.sort(residues), np.sort(ref_residues))
    assert np.allclose(f, fit, rtol=1e-3)


def test_large(large_test_data):
    """Test vectfit with a large set of poles and samples"""
    s, f, weight, init_poles = large_test_data
    poles_fit, residues_fit, _, f_fit, _ = vectfit(f, s, init_poles, weight)
    assert np.allclose(f, f_fit, rtol=1e-3)


def test_evaluate(eval_test_data):
    """Test evaluate function"""
    s, poles, residues, polys = eval_test_data

    # Single signal, no polynomial
    f_ref = np.real(residues[0] / (s - poles[0]) + residues[1] / (s - poles[1]))
    f = evaluate(s, poles, residues)
    assert np.allclose(f[0], f_ref)

    # Single signal, with polynomial
    for n, c in enumerate(polys):
        f_ref += c * np.power(s, n)
    f = evaluate(s, poles, residues, polys)
    assert np.allclose(f[0], f_ref)

    # Multi-signal, multi-residue, multi-poly
    poles = [5.0 + 0.1j, 5.0 - 0.1j]
    residues = [[0.5 - 11.0j, 0.5 + 11.0j], [1.5 - 20.0j, 1.5 + 20.0j]]
    polys = [[1.0, 2.0, 0.3], [4.0, -2.0, -10.0]]
    f_ref = np.zeros((2, len(s)))
    for i in range(2):
        f_ref[i, :] = np.real(
            residues[i][0] / (s - poles[0]) + residues[i][1] / (s - poles[1])
        )
        for n, c in enumerate(polys[i]):
            f_ref[i, :] += c * np.power(s, n)
    f = evaluate(s, poles, residues, polys)
    assert np.allclose(f, f_ref)
