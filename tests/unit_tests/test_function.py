#!/usr/bin/env python

import numpy as np
import pytest

import openmc.data


def test_tabulated1d_scalar_array_consistency():
    """Scalar and array evaluation should agree everywhere, including
    outside the tabulated domain (issue #4041)."""
    f = openmc.data.Tabulated1D([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])

    # Out-of-range scalars return nearest endpoint value
    assert f(-100.0) == 4.0
    assert f(0.999) == 4.0
    assert f(3.001) == 6.0
    assert f(1e300) == 6.0

    xs = np.array([-100.0, 0.0, 0.999, 1.0, 1.5, 2.0, 2.5,
                   3.0, 3.000000001, 100.0])
    assert np.all(f(xs) == np.array([f(x) for x in xs]))


@pytest.mark.parametrize('interp', [1, 2, 3, 4, 5])
def test_tabulated1d_interpolation_schemes(interp):
    """All ENDF interpolation schemes give identical scalar/array results."""
    f = openmc.data.Tabulated1D([1.0, 2.0, 3.0], [1.0, 4.0, 9.0],
                                breakpoints=[3], interpolation=[interp])

    xs = np.array([-10.0, 0.5, 1.0, 1.1, 1.5, 1.9, 2.0, 2.5, 2.99, 3.0, 50.0])
    expected = np.array([f(x) for x in xs])
    assert np.allclose(f(xs), expected, rtol=1e-14)

    # Values within the tabulated range should not be zero-filled
    assert np.all(f(xs[2:-1]) != 0.0)

    if interp == 1:
        # Histogram
        assert f(1.5) == 1.0
    elif interp == 2:
        # Linear-linear
        assert f(1.5) == pytest.approx(2.5)
    elif interp == 4:
        # Log-linear
        assert f(1.5) == pytest.approx(2.0)


def test_tabulated1d_multiple_regions():
    """Interpolation regions are respected by both evaluation paths."""
    f = openmc.data.Tabulated1D(
        [1.0, 2.0, 3.0, 4.0, 5.0], [10.0, 20.0, 30.0, 40.0, 80.0],
        breakpoints=[3, 5], interpolation=[2, 1])

    assert f(1.5) == pytest.approx(15.0)
    assert f(3.5) == 30.0  # histogram region

    xs = np.array([0.0, 1.5, 2.5, 3.0, 3.5, 4.9, 5.0, 6.0])
    assert np.all(f(xs) == np.array([f(x) for x in xs]))


def test_tabulated1d_endpoints():
    """Exact endpoints and precision-level perturbations return endpoint
    values."""
    f = openmc.data.Tabulated1D([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])

    assert np.array_equal(f(np.array([1.0, 3.0])), np.array([4.0, 6.0]))
    assert f(1.0) == 4.0
    assert f(3.0) == 6.0

    # Slightly beyond the upper endpoint due to floating point precision
    assert f(np.array([3.0*(1.0 + 1e-15)]))[0] == 6.0
    assert f(np.array([3.0 + 1e-6]))[0] == 6.0


def test_tabulated1d_multidimensional_input():
    """Array shape is preserved through evaluation."""
    f = openmc.data.Tabulated1D([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])

    x = np.array([[0.0, 1.5], [2.5, 9.0]])
    y = f(x)
    assert y.shape == (2, 2)
    assert np.array_equal(y, np.array([[4.0, 4.5], [5.5, 6.0]]))


def test_sum_functions_partial_domains():
    """Combining tabulated functions with differing domains leaves zeros
    where a component is undefined."""
    f1 = openmc.data.Tabulated1D([1.0, 2.0, 3.0], [10.0, 20.0, 30.0])
    f2 = openmc.data.Tabulated1D([2.0, 3.0, 4.0], [100.0, 200.0, 400.0])

    s = openmc.data.sum_functions([f1, f2])
    assert isinstance(s, openmc.data.Tabulated1D)
    assert np.array_equal(s.x, np.array([1.0, 2.0, 3.0, 4.0]))

    # Where both functions are defined, they add; outside a function's own
    # domain its contribution is zero
    assert s.y[0] == pytest.approx(10.0)   # f2 undefined at 1.0
    assert s.y[1] == pytest.approx(120.0)
    assert s.y[2] == pytest.approx(230.0)
    assert s.y[3] == pytest.approx(400.0)  # f1 undefined at 4.0


def test_sum_functions_polynomial():
    """Polynomial contributions apply over the full union grid."""
    p = openmc.data.Polynomial((1.0, -0.5))
    f = openmc.data.Tabulated1D([2.0, 4.0], [10.0, 20.0])

    s = openmc.data.sum_functions([f, p])
    assert np.array_equal(s.x, np.array([2.0, 4.0]))
    assert s.y[0] == pytest.approx(10.0)
    assert s.y[1] == pytest.approx(19.0)
