"""Tests the ResultsList class"""

from pathlib import Path
from math import inf

import numpy as np
import pytest
import openmc.deplete


@pytest.fixture
def res():
    """Load the reference results"""
    filename = (Path(__file__).parents[1] / 'regression_tests' / 'deplete'
                / 'test_reference.h5')
    return openmc.deplete.ResultsList.from_hdf5(filename)


def test_get_atoms(res):
    """Tests evaluating single nuclide concentration."""
    t, n = res.get_atoms("1", "Xe135")

    t_ref = np.array([0.0, 1296000.0, 2592000.0, 3888000.0])
    n_ref = np.array(
        [6.67473282e+08, 3.76987065e+14, 3.68586723e+14, 3.91338392e+14])

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(n, n_ref)

    # Check alternate units
    volume = res[0].volume["1"]

    t_days, n_cm3 = res.get_atoms("1", "Xe135", nuc_units="atom/cm3", time_units="d")

    assert t_days == pytest.approx(t_ref / (60 * 60 * 24))
    assert n_cm3 == pytest.approx(n_ref / volume)

    t_min, n_bcm = res.get_atoms("1", "Xe135", nuc_units="atom/b-cm", time_units="min")
    assert n_bcm == pytest.approx(n_cm3 * 1e-24)
    assert t_min == pytest.approx(t_ref / 60)

    t_hour, _n = res.get_atoms("1", "Xe135", time_units="h")
    assert t_hour == pytest.approx(t_ref / (60 * 60))


def test_get_reaction_rate(res):
    """Tests evaluating reaction rate."""
    t, r = res.get_reaction_rate("1", "Xe135", "(n,gamma)")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.67473282e+08, 3.76987065e+14, 3.68586723e+14, 3.91338392e+14]
    xs_ref = [3.32282064e-05, 2.76208092e-05, 4.10987995e-05, 3.72454755e-05]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(r, np.array(n_ref) * xs_ref)


def test_get_eigenvalue(res):
    """Tests evaluating eigenvalue."""
    t, k = res.get_eigenvalue()

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    k_ref = [1.16984322, 1.19097429, 1.03012517, 1.20045563]
    u_ref = [0.0375587, 0.03476389, 0.07215969, 0.02839639]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(k[:, 0], k_ref)
    np.testing.assert_allclose(k[:, 1], u_ref)


@pytest.mark.parametrize("unit", ("s", "d", "min", "h"))
def test_get_steps(unit):
    # Make a ResultsList full of near-empty Result instances
    # Just fill out a time schedule
    results = openmc.deplete.ResultsList()
    # Time in units of unit
    times = np.linspace(0, 100, num=5)
    if unit == "d":
        conversion_to_seconds = 60 * 60 * 24
    elif unit == "h":
        conversion_to_seconds = 60 * 60
    elif unit == "min":
        conversion_to_seconds = 60
    else:
        conversion_to_seconds = 1

    for ix in range(times.size):
        res = openmc.deplete.Results()
        res.time = times[ix:ix + 1] * conversion_to_seconds
        results.append(res)

    for expected, value in enumerate(times):
        actual = results.get_step_where(
            value, time_units=unit, atol=0, rtol=0)
        assert actual == expected, (value, results[actual].time[0])

    with pytest.raises(ValueError):
        # Emulate a result file with a non-zero initial point in time
        # as in starting from a restart
        results.get_step_where(times[0] - 1, time_units=unit, atol=0, rtol=0)

    with pytest.raises(ValueError):
        results.get_step_where(times[-1] + 1, time_units=unit, atol=0, rtol=0)

    # Grab intermediate points with a small offset
    delta = (times[1] - times[0])
    offset = delta * 0.1
    for expected, value in enumerate(times[1:-1], start=1):
        # Shoot a little low and a little high
        for mult in (1, -1):
            target = value + mult * offset
            # Compare using absolute and relative tolerances
            actual = results.get_step_where(
                target, time_units=unit, atol=offset * 2, rtol=inf)
            assert actual == expected, (
                target, times[actual], times[expected], offset)

            actual = results.get_step_where(
                target, time_units=unit, atol=inf, rtol=offset / value)
            assert actual == expected, (
                target, times[actual], times[expected], offset)
        # Check that the lower index is returned for the exact mid-point
        target = value + delta * 0.5
        actual = results.get_step_where(
            target, time_units=unit, atol=delta, rtol=delta / value)
        assert actual == expected

    # Shoot way over with no tolerance -> just give closest value
    actual = results.get_step_where(
        times[-1] * 100, time_units=unit, atol=inf, rtol=inf)
    assert actual == times.size - 1
