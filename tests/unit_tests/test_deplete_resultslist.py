"""Tests the Results class"""

from pathlib import Path
from math import inf

import numpy as np
import pytest
import openmc.deplete


@pytest.fixture
def res():
    """Load the reference results"""
    filename = (
        Path(__file__).parents[1]
        / "regression_tests"
        / "deplete_with_transport"
        / "test_reference.h5"
    )
    return openmc.deplete.Results(filename)


def test_get_activity(res):
    """Tests evaluating activity"""
    t, a = res.get_activity("1")

    t_ref = np.array([0.0, 1296000.0, 2592000.0, 3888000.0])
    a_ref = np.array([1.25167956e06, 3.71938527e11, 4.43264300e11, 3.55547176e11])

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(a, a_ref)

    # Check by_nuclide
    a_xe135_ref = np.array(
        [2.106574218e05, 1.227519888e11, 1.177491828e11, 1.031986176e11]
    )
    t_nuc, a_nuc = res.get_activity("1", by_nuclide=True)

    a_xe135 = np.array([a_nuc_i["Xe135"] for a_nuc_i in a_nuc])

    np.testing.assert_allclose(t_nuc, t_ref)
    np.testing.assert_allclose(a_xe135, a_xe135_ref)


def test_get_atoms(res):
    """Tests evaluating single nuclide concentration."""
    t, n = res.get_atoms("1", "Xe135")

    t_ref = np.array([0.0, 1296000.0, 2592000.0, 3888000.0])
    n_ref = np.array([6.67473282e08, 3.88942731e14, 3.73091215e14, 3.26987387e14])

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


def test_get_decay_heat(res):
    """Tests evaluating decay heat."""
    # Set chain file for testing
    openmc.config["chain_file"] = Path(__file__).parents[1] / "chain_simple.xml"

    t_ref = np.array([0.0, 1296000.0, 2592000.0, 3888000.0])
    dh_ref = np.array([1.27933813e-09, 5.85347232e-03, 7.38773010e-03, 5.79954067e-03])

    t, dh = res.get_decay_heat("1")

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(dh, dh_ref)

    # Check by nuclide
    dh_xe135_ref = np.array(
        [1.27933813e-09, 7.45481920e-04, 7.15099509e-04, 6.26732849e-04]
    )
    t_nuc, dh_nuc = res.get_decay_heat("1", by_nuclide=True)

    dh_nuc_xe135 = np.array([dh_nuc_i["Xe135"] for dh_nuc_i in dh_nuc])

    np.testing.assert_allclose(t_nuc, t_ref)
    np.testing.assert_allclose(dh_nuc_xe135, dh_xe135_ref)


def test_get_mass(res):
    """Tests evaluating single nuclide concentration."""
    t, n = res.get_mass("1", "Xe135")

    t_ref = np.array([0.0, 1296000.0, 2592000.0, 3888000.0])
    n_ref = np.array([6.67473282e08, 3.88942731e14, 3.73091215e14, 3.26987387e14])

    # Get g
    n_ref *= openmc.data.atomic_mass("Xe135") / openmc.data.AVOGADRO

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(n, n_ref)

    # Check alternate units
    volume = res[0].volume["1"]
    t_days, n_cm3 = res.get_mass("1", "Xe135", mass_units="g/cm3", time_units="d")

    assert t_days == pytest.approx(t_ref / (60 * 60 * 24))
    assert n_cm3 == pytest.approx(n_ref / volume)

    t_min, n_bcm = res.get_mass("1", "Xe135", mass_units="kg", time_units="min")
    assert n_bcm == pytest.approx(n_ref / 1e3)
    assert t_min == pytest.approx(t_ref / 60)

    t_hour, _n = res.get_mass("1", "Xe135", time_units="h")
    assert t_hour == pytest.approx(t_ref / (60 * 60))


def test_get_reaction_rate(res):
    """Tests evaluating reaction rate."""
    t, r = res.get_reaction_rate("1", "Xe135", "(n,gamma)")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.67473282e08, 3.88942731e14, 3.73091215e14, 3.26987387e14]
    xs_ref = [2.53336104e-05, 4.21747011e-05, 3.48616127e-05, 3.61775563e-05]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(r, np.array(n_ref) * xs_ref)


def test_get_keff(res):
    """Tests evaluating keff."""
    t, k = res.get_keff()
    t_min, k = res.get_keff(time_units="min")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    k_ref = [1.1596402556, 1.1914183335, 1.2292570871, 1.1797030302]
    u_ref = [0.0270680649, 0.0219163444, 0.024268508, 0.0221401194]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(t_min * 60, t_ref)
    np.testing.assert_allclose(k[:, 0], k_ref)
    np.testing.assert_allclose(k[:, 1], u_ref)


@pytest.mark.parametrize("unit", ("s", "d", "min", "h", "a"))
def test_get_steps(unit):
    # Make a Results full of near-empty Result instances
    # Just fill out a time schedule
    results = openmc.deplete.Results(filename=None)
    # Time in units of unit
    times = np.linspace(0, 100, num=5)
    if unit == "a":
        conversion_to_seconds = 60 * 60 * 24 * 365.25
    elif unit == "d":
        conversion_to_seconds = 60 * 60 * 24
    elif unit == "h":
        conversion_to_seconds = 60 * 60
    elif unit == "min":
        conversion_to_seconds = 60
    else:
        conversion_to_seconds = 1

    for ix in range(times.size):
        res = openmc.deplete.StepResult()
        res.time = times[ix : ix + 1] * conversion_to_seconds
        results.append(res)

    for expected, value in enumerate(times):
        actual = results.get_step_where(value, time_units=unit, atol=0, rtol=0)
        assert actual == expected, (value, results[actual].time[0])

    with pytest.raises(ValueError):
        # Emulate a result file with a non-zero initial point in time
        # as in starting from a restart
        results.get_step_where(times[0] - 1, time_units=unit, atol=0, rtol=0)

    with pytest.raises(ValueError):
        results.get_step_where(times[-1] + 1, time_units=unit, atol=0, rtol=0)

    # Grab intermediate points with a small offset
    delta = times[1] - times[0]
    offset = delta * 0.1
    for expected, value in enumerate(times[1:-1], start=1):
        # Shoot a little low and a little high
        for mult in (1, -1):
            target = value + mult * offset
            # Compare using absolute and relative tolerances
            actual = results.get_step_where(
                target, time_units=unit, atol=offset * 2, rtol=inf
            )
            assert actual == expected, (target, times[actual], times[expected], offset)

            actual = results.get_step_where(
                target, time_units=unit, atol=inf, rtol=offset / value
            )
            assert actual == expected, (target, times[actual], times[expected], offset)
        # Check that the lower index is returned for the exact mid-point
        target = value + delta * 0.5
        actual = results.get_step_where(
            target, time_units=unit, atol=delta, rtol=delta / value
        )
        assert actual == expected

    # Shoot way over with no tolerance -> just give closest value
    actual = results.get_step_where(
        times[-1] * 100, time_units=unit, atol=inf, rtol=inf
    )
    assert actual == times.size - 1


def test_stepresult_get_material(res):
    # Get material at first timestep
    step_result = res[0]
    mat1 = step_result.get_material("1")
    assert mat1.id == 1
    assert mat1.volume == step_result.volume["1"]

    # Spot check number densities
    densities = mat1.get_nuclide_atom_densities()
    assert densities["Xe135"] == pytest.approx(1e-14)
    assert densities["I135"] == pytest.approx(1e-21)
    assert densities["U234"] == pytest.approx(1.00506e-05)
