"""Regression tests for openmc.deplete restart capability.

These tests run in two steps, a first run then a restart run, a simple test
problem described in dummy_geometry.py.
"""

from pytest import approx
import openmc.deplete

from tests import dummy_operator


def test_restart_predictor(run_in_tmpdir):
    """Integral regression test of integrator algorithm using predictor."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_predictor"
    op.output_dir = output_dir

    # Perform simulation using the predictor algorithm
    dt = [0.75]
    power = 1.0
    openmc.deplete.predictor(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation using the predictor algorithm
    openmc.deplete.predictor(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Mathematica solution
    s1 = [2.46847546272295, 0.986431226850467]
    s2 = [4.11525874568034, -0.0581692232513460]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[2] == approx(s2[0])
    assert y2[2] == approx(s2[1])


def test_restart_cecm(run_in_tmpdir):
    """Integral regression test of integrator algorithm using CE/CM."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_cecm"
    op.output_dir = output_dir

    # Perform simulation using the MCNPX/MCNP6 algorithm
    dt = [0.75]
    power = 1.0
    openmc.deplete.cecm(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation using the MCNPX/MCNP6 algorithm
    openmc.deplete.cecm(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Mathematica solution
    s1 = [1.86872629872102, 1.395525772416039]
    s2 = [2.18097439443550, 2.69429754646747]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[3] == approx(s2[0])
    assert y2[3] == approx(s2[1])


def test_restart_predictor_cecm(run_in_tmpdir):
    """Integral regression test of integrator algorithm using predictor
     for the first run then CE/CM for the restart run."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_predictor_cecm"
    op.output_dir = output_dir

    # Perform simulation using the predictor algorithm
    dt = [0.75]
    power = 1.0
    openmc.deplete.predictor(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation using the MCNPX/MCNP6 algorithm
    openmc.deplete.cecm(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Test solution
    s1 = [2.46847546272295, 0.986431226850467]
    s2 = [3.09106948392, 0.607102912398]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[2] == approx(s2[0])
    assert y2[2] == approx(s2[1])


def test_restart_cecm_predictor(run_in_tmpdir):
    """Integral regression test of integrator algorithm using CE/CM for the
    first run then predictor for the restart run."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_cecm_predictor"
    op.output_dir = output_dir

    # Perform simulation using the MCNPX/MCNP6 algorithm
    dt = [0.75]
    power = 1.0
    openmc.deplete.cecm(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation using the predictor algorithm
    openmc.deplete.predictor(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Test solution
    s1 = [1.86872629872102, 1.395525772416039]
    s2 = [3.32776806576, 2.391425905]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[2] == approx(s2[0])
    assert y2[2] == approx(s2[1])
