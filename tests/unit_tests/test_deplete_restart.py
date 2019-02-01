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


def test_restart_cf4(run_in_tmpdir):
    """Integral regression test of integrator algorithm using CF4."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_cf4"
    op.output_dir = output_dir

    # Perform simulation
    dt = [0.75]
    power = 1.0
    openmc.deplete.cf4(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation
    openmc.deplete.cf4(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Reference solution
    s1 = [2.06101629, 1.37783588]
    s2 = [2.57241318, 2.63731630]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[3] == approx(s2[0])
    assert y2[3] == approx(s2[1])


def test_restart_epc_rk4(run_in_tmpdir):
    """Integral regression test of integrator algorithm using EPC-RK4."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_epc_rk4"
    op.output_dir = output_dir

    # Perform simulation
    dt = [0.75]
    power = 1.0
    openmc.deplete.epc_rk4(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation
    openmc.deplete.epc_rk4(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Reference solution
    s1 = [2.01978516, 1.42038037]
    s2 = [2.05246421, 3.06177191]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[3] == approx(s2[0])
    assert y2[3] == approx(s2[1])


def test_restart_celi(run_in_tmpdir):
    """Integral regression test of integrator algorithm using CELI."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_celi"
    op.output_dir = output_dir

    # Perform simulation
    dt = [0.75]
    power = 1.0
    openmc.deplete.celi(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation
    openmc.deplete.celi(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Reference solution
    s1 = [1.82078767, 0.97122898]
    s2 = [2.68441779, 0.05125966]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[3] == approx(s2[0])
    assert y2[3] == approx(s2[1])


def test_restart_leqi(run_in_tmpdir):
    """Integral regression test of integrator algorithm using LEQI."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_leqi"
    op.output_dir = output_dir

    # Perform simulation
    dt = [0.75]
    power = 1.0
    openmc.deplete.leqi(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation
    openmc.deplete.leqi(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Reference solution
    s1 = [1.82078767, 0.97122898]
    s2 = [2.74526197, 0.23339915]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[3] == approx(s2[0])
    assert y2[3] == approx(s2[1])

def test_restart_si_celi(run_in_tmpdir):
    """Integral regression test of integrator algorithm using SI-CELI."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_si_celi"
    op.output_dir = output_dir

    # Perform simulation
    dt = [0.75]
    power = 1.0
    openmc.deplete.si_celi(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation
    openmc.deplete.si_celi(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Reference solution
    s1 = [2.03325094, 1.16826254]
    s2 = [2.69291933, 0.37907772]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[3] == approx(s2[0])
    assert y2[3] == approx(s2[1])


def test_restart_si_leqi(run_in_tmpdir):
    """Integral regression test of integrator algorithm using SI-LEQI."""

    op = dummy_operator.DummyOperator()
    output_dir = "test_restart_si_leqi"
    op.output_dir = output_dir

    # Perform simulation
    dt = [0.75]
    power = 1.0
    openmc.deplete.si_leqi(op, dt, power, print_out=False)

    # Load the files
    prev_res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    # Re-create depletion operator and load previous results
    op = dummy_operator.DummyOperator(prev_res)
    op.output_dir = output_dir

    # Perform restarts simulation
    openmc.deplete.si_leqi(op, dt, power, print_out=False)

    # Load the files
    res = openmc.deplete.ResultsList(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Reference solution
    s1 = [2.03325094, 1.16826254]
    s2 = [2.92711288, 0.53753236]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[3] == approx(s2[0])
    assert y2[3] == approx(s2[1])
