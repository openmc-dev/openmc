"""Regression tests for openmc.deplete.integrator.si_celi algorithm.

These tests integrate a simple test problem described in dummy_geometry.py.
"""

from pytest import approx
from openmc.deplete import SICELIIntegrator, ResultsList

from tests import dummy_operator


def test_si_celi(run_in_tmpdir):
    """Integral regression test of integrator algorithm using si_celi"""

    op = dummy_operator.DummyOperator()
    op.output_dir = "test_integrator_regression"

    # Perform simulation using the si_celi algorithm
    dt = [0.75, 0.75]
    power = 1.0
    SICELIIntegrator(op, dt, power).integrate()

    # Load the files
    res = ResultsList.from_hdf5(op.output_dir / "depletion_results.h5")

    _, y1 = res.get_atoms("1", "1")
    _, y2 = res.get_atoms("1", "2")

    # Reference solution
    s1 = [2.03325094, 1.16826254]
    s2 = [2.69291933, 0.37907772]

    assert y1[1] == approx(s1[0])
    assert y2[1] == approx(s1[1])

    assert y1[2] == approx(s2[0])
    assert y2[2] == approx(s2[1])

    # Test structure of depletion time dataset

    dep_time = res.get_depletion_time()
    assert dep_time.shape == (len(dt), )
    assert all(dep_time > 0)
