from tests.testing_harness import TestHarness


def test_dnp_drift_regular_mesh_nodal():
    harness = TestHarness('statepoint.20.h5')
    harness.main()
