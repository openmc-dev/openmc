from tests.testing_harness import TestHarness


def test_lattice_corner_crossing():
    # Ensure we account for potential corner crossings
    # in floating point precision.
    harness = TestHarness('statepoint.10.h5')
    harness.main()
