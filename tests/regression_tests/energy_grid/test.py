from tests.testing_harness import TestHarness


def test_energy_grid():
    harness = TestHarness('statepoint.10.h5')
    harness.main()
