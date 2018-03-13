from tests.testing_harness import TestHarness


def test_seed():
    harness = TestHarness('statepoint.10.h5')
    harness.main()
