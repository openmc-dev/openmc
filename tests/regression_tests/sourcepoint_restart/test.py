from tests.testing_harness import TestHarness


def test_sourcepoint_restart():
    harness = TestHarness('statepoint.10.h5')
    harness.main()
