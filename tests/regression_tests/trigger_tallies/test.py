from tests.testing_harness import TestHarness


def test_trigger_tallies():
    harness = TestHarness('statepoint.15.h5')
    harness.main()
