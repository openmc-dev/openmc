from tests.testing_harness import TestHarness


def test_trigger_no_status():
    harness = TestHarness('statepoint.10.h5')
    harness.main()
