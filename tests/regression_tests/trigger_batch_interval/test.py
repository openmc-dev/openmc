from tests.testing_harness import TestHarness


def test_trigger_batch_interval(request):
    harness = TestHarness('statepoint.15.h5')
    harness.request = request
    harness.main()
