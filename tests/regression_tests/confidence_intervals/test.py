from tests.testing_harness import TestHarness


def test_confidence_intervals(request):
    harness = TestHarness('statepoint.10.h5')
    harness.request = request
    harness.main()
