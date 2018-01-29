from tests.testing_harness import TestHarness


def test_survival_biasing(request):
    harness = TestHarness('statepoint.10.h5')
    harness.request = request
    harness.main()
