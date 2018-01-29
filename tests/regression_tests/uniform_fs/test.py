from tests.testing_harness import TestHarness


def test_uniform_fs(request):
    harness = TestHarness('statepoint.10.h5')
    harness.request = request
    harness.main()
