from tests.testing_harness import TestHarness


def test_lattice_mixed(request):
    harness = TestHarness('statepoint.10.h5')
    harness.request = request
    harness.main()
