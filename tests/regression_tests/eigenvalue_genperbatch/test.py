from tests.testing_harness import TestHarness


def test_eigenvalue_genperbatch():
    harness = TestHarness('statepoint.7.h5')
    harness.main()
