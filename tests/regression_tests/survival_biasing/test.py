from tests.testing_harness import TestHarness


def test_survival_biasing():
    harness = TestHarness('statepoint.10.h5')
    harness.main()
