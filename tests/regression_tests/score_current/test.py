from tests.testing_harness import HashedTestHarness


def test_score_current():
    harness = HashedTestHarness('statepoint.10.h5')
    harness.main()
