from tests.testing_harness import TestHarness


def test_white_plane():
    harness = TestHarness('statepoint.10.h5')
    harness.main()
