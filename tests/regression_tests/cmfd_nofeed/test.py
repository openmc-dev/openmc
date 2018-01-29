from tests.testing_harness import CMFDTestHarness


def test_cmfd_nofeed():
    harness = CMFDTestHarness('statepoint.20.h5')
    harness.main()
