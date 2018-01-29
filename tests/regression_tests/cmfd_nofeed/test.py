from tests.testing_harness import CMFDTestHarness


def test_cmfd_nofeed(request):
    harness = CMFDTestHarness('statepoint.20.h5')
    harness.request = request
    harness.main()
