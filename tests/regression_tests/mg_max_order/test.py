from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def test_mg_max_order(request):
    model = slab_mg(reps=['iso'])
    model.settings.max_order = 1
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.request = request
    harness.main()
