from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def test_mg_survival_biasing(request):
    model = slab_mg()
    model.settings.survival_biasing = True
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.request = request
    harness.main()
