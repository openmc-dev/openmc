from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def test_mg_legendre():
    model = slab_mg(reps=['iso'])
    model.settings.tabular_legendre = {'enable': False}

    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
