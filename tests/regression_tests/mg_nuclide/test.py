from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def test_mg_nuclide():
    model = slab_mg(as_macro=False)
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
