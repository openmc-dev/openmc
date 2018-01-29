from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def test_mg_basic():
    model = slab_mg()
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
