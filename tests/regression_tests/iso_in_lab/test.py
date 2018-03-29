from tests.testing_harness import PyAPITestHarness


def test_iso_in_lab():
    # Force iso-in-lab scattering.
    harness = PyAPITestHarness('statepoint.10.h5')
    harness._model.materials.make_isotropic_in_lab()
    harness.main()
