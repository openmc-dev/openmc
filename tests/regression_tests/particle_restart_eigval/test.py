from tests.testing_harness import ParticleRestartTestHarness


def test_particle_restart_eigval():
    harness = ParticleRestartTestHarness('particle_11_902.h5')
    harness.main()
