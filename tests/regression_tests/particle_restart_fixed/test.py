from tests.testing_harness import ParticleRestartTestHarness


def test_particle_restart_fixed():
    harness = ParticleRestartTestHarness('particle_7_144.h5')
    harness.main()
