from tests.testing_harness import ParticleRestartTestHarness


def test_particle_restart_fixed_shared_secondary():
    harness = ParticleRestartTestHarness('particle_4_3241.h5')
    harness.main()
