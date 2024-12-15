from tests.testing_harness import ParticleRestartTestHarness


def test_particle_restart_fixed():
    harness = ParticleRestartTestHarness("particle_4_241.h5")
    harness.main()
