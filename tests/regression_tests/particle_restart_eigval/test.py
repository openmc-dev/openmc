from tests.testing_harness import ParticleRestartTestHarness


def test_particle_restart_eigval(request):
    harness = ParticleRestartTestHarness('particle_10_1030.h5')
    harness.request = request
    harness.main()
