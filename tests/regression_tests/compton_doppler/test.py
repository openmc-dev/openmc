import openmc

from tests.testing_harness import PyAPITestHarness


def test_compton_doppler():
    hydrogen = openmc.Material()
    hydrogen.add_nuclide('H1', 1.0)
    hydrogen.set_density('g/cm3', 1.0)

    sphere = openmc.Sphere(r=1.0e9, boundary_type='reflective')
    model = openmc.Model()
    model.geometry = openmc.Geometry([
        openmc.Cell(fill=hydrogen, region=-sphere)
    ])

    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.electron_treatment = 'led'
    model.settings.atomic_relaxation = False
    model.settings.batches = 1
    model.settings.particles = 10000
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        energy=openmc.stats.delta_function(2000.0),
    )

    tally = openmc.Tally()
    tally.filters = [
        openmc.ParticleFilter(['photon', 'electron']),
        openmc.EnergyFilter([0.0, 500.0, 1000.0, 1500.0, 2000.0]),
    ]
    tally.scores = ['flux', 'heating']
    model.tallies = [tally]

    harness = PyAPITestHarness('statepoint.1.h5', model=model)
    harness.main()
