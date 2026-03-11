import numpy as np

import openmc


def test_compton_relaxation_heating():
    """Test that Compton scattering with atomic relaxation does not produce
    negative photon heating scores. When Doppler broadening causes the photon
    to transfer less energy than the shell binding energy, the electron is not
    ejected and atomic relaxation should not occur. Previously, relaxation ran
    unconditionally, creating secondaries whose energy exceeded the photon's
    energy transfer, producing negative heating.
    """
    # Lead has high-Z (82) with many inner shells and large binding energies,
    # maximizing the number of Compton events where the energy transfer is
    # less than the shell binding energy.
    lead = openmc.Material()
    lead.add_element('Pb', 1.0)
    lead.set_density('g/cm3', 11.35)

    box = openmc.model.RectangularParallelepiped(
        -10, 10, -10, 10, -10, 10, boundary_type='reflective')
    cell = openmc.Cell(region=-box, fill=lead)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])

    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        energy=openmc.stats.delta_function(200e3),
        particle='photon')
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 500

    mesh = openmc.RegularMesh.from_domain(model.geometry, dimension=(3, 3, 3))
    tally = openmc.Tally()
    tally.scores = ['heating']
    tally.filters = [
        openmc.ParticleFilter(['photon']),
        openmc.MeshFilter(mesh)
    ]
    model.tallies = [tally]

    sp_file = model.run()
    with openmc.StatePoint(sp_file) as sp:
        tally_mean = sp.tallies[tally.id].mean

    assert np.all(tally_mean >= 0), \
        "Negative photon heating from Compton + atomic relaxation"


def test_negative_positron_heating():
    m = openmc.Material()
    m.add_element('Li', 1.0)
    m.set_density('g/cm3', 10.0)

    surf = openmc.Sphere(r=100.0, boundary_type='reflective')
    cell = openmc.Cell(fill=m, region=-surf)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([5.0e6], [1.0]),
        particle='photon',
    )
    model.settings.particles = 7
    model.settings.batches = 1
    model.settings.electron_treatment = 'led'
    model.settings.seed = 513836

    tally = openmc.Tally()
    tally.filters = [openmc.ParticleFilter(['photon', 'electron', 'positron'])]
    tally.scores = ['heating']
    model.tallies = openmc.Tallies([tally])
    model.run(apply_tally_results=True)

    assert (tally.mean >= 0.0).all(), "Negative heating detected"
