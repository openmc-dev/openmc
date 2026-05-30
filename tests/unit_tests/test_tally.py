import openmc


def test_tally_init_args():
    """Test that Tally constructor kwargs are applied correctly."""
    filter = openmc.EnergyFilter([0.0, 1.0, 20.0e6])
    tally = openmc.Tally(
        name='my tally',
        scores=['flux', 'fission'],
        filters=[filter],
        nuclides=['U235'],
        estimator='tracklength',
    )

    assert tally.name == 'my tally'
    assert tally.scores == ['flux', 'fission']
    assert tally.filters == [filter]
    assert tally.nuclides == ['U235']
    assert tally.estimator == 'tracklength'


def test_heating_estimator_with_photon_transport(run_in_tmpdir):
    """Test that a neutron-only heating tally uses the tracklength estimator
    even when photon transport is enabled.

    Without a neutron particle filter, photon heating requires the collision
    estimator (analog energy balance). But neutron heating has kerma
    coefficients that support tracklength scoring, so a neutron-only tally
    should keep the tracklength estimator.
    """
    mat = openmc.Material()
    mat.add_nuclide('Si28', 1.0)
    mat.set_density('g/cm3', 2.3)
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)

    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 100
    model.settings.batches = 1
    model.settings.photon_transport = True
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point()
    )

    neutron_filter = openmc.ParticleFilter(['neutron'])

    # Neutron-only heating — should get tracklength estimator
    t_neutron = openmc.Tally(name='heating_neutron')
    t_neutron.filters = [neutron_filter]
    t_neutron.scores = ['heating']

    # No particle filter — should get collision estimator
    t_all = openmc.Tally(name='heating_all')
    t_all.scores = ['heating']

    model.tallies = [t_neutron, t_all]

    sp_filename = model.run()
    with openmc.StatePoint(sp_filename) as sp:
        assert sp.tallies[t_neutron.id].estimator == 'tracklength'
        assert sp.tallies[t_all.id].estimator == 'collision'
