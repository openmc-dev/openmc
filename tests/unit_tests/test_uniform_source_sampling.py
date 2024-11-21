import openmc
import pytest


@pytest.fixture
def sphere_model():
    mat = openmc.Material()
    mat.add_nuclide('Li6', 1.0)
    mat.set_density('g/cm3', 1.0)
    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=mat)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])

    model.settings.particles = 100
    model.settings.batches = 1
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(1.0e3),
        strength=100.0
    )
    model.settings.run_mode = "fixed source"
    model.settings.surf_source_write = {
        "max_particles": 100,
    }

    tally = openmc.Tally()
    tally.scores = ['flux']
    model.tallies = [tally]
    return model


def test_source_weight(run_in_tmpdir, sphere_model):
    # Run OpenMC without uniform source sampling and check that banked particles
    # have weight 1
    sphere_model.settings.uniform_source_sampling = False
    sphere_model.run()
    particles = openmc.ParticleList.from_hdf5('surface_source.h5')
    assert set(p.wgt for p in particles) == {1.0}

    # Run with uniform source sampling and check that banked particles have
    # weight == strength
    sphere_model.settings.uniform_source_sampling = True
    sphere_model.run()
    particles = openmc.ParticleList.from_hdf5('surface_source.h5')
    strength = sphere_model.settings.source[0].strength
    assert set(p.wgt for p in particles) == {strength}


def test_tally_mean(run_in_tmpdir, sphere_model):
    # Run without uniform source sampling
    sphere_model.settings.uniform_source_sampling = False
    sp_file = sphere_model.run()
    with openmc.StatePoint(sp_file) as sp:
        reference_mean = sp.tallies[sphere_model.tallies[0].id].mean

    # Run with uniform source sampling
    sphere_model.settings.uniform_source_sampling = True
    sp_file = sphere_model.run()
    with openmc.StatePoint(sp_file) as sp:
        mean = sp.tallies[sphere_model.tallies[0].id].mean

    # Check that tally means match
    assert mean == pytest.approx(reference_mean)
