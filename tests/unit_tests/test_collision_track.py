"""Test the 'collision_track' setting used to store particle information
during specified collision conditions in a file for a given simulation."""

import openmc
import pytest
import h5py
import numpy as np
import shutil

from tests.testing_harness import CollisionTrackTestHarness as ctt


@pytest.fixture(scope="module")
def geometry():
    """Simple hydrogen sphere geometry"""
    openmc.reset_auto_ids()
    material = openmc.Material(name="H1")
    material.add_element("H", 1.0)
    sphere = openmc.Sphere(r=1.0, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=material)
    return openmc.Geometry([cell])


@pytest.mark.parametrize(
    "parameter",
    [
        {"max_collisions": 200},
        {"max_collisions": 200, "reactions": ["(n,disappear)"]},
        {"max_collisions": 200, "cell_ids": [1]},
        {"max_collisions": 200, "material_ids": [1]},
        {"max_collisions": 200, "universe_ids": [1]},
        {"max_collisions": 200, "nuclides": ["H1"]},
        {"max_collisions": 200, "deposited_E_threshold": 200000.0},
        {"max_collisions": 200, "mcpl": True}

    ],
    ids=str
)
def test_xml_serialization(parameter, run_in_tmpdir):
    """Check that the different use cases can be written and read in XML."""
    settings = openmc.Settings()
    settings.collision_track = parameter
    settings.export_to_xml()

    read_settings = openmc.Settings.from_xml()
    assert read_settings.collision_track == parameter


@pytest.fixture
def model():
    """Simple hydrogen sphere divided in two hemispheres
    by a z-plane to form 2 cells."""
    openmc.reset_auto_ids()
    model = openmc.Model()

    # Material
    material = openmc.Material(name="H1")
    material.add_element("H", 1.0)

    # Geometry
    radius = 1.0
    sphere = openmc.Sphere(r=radius, boundary_type="reflective")
    plane = openmc.ZPlane(0.0)
    cell_1 = openmc.Cell(region=-sphere & -plane, fill=material, cell_id=1)
    cell_2 = openmc.Cell(region=-sphere & +plane, fill=material, cell_id=2)
    root = openmc.Universe(cells=[cell_1, cell_2])
    model.geometry = openmc.Geometry(root)

    # Settings
    model.settings = openmc.Settings()
    model.settings.run_mode = "fixed source"
    model.settings.particles = 1
    model.settings.batches = 1
    model.settings.seed = 2

    bounds = [-radius, -radius, -radius, radius, radius, radius]
    distribution = openmc.stats.Box(bounds[:3], bounds[3:])
    model.settings.source = openmc.IndependentSource(space=distribution)

    return model


def test_particle_location(run_in_tmpdir, model):
    """Test the location of particles with respected to the "cell_ids"
    and the location x, y, z of the particle itself. the upper sphere will
    have positive z component and the bottom sphere a negative z compnent.

    """
    model.settings.collision_track = {
        "max_collisions": 200,
        "reactions": ["elastic"],
        "cell_ids": [1, 2]
    }
    model.run()

    with h5py.File("collision_track.h5", "r") as f:
        source = f["collision_track_bank"]

        assert len(source) == 60

        # We want to verify that the collisions happenening are in the right cells
        # and the position of the particle is either positive or negative relative
        # to the z plane. In this case, we track the position of the particle
        # relative to the cell_id already set.
        for point in source:
            if point['cell_id'] == 1:
                assert point['r'][2] < 0.0  # z component negative
            elif point['cell_id'] == 2:
                assert point["r"][2] > 0.0  # z component positive
            else:
                assert False


@pytest.mark.skipif(shutil.which("mcpl-config") is None, reason="MCPL is not available.")
def test_format_similarity(run_in_tmpdir, model):
    model.settings.collision_track = {"max_collisions": 200, "reactions": ['elastic'],
                                      "cell_ids": [1, 2], "mcpl": False}
    model.run()
    data_h5 = ctt._return_collision_track_data('collision_track.h5')

    model.settings.collision_track["mcpl"] = True
    model.run()
    data_mcpl = ctt._return_collision_track_data('collision_track.mcpl')

    assert len(data_h5) == 60
    assert len(data_mcpl) == 60

    np.testing.assert_allclose(data_h5, data_mcpl, rtol=1e-05)
    # tolerance not that low due to the strings that is saved in MCPL,
    # not enough precision!


def test_photon_particles(run_in_tmpdir, model):
    """Test that the collision track can be used to track photon particles."""
    model.settings.collision_track = {"max_collisions": 200, "cell_ids": [1, 2]}

    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(*model.geometry.bounding_box),
        energy=openmc.stats.delta_function(1e5),
        particle='photon'
    )
    model.run()

    with h5py.File("collision_track.h5", "r") as f:
        source = f["collision_track_bank"]

        assert len(source) < 200

        allowed_particles = (openmc.ParticleType.PHOTON, openmc.ParticleType.ELECTRON)

        for point in source:
            particle_type = openmc.ParticleType(point['particle'])
            assert particle_type in allowed_particles

            if particle_type == openmc.ParticleType.ELECTRON:
                assert point['nuclide_id'] == 0


def test_collision_track_two_threads(model, run_in_tmpdir):
    # This test checks that the `max_collisions` setting is honored:
    # no collisions beyond the specified limit should be recorded.
    #
    # The exact set of events in the capped bank is not reproducible with
    # multiple threads because the bank stores whichever thread appends first
    # until capacity is reached.
    model.settings.collision_track = {"max_collisions": 200}
    model.run(threads=2, particles=500)

    collision_track = openmc.read_collision_track_hdf5("collision_track.h5")
    assert len(collision_track) == 200
