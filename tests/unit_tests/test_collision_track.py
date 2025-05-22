"""Test the 'collision_track' setting used to store particle information 
during specified collision conditions in a file for a given simulation."""

from pathlib import Path

import openmc
import openmc.lib
import pytest
import h5py
import numpy as np


def return_collision_track_data(filepath):
    """Read a collision_track file and return a sorted array composed
    of flatten arrays of collision information.

    Parameters
    ----------
    filepath : str
        Path to the collision_track file

    Returns
    -------
    data : np.array
        Sorted array composed of flatten arrays of source data for
        each collision information

    """
    data = []
    keys = []

    # Read source file
    source = openmc.read_collision_track_file(filepath)
    N = len(source)
    for j in range(N):
        r = source['r'][j]
        u = source['u'][j]
        e = source['E'][j]
        de = source['dE'][j]
        time = source['time'][j]
        wgt = source['wgt'][j]
        delayed_group = source['delayed_group'][j]
        cell_id = source['cell_id'][j]
        nuclide_id = source['nuclide_id'][j]
        material_id = source['material_id'][j]
        univ_id = source['universe_id'][j]
        event_mt = source['event_mt'][j]
        # no particle type because the code for MCPL is different from h5
        key = (
            f"{r[0]:.10e} {r[1]:.10e} {r[2]:.10e} {u[0]:.10e} {u[1]:.10e} {u[2]:.10e}"
            f"{e:.10e} {de:.10e}  {time:.10e} {wgt:.10e} {event_mt} {delayed_group} {cell_id}"
            f"{nuclide_id} {material_id} {univ_id} "
        )
        keys.append(key)
        values = [*r, *u, e, de, time, wgt, event_mt,
                  delayed_group, cell_id, nuclide_id, material_id, univ_id,]
        assert len(values) == 16
        data.append(values)

    data = np.array(data)
    keys = np.array(keys)
    sorted_idx = np.argsort(keys)

    return data[sorted_idx]


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
        {"max_collisions": 200, "mt_numbers": [1]},
        {"max_collisions": 200, "cell_ids": [1]},
        {"max_collisions": 200, "material_ids": [1]},
        {"max_collisions": 200, "universe_ids": [1]},
        {"max_collisions": 200, "nuclide_ids": [1001]},
        {"max_collisions": 200, "deposited_E_threshold": 200000},
        {"max_collisions": 200, "mcpl": True}

    ],
)
def test_xml_serialization(parameter, run_in_tmpdir):
    """Check that the different use cases can be written and read in XML."""
    settings = openmc.Settings()
    settings.collision_track = parameter
    settings.export_to_xml()

    read_settings = openmc.Settings.from_xml()
    assert read_settings.collision_track == parameter


@pytest.fixture(scope="module")
def model():
    """Simple hydrogen sphere geometry"""
    openmc.reset_auto_ids()
    model = openmc.Model()

    # Material
    h1 = openmc.Material(name="H1")
    h1.add_nuclide("H1", 1.0)
    h1.set_density('g/cm3', 1e-3)

    # Geometry
    radius = 1.0
    sphere = openmc.Sphere(r=radius, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere, fill=h1)
    model.geometry = openmc.Geometry([cell])

    # Settings
    model.settings = openmc.Settings()
    model.settings.run_mode = "fixed source"
    model.settings.particles = 100
    model.settings.batches = 3
    model.settings.seed = 1

    distribution = openmc.stats.Point()
    model.settings.source = openmc.IndependentSource(space=distribution)
    return model


@pytest.fixture(scope="module")
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
    model.settings.particles = 100
    model.settings.batches = 3
    model.settings.seed = 1

    bounds = [-radius, -radius, -radius, radius, radius, radius]
    distribution = openmc.stats.Box(bounds[:3], bounds[3:])
    model.settings.source = openmc.IndependentSource(space=distribution)

    return model


@pytest.mark.parametrize(
    "parameter",
    [
        {"max_collisions": 200, "mt_numbers": [2], "cell_ids": [1, 2]},

    ],
)
def test_particle_location(parameter, run_in_tmpdir, model):
    """Test the location of particles with respected to the "cell_ids" 
    and the location x, y, z of the particle itself. the upper sphere will
    have positive z component and the bottom sphere a negative z compnent.

    """
    model.settings.collision_track = parameter
    model.run()

    with h5py.File("collision_track.h5", "r") as f:
        source = f["collision_track_bank"]

        assert len(source) == 200

        # We want to verify that the collisions happenening are in the right cells
        # and the position of the particle is either positive or negative relative
        # to the z plane. In this case, we track the position of the particle
        # relative to the cell_id already set.
        for point in source:
            if point['cell_id'] == 1:
                assert point['r'][2] < 0.0
            elif point['cell_id'] == 2:
                assert point["r"][2] > 0.0
            else:
                assert False


@pytest.mark.skipif(
    not openmc.lib._mcpl_enabled(), reason="MCPL format is not enabled."
)
def test_format_similarity(run_in_tmpdir, model):

    model.settings.collision_track = {"max_collisions": 200, "mt_numbers": [2],
                                      "cell_ids": [1, 2], "mcpl": False}

    model.run(threads=1)

    data_h5 = return_collision_track_data('collision_track.h5')

    model.settings.collision_track = {"max_collisions": 200, "mt_numbers": [2],
                                      "cell_ids": [1, 2], "mcpl": True}

    model.run(threads=1)

    data_mcpl = return_collision_track_data('collision_track.mcpl')

    assert len(data_h5) == 200
    assert len(data_mcpl) == 200

    np.testing.assert_allclose(data_h5, data_mcpl, rtol=1e-05)
    # tolerance not that low due to the strings that is saved in MCPL,
    # not enough precision!
