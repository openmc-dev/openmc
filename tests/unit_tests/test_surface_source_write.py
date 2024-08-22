"""Test the 'surf_source_write' setting used to store particles that cross
surfaces in a file for a given simulation."""

from pathlib import Path

import openmc
import openmc.lib
import pytest
import h5py
import numpy as np


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
        {"max_particles": 200},
        {"max_particles": 200, "cell": 1},
        {"max_particles": 200, "cellto": 1},
        {"max_particles": 200, "cellfrom": 1},
        {"max_particles": 200, "surface_ids": [2]},
        {"max_particles": 200, "surface_ids": [2], "cell": 1},
        {"max_particles": 200, "surface_ids": [2], "cellto": 1},
        {"max_particles": 200, "surface_ids": [2], "cellfrom": 1},
    ],
)
def test_xml_serialization(parameter, run_in_tmpdir):
    """Check that the different use cases can be written and read in XML."""
    settings = openmc.Settings()
    settings.surf_source_write = parameter
    settings.export_to_xml()

    read_settings = openmc.Settings.from_xml()
    assert read_settings.surf_source_write == parameter


ERROR_MSG_1 = (
    "A maximum number of particles needs to be specified "
    "using the 'max_particles' parameter to store surface "
    "source points."
)
ERROR_MSG_2 = "'cell', 'cellfrom' and 'cellto' cannot be used at the same time."


@pytest.mark.parametrize(
    "parameter, error",
    [
        ({"cell": 1}, ERROR_MSG_1),
        ({"max_particles": 200, "cell": 1, "cellto": 1}, ERROR_MSG_2),
        ({"max_particles": 200, "cell": 1, "cellfrom": 1}, ERROR_MSG_2),
        ({"max_particles": 200, "cellto": 1, "cellfrom": 1}, ERROR_MSG_2),
        ({"max_particles": 200, "cell": 1, "cellto": 1, "cellfrom": 1}, ERROR_MSG_2),
    ],
)
def test_exceptions(parameter, error, run_in_tmpdir, geometry):
    """Test parameters configuration that should return an error."""
    settings = openmc.Settings(run_mode="fixed source", batches=5, particles=100)
    settings.surf_source_write = parameter
    model = openmc.Model(geometry=geometry, settings=settings)
    with pytest.raises(RuntimeError, match=error):
        model.run()


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
    cell_1 = openmc.Cell(region=-sphere & -plane, fill=material)
    cell_2 = openmc.Cell(region=-sphere & +plane, fill=material)
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
        {"max_particles": 200, "cellto": 2, "surface_ids": [2]},
        {"max_particles": 200, "cellfrom": 2, "surface_ids": [2]},
    ],
)
def test_particle_direction(parameter, run_in_tmpdir, model):
    """Test the direction of particles with the 'cellfrom' and 'cellto' parameters
    on a simple model with only one surface of interest.

    Cell 2 is the upper hemisphere and surface 2 is the plane dividing the sphere
    into two hemispheres.

    """
    model.settings.surf_source_write = parameter
    model.run()
    with h5py.File("surface_source.h5", "r") as f:
        source = f["source_bank"]

        assert len(source) == 200

        # We want to verify that the dot product of the surface's normal vector
        # and the direction of the particle is either positive or negative
        # depending on cellfrom or cellto. In this case, it is equivalent
        # to just compare the z component of the direction of the particle.
        for point in source:
            if "cellto" in parameter.keys():
                assert point["u"]["z"] > 0.0
            elif "cellfrom" in parameter.keys():
                assert point["u"]["z"] < 0.0
            else:
                assert False


@pytest.fixture
def model_dagmc(request):
    """Model based on the mesh file 'dagmc.h5m' available from
    tests/regression_tests/dagmc/legacy.

    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    u235 = openmc.Material(name="no-void fuel")
    u235.add_nuclide("U235", 1.0, "ao")
    u235.set_density("g/cc", 11)
    u235.id = 40

    water = openmc.Material(name="water")
    water.add_nuclide("H1", 2.0, "ao")
    water.add_nuclide("O16", 1.0, "ao")
    water.set_density("g/cc", 1.0)
    water.add_s_alpha_beta("c_H_in_H2O")
    water.id = 41

    model.materials = openmc.Materials([u235, water])

    # =============================================================================
    # Geometry
    # =============================================================================
    dagmc_path = Path(request.fspath).parent / "../regression_tests/dagmc/legacy/dagmc.h5m"
    dagmc_univ = openmc.DAGMCUniverse(dagmc_path)
    model.geometry = openmc.Geometry(dagmc_univ)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.particles = 300
    model.settings.batches = 5
    model.settings.inactive = 1
    model.settings.seed = 1

    source_box = openmc.stats.Box([-4, -4, -20], [4, 4, 20])
    model.settings.source = openmc.IndependentSource(space=source_box)

    return model


@pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(), reason="DAGMC CAD geometry is not enabled."
)
@pytest.mark.parametrize(
    "parameter",
    [
        {"max_particles": 200, "cellto": 1},
        {"max_particles": 200, "cellfrom": 1},
    ],
)
def test_particle_direction_dagmc(parameter, run_in_tmpdir, model_dagmc):
    """Test the direction of particles with the 'cellfrom' and 'cellto' parameters
    on a DAGMC model."""
    model_dagmc.settings.surf_source_write = parameter
    model_dagmc.run()

    r = 7.0
    h = 20.0

    with h5py.File("surface_source.h5", "r") as f:
        source = f["source_bank"]

        assert len(source) == 200

        for point in source:

            x, y, z = point["r"]
            ux, uy, uz = point["u"]

            # If the point is on the upper or lower circle
            if np.allclose(abs(z), h):
                # If the point is also on the cylindrical surface
                if np.allclose(np.sqrt(x**2 + y**2), r):
                    if "cellfrom" in parameter.keys():
                        assert (uz * z > 0) or (ux * x + uy * y > 0)
                    elif "cellto" in parameter.keys():
                        assert (uz * z < 0) or (ux * x + uy * y < 0)
                    else:
                        assert False
                # If the point is not on the cylindrical surface
                else:
                    if "cellfrom" in parameter.keys():
                        assert uz * z > 0
                    elif "cellto" in parameter.keys():
                        assert uz * z < 0
                    else:
                        assert False
            # If the point is not on the upper or lower circle,
            # meaning it is on the cylindrical surface
            else:
                if "cellfrom" in parameter.keys():
                    assert ux * x + uy * y > 0
                elif "cellto" in parameter.keys():
                    assert ux * x + uy * y < 0
                else:
                    assert False
