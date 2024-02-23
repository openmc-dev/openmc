"""Test the meshborn filter using a fixed source calculation on a H1 sphere.

"""

import numpy as np
from uncertainties import unumpy
import openmc
import pytest


@pytest.fixture
def model():
    """Sphere of H1 with one hemisphere containing the source (x>0) and one
    hemisphere with no source (x<0).

    """
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    h1 = openmc.Material()
    h1.add_nuclide("H1", 1.0)
    h1.set_density("g/cm3", 1.0)

    model.materials = openmc.Materials([h1])

    # =============================================================================
    # Geometry
    # =============================================================================

    # Core geometry
    radius = 10.0
    sphere = openmc.Sphere(r=radius, boundary_type="reflective")
    core_region = -sphere
    core = openmc.Cell(fill=h1, region=core_region)

    # Root universe
    root = openmc.Universe(cells=[core])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 2000
    model.settings.batches = 8
    model.settings.photon_transport = False

    bounds = [0., -radius, -radius, radius, radius, radius]
    distribution = openmc.stats.Box(bounds[:3], bounds[3:])
    model.settings.source = openmc.source.IndependentSource(space=distribution)

    # =============================================================================
    # Tallies
    # =============================================================================

    mesh_1 = openmc.RegularMesh()
    mesh_1.dimension = [2,2,1]
    mesh_1.lower_left = [-radius, -radius, -radius]
    mesh_1.upper_right = [radius, radius, radius]

    f_1 = openmc.MeshBornFilter(mesh_1)

    t_1 = openmc.Tally(name="scatter-collision")
    t_1.filters = [f_1]
    t_1.scores = ["scatter"]
    t_1.estimator = "collision"

    mesh_2 = openmc.RegularMesh()
    mesh_2.dimension = [2,2,1]
    mesh_2.lower_left = [-radius, -radius, -radius]
    mesh_2.upper_right = [radius, radius, radius]

    f_2 = openmc.MeshBornFilter(mesh_2)

    t_2 = openmc.Tally(name="scatter-tracklength")
    t_2.filters = [f_2]
    t_2.scores = ["scatter"]
    t_2.estimator = "tracklength"

    model.tallies += [t_1, t_2]

    return model


def test_estimator_consistency(model, run_in_tmpdir):
    """Test that resuts obtained from a tracklength estimator are
    consistent with results obtained from a collision estimator.

    """
    # Run OpenMC
    sp_filename = model.run()

    # Get radial flux distribution
    with openmc.StatePoint(sp_filename) as sp:
        scatter_collision = sp.get_tally(name="scatter-collision").mean.ravel()
        scatter_collision_std_dev = sp.get_tally(name="scatter-collision").std_dev.ravel()
        scatter_tracklength = sp.get_tally(name="scatter-tracklength").mean.ravel()
        scatter_tracklength_std_dev = sp.get_tally(name="scatter-tracklength").std_dev.ravel()

    collision = unumpy.uarray(scatter_collision, scatter_collision_std_dev)
    tracklength = unumpy.uarray(scatter_tracklength, scatter_tracklength_std_dev)
    delta = abs(collision - tracklength)

    diff = unumpy.nominal_values(delta)
    std_dev = unumpy.std_devs(delta)
    assert np.all(diff <= 3 * std_dev)


def test_xml_serialization():
    """Test xml serialization of the meshborn filter."""
    openmc.reset_auto_ids()

    mesh = openmc.RegularMesh()
    mesh.dimension = [1, 1, 1]
    mesh.lower_left = [0.0, 0.0, 0.0]
    mesh.upper_right = [1.0, 1.0, 1.0]

    filter = openmc.MeshBornFilter(mesh)
    filter.translation = [2.0 ,2.0 ,2.0]
    assert filter.mesh.id == 1
    assert filter.mesh.dimension == (1, 1 ,1)
    assert filter.mesh.lower_left == [0.0, 0.0, 0.0]
    assert filter.mesh.upper_right == [1.0 ,1.0 ,1.0]

    repr(filter)

    elem = filter.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'meshborn'
    assert elem[0].text == "1"
    assert elem.get("translation") == "2.0 2.0 2.0"

    meshes = {1: mesh}
    new_filter = openmc.Filter.from_xml_element(elem, meshes=meshes)
    assert new_filter.bins == filter.bins
    np.testing.assert_equal(new_filter.translation, [2.0, 2.0, 2.0])
