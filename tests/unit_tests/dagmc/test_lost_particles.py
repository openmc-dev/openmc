import numpy as np
from pathlib import Path

import openmc
import openmc.lib

import pytest

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")


@pytest.fixture
def broken_dagmc_model(request):
    openmc.reset_auto_ids()
    model = openmc.Model()

    ### MATERIALS ###
    fuel = openmc.Material(name='no-void fuel')
    fuel.set_density('g/cc', 10.29769)
    fuel.add_nuclide('U233', 1.0)

    cladding = openmc.Material(name='clad')
    cladding.set_density('g/cc', 6.55)
    cladding.add_nuclide('Zr90', 1.0)

    h1 = openmc.Material(name='water')
    h1.set_density('g/cc', 0.75)
    h1.add_nuclide('H1', 1.0)

    model.materials = openmc.Materials([fuel, cladding, h1])

    ### GEOMETRY ###
    # create the DAGMC universe using a model that has many triangles
    # removed
    dagmc_file = Path(request.fspath).parent / "broken_model.h5m"
    pincell_univ = openmc.DAGMCUniverse(filename=dagmc_file, auto_geom_ids=True)

    # create a 2 x 2 lattice using the DAGMC pincell
    pitch = np.asarray((24.0, 24.0))
    lattice = openmc.RectLattice()
    lattice.pitch = pitch
    lattice.universes = [[pincell_univ] * 2] * 2
    lattice.lower_left = -pitch

    # clip the DAGMC geometry at +/- 10 cm w/ CSG planes
    rpp = openmc.model.RectangularParallelepiped(
        -pitch[0], pitch[0], -pitch[1], pitch[1], -10.0, 10.0, boundary_type='reflective')
    bounding_cell = openmc.Cell(fill=lattice, region=-rpp)

    model.geometry = openmc.Geometry(root=[bounding_cell])

    # settings
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.inactive = 2
    model.settings.output = {'summary': False}

    model.export_to_xml()

    return model


def test_lost_particles(broken_dagmc_model):
    broken_dagmc_model.export_to_xml()
    # ensure that particles will be lost when cell intersections can't be found
    # due to the removed triangles in this model
    with pytest.raises(RuntimeError, match='Maximum number of lost particles has been reached.'):
        openmc.run()

    # run this again, but with the dagmc universe as the root unvierse
    # to ensure that lost particles are still caught in this case
    for univ in broken_dagmc_model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            broken_dagmc_model.geometry.root_universe = univ
            break

    broken_dagmc_model.export_to_xml()
    with pytest.raises(RuntimeError, match='Maximum number of lost particles has been reached.'):
        openmc.run()
