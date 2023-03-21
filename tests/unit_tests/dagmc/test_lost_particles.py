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
    model = openmc.model.Model()

    ### MATERIALS ###
    fuel = openmc.Material(name='no-void fuel')
    fuel.set_density('g/cc', 10.29769)
    fuel.add_nuclide('U234', 0.93120485)
    fuel.add_nuclide('U235', 0.00055815)
    fuel.add_nuclide('U238', 0.022408)
    fuel.add_nuclide('O16', 0.045829)

    cladding = openmc.Material(name='clad')
    cladding.set_density('g/cc', 6.55)
    cladding.add_nuclide('Zr90', 0.021827)
    cladding.add_nuclide('Zr91', 0.00476)
    cladding.add_nuclide('Zr92', 0.0072758)
    cladding.add_nuclide('Zr94', 0.0073734)
    cladding.add_nuclide('Zr96', 0.0011879)

    water = openmc.Material(name='water')
    water.set_density('g/cc', 0.740582)
    water.add_nuclide('H1', 0.049457)
    water.add_nuclide('O16', 0.024672)
    water.add_nuclide('B10', 8.0042e-06)
    water.add_nuclide('B11', 3.2218e-05)
    water.add_s_alpha_beta('c_H_in_H2O')

    model.materials = openmc.Materials([fuel, cladding, water])

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

    left = openmc.XPlane(x0=-pitch[0], name='left', boundary_type='reflective')
    right = openmc.XPlane(x0=pitch[0], name='right', boundary_type='reflective')
    front = openmc.YPlane(y0=-pitch[1], name='front', boundary_type='reflective')
    back = openmc.YPlane(y0=pitch[1], name='back', boundary_type='reflective')
    # clip the DAGMC geometry at +/- 10 cm w/ CSG planes
    bottom = openmc.ZPlane(z0=-10.0, name='bottom', boundary_type='reflective')
    top = openmc.ZPlane(z0=10.0, name='top', boundary_type='reflective')

    bounding_region = +left & -right & +front & -back & +bottom & -top
    bounding_cell = openmc.Cell(fill=lattice, region=bounding_region)

    model.geometry = openmc.Geometry(root=[bounding_cell])

    # settings
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.inactive = 2
    model.settings.output = {'summary' : False}

    model.export_to_xml()

    return model


def test_lost_particles(run_in_tmpdir, broken_dagmc_model):
    broken_dagmc_model.export_to_xml()
    # ensure that particles will be lost when cell intersections can't be found
    # due to the removed triangles in this model
    with pytest.raises(RuntimeError, match='Maximum number of lost particles has been reached.'):
        openmc.run()

    # run this again, but with the dagmc universe as the root unvierse
    for univ in broken_dagmc_model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            broken_dagmc_model.geometry.root_unvierse = univ
            break

    broken_dagmc_model.export_to_xml()
    with pytest.raises(RuntimeError, match='Maximum number of lost particles has been reached.'):
        openmc.run()

