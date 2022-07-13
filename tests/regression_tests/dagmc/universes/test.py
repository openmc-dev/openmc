import numpy as np

import openmc
import openmc.lib

import pytest
from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")


class DAGMCUniverseTest(PyAPITestHarness):

    def _build_inputs(self):
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
        # create the DAGMC universe
        pincell_univ = openmc.DAGMCUniverse(filename='dagmc.h5m', auto_geom_ids=True)

        # checks that the bounding box is calculated correctly
        bounding_box = pincell_univ.bounding_box
        assert bounding_box[0].tolist() == [-25., -25., -25.]
        assert bounding_box[1].tolist() == [25., 25., 25.]

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

        model.geometry = openmc.Geometry([bounding_cell])

        # settings
        model.settings.particles = 100
        model.settings.batches = 10
        model.settings.inactive = 2
        model.settings.output = {'summary' : False}

        model.export_to_xml()


def test_univ():
    harness = DAGMCUniverseTest('statepoint.10.h5')
    harness.main()
