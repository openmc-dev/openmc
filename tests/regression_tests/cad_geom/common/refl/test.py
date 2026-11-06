import openmc
from openmc.stats import Box

import pytest
from tests.testing_harness import PyAPITestHarness

test_libs = []
if openmc.lib.feature_enabled('dagmc'):
    test_libs.append('dagmc')
if openmc.lib.feature_enabled('xdg'):
    test_libs.append('xdg')

@pytest.fixture()
def refl_pincell():

    def _refl_pincell(interface='dagmc'):
        openmc.reset_auto_ids()

        model = openmc.Model()
        # settings
        model.settings.batches = 5
        model.settings.inactive = 0
        model.settings.particles = 100

        source = openmc.IndependentSource(space=Box([-4, -4, -4],
                                                    [ 4,  4,  4]))
        model.settings.source = source

        # geometry
        if interface == 'dagmc':
            cad_univ = openmc.DAGMCUniverse("dagmc.h5m", auto_geom_ids=True)
        elif interface == 'xdg':
            xdg_mesh = openmc.UnstructuredMesh('dagmc.h5m', library='moab')
            xdg_mesh.interface = 'xdg'
            cad_univ = openmc.XDGUniverse(xdg_mesh, auto_geom_ids=True)

        model.geometry = openmc.Geometry(cad_univ)

        # materials
        u235 = openmc.Material(name="fuel")
        u235.add_nuclide('U235', 1.0, 'ao')
        u235.set_density('g/cc', 11)
        u235.id = 40

        water = openmc.Material(name="moderator")
        water.add_nuclide('H1', 2.0, 'ao')
        water.add_nuclide('O16', 1.0, 'ao')
        water.set_density('g/cc', 1.0)
        water.add_s_alpha_beta('c_H_in_H2O')
        water.id = 41

        model.materials = openmc.Materials([u235, water])

        # tally
        tally = openmc.Tally()
        tally.scores = ['total']
        tally.filters = [openmc.CellFilter(2)]
        model.tallies = [tally]

        return model

    return _refl_pincell

@pytest.mark.parametrize("interface", test_libs)
def test_refl(refl_pincell, interface):
    input = f'inputs_{interface}.dat'
    harness = PyAPITestHarness('statepoint.5.h5', refl_pincell(interface), input)
    harness.main()
