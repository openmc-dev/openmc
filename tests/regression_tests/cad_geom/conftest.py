from pathlib import Path

import openmc
import pytest


@pytest.fixture()
def cad_pincell_model():

    def _pincell_model(interface='dagmc'):
        openmc.reset_auto_ids()

        model = openmc.Model()

        # settings
        model.settings.batches = 5
        model.settings.inactive = 0
        model.settings.particles = 100

        source_box = openmc.stats.Box([-4, -4, -4],
                                    [ 4,  4,  4])
        source = openmc.IndependentSource(space=source_box)

        model.settings.source = source

        # geometry
        if interface == 'dagmc':
            dagmc_univ = openmc.DAGMCUniverse(Path("dagmc.h5m"))
            model.geometry = openmc.Geometry(dagmc_univ)
        elif interface == 'xdg':
            mesh = openmc.UnstructuredMesh("dagmc.h5m", library='moab')
            mesh.interface = 'xdg'
            xdg_univ = openmc.XDGUniverse(mesh)
            model.geometry = openmc.Geometry(xdg_univ)

        # tally
        tally = openmc.Tally()
        tally.scores = ['total']
        tally.filters = [openmc.CellFilter(1)]
        model.tallies = [tally]

        # materials
        u235 = openmc.Material(name="no-void fuel")
        u235.add_nuclide('U235', 1.0, 'ao')
        u235.set_density('g/cc', 11)
        u235.id = 40

        water = openmc.Material(name="water")
        water.add_nuclide('H1', 2.0, 'ao')
        water.add_nuclide('O16', 1.0, 'ao')
        water.set_density('g/cc', 1.0)
        water.add_s_alpha_beta('c_H_in_H2O')
        water.id = 41

        mats = openmc.Materials([u235, water])
        model.materials = mats

        return model

    return _pincell_model
