import shutil

import numpy as np
import pytest

import openmc
import openmc.lib

from tests import cdtemp

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")


@pytest.fixture(scope="module", autouse=True)
def dagmc_model(request):

    model = openmc.model.Model()

    # settings
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 100
    model.settings.temperature = {'tolerance': 50.0}
    model.settings.verbosity = 1
    source_box = openmc.stats.Box([ -4, -4, -4 ],
                                  [  4,  4,  4 ])
    source = openmc.Source(space=source_box)
    model.settings.source = source

    model.settings.dagmc = True

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
    u235.temperature = 320

    water = openmc.Material(name="water")
    water.add_nuclide('H1', 2.0, 'ao')
    water.add_nuclide('O16', 1.0, 'ao')
    water.set_density('g/cc', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.id = 41

    mats = openmc.Materials([u235, water])
    model.materials = mats

    # location of  dagmc file in test directory
    dagmc_file = request.fspath.dirpath() + "/dagmc.h5m"
    # move to a temporary directory
    with cdtemp():
        shutil.copyfile(dagmc_file, "./dagmc.h5m")
        model.export_to_xml()
        openmc.lib.init()
        yield

    openmc.lib.finalize()


@pytest.mark.parametrize("cell_id,exp_temp", ((1, 320.0),   # assigned by material
                                              (2, 300.0),   # assigned in dagmc file
                                              (3, 293.6)))  # assigned by default
def test_dagmc_temperatures(cell_id, exp_temp):
    cell = openmc.lib.cells[cell_id]
    assert np.isclose(cell.get_temperature(), exp_temp)
