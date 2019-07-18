import glob
import os

import numpy as np
import pytest

import openmc
import openmc.capi

from tests import cdtemp

pytestmark = pytest.mark.skipif(
    not openmc.capi._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")

def test_dagmc_temperatures():
    model = openmc.model.Model()

    # settings
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 100
    model.settings.temperature = {'tolerance': 50.0}
    model.settings.verbosity = 1
    source_box = openmc.stats.Box([-4, -4, -4],
                                  [ 4,  4,  4])
    source = openmc.Source(space=source_box)
    model.settings.source = source

    model.settings.dagmc = True

    # tally
    tally = openmc.Tally()
    tally.scores = ['total']
    tally.filters = [openmc.CellFilter(1)]
    model.tallies = [tally]

    # materials
    u235 = openmc.Material(name="fuel")
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

    model.export_to_xml()

    # check cell temps as well here
    openmc.capi.init()

    expected_temps = { 1 : 320.0,  # assigned by material
                       2 : 300.0,  # assigned in dagmc file
                       3 : 293.6 } # assigned by default

    for cell_id, temp in expected_temps.items():
        capi_cell = openmc.capi.cells[cell_id]
        assert np.isclose(capi_cell.get_temperature(), temp)

    openmc.capi.finalize()

    # cleanup
    input_files = glob.glob("*.xml")
    input_files += glob.glob("*.h5")
    for f in input_files:
        if os.path.exists(f):
            os.remove(f)
