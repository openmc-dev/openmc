import pkg_resources
import shutil

import lxml.etree as ET
import numpy as np
from pathlib import Path
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
    source = openmc.IndependentSource(space=source_box)
    model.settings.source = source

    # geometry
    dagmc_file = Path(request.fspath).parent / 'dagmc.h5m'
    dagmc_universe = openmc.DAGMCUniverse(dagmc_file)
    model.geometry = openmc.Geometry(dagmc_universe)

    # check number of surfaces and volumes for this pincell model there should
    # be 5 volumes: two fuel regions, water, graveyard, implicit complement (the
    # implicit complement cell is created automatically at runtime)
    # and 21 surfaces: 3 cylinders (9 surfaces) and a bounding cubic shell
    # (12 surfaces)
    assert dagmc_universe.n_cells == 5
    assert dagmc_universe.n_surfaces == 21

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


def test_dagmc_xml():

    # Set the environment
    mats = {}
    mats["no-void fuel"] = openmc.Material(1, name="no-void fuel")
    mats["no-void fuel"].add_nuclide("U235", 0.03)
    mats["no-void fuel"].add_nuclide("U238", 0.97)
    mats["no-void fuel"].add_nuclide("O16", 2.0)
    mats["no-void fuel"].set_density("g/cm3", 10.0)

    mats["41"] = openmc.Material(name="41")
    mats["41"].add_nuclide("H1", 2.0)
    mats["41"].add_element("O", 1.0)
    mats["41"].set_density("g/cm3", 1.0)
    mats["41"].add_s_alpha_beta("c_H_in_H2O")

    p = pkg_resources.resource_filename(__name__, "dagmc.h5m")
    daguniv = openmc.DAGMCUniverse(p,auto_geom_ids=True,)
    daguniv._material_overrides = {40: [mats["no-void fuel"]], 52: [mats["41"]]}

    root = ET.Element('dagmc_universe')
    daguniv.create_xml_subelement(root)
    dagmc_ele = root.find('dagmc_universe')

    assert dagmc_ele.get('id') == str(daguniv.id)
    assert dagmc_ele.get('filename') == str(daguniv.filename)
    assert dagmc_ele.get('auto_geom_ids') == str(daguniv.auto_geom_ids).lower()

    # Verify the material overrides element in the XML tree
    override_eles = dagmc_ele.find('material_overrides').findall('cell')
    assert len(override_eles) == 2
    assert override_eles[0].get('id') == '40'
    assert override_eles[0].get('material') == str(mats["no-void fuel"].id)
    assert override_eles[1].get('id') == '52'
    assert override_eles[1].get('material') == str(mats["41"].id)

    # Create a dictionary of materials to pass to from_xml_element indexing by
    # material ID strings
    dict_mats = {str(m.id): m for m in mats.values()}
    dagmcuni_2 = openmc.DAGMCUniverse.from_xml_element(dagmc_ele, dict_mats)
    assert daguniv.id == dagmcuni_2.id
    assert daguniv.filename == dagmcuni_2.filename
    assert daguniv.auto_geom_ids == dagmcuni_2.auto_geom_ids
    assert daguniv._material_overrides == dagmcuni_2._material_overrides
