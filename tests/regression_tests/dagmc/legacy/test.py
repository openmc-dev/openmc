import openmc
import openmc.lib

from pathlib import Path
import pytest
from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")

@pytest.fixture
def model():
    openmc.reset_auto_ids()

    model = openmc.model.Model()

    # settings
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 100

    source_box = openmc.stats.Box([-4, -4, -4],
                                  [ 4,  4,  4])
    source = openmc.IndependentSource(space=source_box)

    model.settings.source = source

    # geometry
    dag_univ = openmc.DAGMCUniverse(Path("dagmc.h5m"))
    model.geometry = openmc.Geometry(dag_univ)

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


def test_missing_material_id(model):
    # remove the last material, which is identified by ID in the DAGMC file
    model.materials = model.materials[:-1]
    with pytest.raises(RuntimeError) as exec_info:
        model.run()
    exp_error_msg = "Material with name/ID '41' not found for volume (cell) 3"
    assert exp_error_msg in str(exec_info.value)


def test_missing_material_name(model):
    # remove the first material, which is identified by name in the DAGMC file
    model.materials = model.materials[1:]
    with pytest.raises(RuntimeError) as exec_info:
        model.run()
    exp_error_msg = "Material with name/ID 'no-void fuel' not found for volume (cell) 1"
    assert exp_error_msg in str(exec_info.value)


def test_dagmc(model):
    harness = PyAPITestHarness('statepoint.5.h5', model)
    harness.main()
