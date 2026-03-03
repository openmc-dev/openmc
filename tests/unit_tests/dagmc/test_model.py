from pathlib import Path

import lxml.etree as ET
import numpy as np
import pytest
import openmc
import openmc.lib
from openmc.utility_funcs import change_directory

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")


@pytest.fixture()
def model(request):
    pitch = 1.26

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

    p = Path(request.fspath).parent / "dagmc.h5m"

    daguniv = openmc.DAGMCUniverse(p, name='simple-dagmc', auto_geom_ids=True)

    lattice = openmc.RectLattice()
    lattice.dimension = [2, 2]
    lattice.lower_left = [-pitch, -pitch]
    lattice.pitch = [pitch, pitch]
    lattice.universes = [
        [daguniv, daguniv],
        [daguniv, daguniv]]

    box = openmc.model.RectangularParallelepiped(-pitch, pitch, -pitch, pitch, -5, 5)

    root = openmc.Universe(cells=[openmc.Cell(region=-box, fill=lattice)])

    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 10
    settings.particles = 1000

    ll, ur = root.bounding_box
    mat_vol = openmc.VolumeCalculation([mats["no-void fuel"]], 1000000, ll, ur)
    cell_vol = openmc.VolumeCalculation(list(root.cells.values()), 1000000, ll, ur)
    settings.volume_calculations = [mat_vol, cell_vol]

    model = openmc.Model()
    model.materials = openmc.Materials(mats.values())
    model.geometry = openmc.Geometry(root=root)
    model.settings = settings

    with change_directory(tmpdir=True):
        try:
            model.init_lib()
            model.sync_dagmc_universes()
            yield model
        finally:
            model.finalize_lib()
            openmc.reset_auto_ids()


def test_dagmc_replace_material_assignment(model):
    mats = {}

    mats["foo"] = openmc.Material(name="foo")
    mats["foo"].add_nuclide("H1", 2.0)
    mats["foo"].add_element("O", 1.0)
    mats["foo"].set_density("g/cm3", 1.0)
    mats["foo"].add_s_alpha_beta("c_H_in_H2O")

    for univ in model.geometry.get_all_universes().values():
        if not isinstance(univ, openmc.DAGMCUniverse):
            break

        cells_with_41 = []
        for cell in univ.cells.values():
            if cell.fill is None:
                continue
            if cell.fill.name == "41":
                cells_with_41.append(cell.id)
        univ.replace_material_assignment("41", mats["foo"])
        for cell_id in cells_with_41:
            assert univ.cells[cell_id].fill == mats["foo"]


def test_dagmc_sync_cell_names(model):
    dag_univ = None
    for univ in model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            dag_univ = univ
            break

    assert dag_univ is not None

    for cell_id, cell in dag_univ.cells.items():
        assert cell.name == openmc.lib.cells[cell_id].name

    assert any(cell.name == "implicit complement"
               for cell in dag_univ.cells.values())


def test_dagmc_add_material_override_with_id(model):
    mats = {}
    mats["foo"] = openmc.Material(name="foo")
    mats["foo"].add_nuclide("H1", 2.0)
    mats["foo"].add_element("O", 1.0)
    mats["foo"].set_density("g/cm3", 1.0)
    mats["foo"].add_s_alpha_beta("c_H_in_H2O")

    for univ in model.geometry.get_all_universes().values():
        if not isinstance(univ, openmc.DAGMCUniverse):
            break

        cells_with_41 = []
        for cell in univ.cells.values():
            if cell.fill is None:
                continue
            if cell.fill.name == "41":
                cells_with_41.append(cell.id)
                univ.add_material_override(cell.id, mats["foo"])
        for cell_id in cells_with_41:
            assert univ.cells[cell_id].fill == mats["foo"]


def test_dagmc_add_material_override_with_cell(model):
    mats = {}
    mats["foo"] = openmc.Material(name="foo")
    mats["foo"].add_nuclide("H1", 2.0)
    mats["foo"].add_element("O", 1.0)
    mats["foo"].set_density("g/cm3", 1.0)
    mats["foo"].add_s_alpha_beta("c_H_in_H2O")

    for univ in model.geometry.get_all_universes().values():
        if not isinstance(univ, openmc.DAGMCUniverse):
            break

        cells_with_41 = []
        for cell in univ.cells.values():
            if cell.fill is None:
                continue
            if cell.fill.name == "41":
                cells_with_41.append(cell.id)
                univ.add_material_override(cell, mats["foo"])
        for cell_id in cells_with_41:
            assert univ.cells[cell_id].fill == mats["foo"]


def test_model_differentiate_depletable_with_dagmc(model, run_in_tmpdir):
    model.calculate_volumes()

    # Get the volume of the no-void fuel material before differentiation
    volume_before = np.sum([m.volume for m in model.materials if m.name == "no-void fuel"])

    # Differentiate the depletable materials
    model.differentiate_depletable_mats(diff_volume_method="divide equally")
    # Get the volume of the no-void fuel material after differentiation
    volume_after = np.sum([m.volume for m in model.materials if "fuel" in m.name])
    assert np.isclose(volume_before, volume_after)
    assert len(model.materials) == 4*2 +1


def test_model_differentiate_with_dagmc(model):
    root = model.geometry.root_universe
    ll, ur = root.bounding_box
    model.calculate_volumes()
    # Get the volume of the no-void fuel material before differentiation
    volume_before = np.sum([m.volume for m in model.materials if m.name == "no-void fuel"])

    # Differentiate all the materials
    model.differentiate_mats(depletable_only=False)

    # Get the volume of the no-void fuel material after differentiation
    mat_vol = openmc.VolumeCalculation(model.materials, 1000000, ll, ur)
    model.settings.volume_calculations = [mat_vol]
    model.init_lib()  # need to reinitialize the lib after differentiating the materials
    model.calculate_volumes()
    volume_after = np.sum([m.volume for m in model.materials if "fuel" in m.name])
    assert np.isclose(volume_before, volume_after)
    assert len(model.materials) == 4*2 + 4


def test_bad_override_cell_id(model):
    for univ in model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            break
    with pytest.raises(ValueError, match="Cell ID '1' not found in DAGMC universe"):
        univ.material_overrides = {1: model.materials[0]}


def test_bad_override_type(model):
    not_a_dag_cell = openmc.Cell()
    for univ in model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            break
    with pytest.raises(ValueError, match="Unrecognized key type. Must be an integer or openmc.DAGMCCell object"):
        univ.material_overrides = {not_a_dag_cell: model.materials[0]}


def test_bad_replacement_mat_name(model):
    for univ in model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            break
    with pytest.raises(ValueError, match="No material with name 'not_a_mat' found in the DAGMC universe"):
        univ.replace_material_assignment("not_a_mat", model.materials[0])


def test_dagmc_xml(model):
    # Set the environment
    mats = {}
    mats["no-void fuel"] = openmc.Material(1, name="no-void fuel")
    mats["no-void fuel"].add_nuclide("U235", 0.03)
    mats["no-void fuel"].add_nuclide("U238", 0.97)
    mats["no-void fuel"].add_nuclide("O16", 2.0)
    mats["no-void fuel"].set_density("g/cm3", 10.0)

    mats[5] = openmc.Material(name="41")
    mats[5].add_nuclide("H1", 2.0)
    mats[5].add_element("O", 1.0)
    mats[5].set_density("g/cm3", 1.0)
    mats[5].add_s_alpha_beta("c_H_in_H2O")

    for univ in model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            dag_univ = univ
            break

    for k, v in mats.items():
        if isinstance(k, int):
            dag_univ.add_material_override(k, v)
            model.materials.append(v)
        elif isinstance(k, str):
            dag_univ.replace_material_assignment(k, v)

    # Tesing the XML subelement generation
    root = ET.Element('dagmc_universe')
    dag_univ.create_xml_subelement(root)
    dagmc_ele = root.find('dagmc_universe')

    assert dagmc_ele.get('id') == str(dag_univ.id)
    assert dagmc_ele.get('name') == str(dag_univ.name)
    assert dagmc_ele.get('filename') == str(dag_univ.filename)
    assert dagmc_ele.get('auto_geom_ids') == str(dag_univ.auto_geom_ids).lower()

    assert dagmc_ele.find('material_overrides') is None

    override_elements = dagmc_ele.findall('cell')
    assert len(override_elements) == len(dag_univ.cells)
    xml_cells = {int(elem.get('id')): elem for elem in override_elements}
    for cell_id, cell in dag_univ.cells.items():
        assert cell_id in xml_cells
        xml_cell = xml_cells[cell_id]
        if cell.fill_type == 'void':
            assert xml_cell.get('material') == 'void'
        elif cell.fill_type == 'material':
            assert xml_cell.get('material') == str(cell.fill.id)
        elif cell.fill_type == 'distribmat':
            mat_list = xml_cell.find('material').text.split()
            expected = ["void" if m is None else str(m.id) for m in cell.fill]
            assert mat_list == expected
        else:
            pytest.fail(f"Unexpected DAGMC cell fill type: {cell.fill_type}")

    model.export_to_model_xml()

    xml_model = openmc.Model.from_model_xml()

    for univ in xml_model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            xml_dagmc_univ = univ
            break

    assert xml_dagmc_univ.cells.keys() == dag_univ.cells.keys()

    for cell_id, cell in dag_univ.cells.items():
        xml_cell = xml_dagmc_univ.cells[cell_id]
        assert xml_cell.fill_type == cell.fill_type
        if cell.fill_type == 'void':
            assert xml_cell.fill is None
        elif cell.fill_type == 'material':
            assert xml_cell.fill.id == cell.fill.id
        elif cell.fill_type == 'distribmat':
            xml_ids = [m.id if m is not None else None for m in xml_cell.fill]
            model_ids = [m.id if m is not None else None for m in cell.fill]
            assert xml_ids == model_ids
        else:
            pytest.fail(f"Unexpected DAGMC cell fill type: {cell.fill_type}")


def test_dagmc_xml_reject_fill_override():
    mats = {'1': openmc.Material(1), 'void': None}
    elem = ET.fromstring(
        '<dagmc_universe id="1" filename="dagmc.h5m">'
        '<cell id="1" fill="2"/>'
        '</dagmc_universe>'
    )
    with pytest.raises(ValueError, match="only support material fills"):
        openmc.DAGMCUniverse.from_xml_element(elem, mats)


def test_dagmc_xml_reject_region_override():
    mats = {'1': openmc.Material(1), 'void': None}
    elem = ET.fromstring(
        '<dagmc_universe id="1" filename="dagmc.h5m">'
        '<cell id="1" material="1" region="-1"/>'
        '</dagmc_universe>'
    )
    with pytest.raises(ValueError, match="cannot include a region"):
        openmc.DAGMCUniverse.from_xml_element(elem, mats)


def test_dagmc_xml_reject_legacy_material_overrides():
    mats = {'1': openmc.Material(1), 'void': None}
    elem = ET.fromstring(
        '<dagmc_universe id="1" filename="dagmc.h5m">'
        '<material_overrides>'
        '<cell_override id="1"><material_ids>1</material_ids></cell_override>'
        '</material_overrides>'
        '</dagmc_universe>'
    )
    with pytest.raises(ValueError, match="no longer supported"):
        openmc.DAGMCUniverse.from_xml_element(elem, mats)


def test_dagmc_xml_temperature_roundtrip():
    mat = openmc.Material(1)
    mats = {'1': mat, 'void': None}

    elem = ET.fromstring(
        '<dagmc_universe id="10" filename="dagmc.h5m">'
        '<cell id="7" material="1" temperature="825.0"/>'
        '</dagmc_universe>'
    )

    dag_univ = openmc.DAGMCUniverse.from_xml_element(elem, mats)
    assert dag_univ.cells[7].fill.id == 1
    assert dag_univ.cells[7].temperature == pytest.approx(825.0)

    root = ET.Element('geometry')
    dag_univ.create_xml_subelement(root)
    dagmc_elem = root.find('dagmc_universe')
    xml_cell = dagmc_elem.find('cell')
    assert xml_cell.get('temperature') == '825.0'

    dag_univ_roundtrip = openmc.DAGMCUniverse.from_xml_element(dagmc_elem, mats)
    assert dag_univ_roundtrip.cells[7].fill.id == 1
    assert dag_univ_roundtrip.cells[7].temperature == pytest.approx(825.0)
