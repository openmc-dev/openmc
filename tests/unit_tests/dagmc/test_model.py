from pathlib import Path

import os
import h5py
import shutil
import lxml.etree as ET
import numpy as np
import pytest
import openmc
from openmc.utility_funcs import change_directory
import openmc.lib
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

    daguniv = openmc.DAGMCUniverse(p, auto_geom_ids=True)

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
            assert univ.cells[cell_id] == mats["foo"]


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
            assert univ.cells[cell_id] == mats["foo"]


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
            assert univ.cells[cell_id] == mats["foo"]


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
    assert dagmc_ele.get('filename') == str(dag_univ.filename)
    assert dagmc_ele.get('auto_geom_ids') == str(dag_univ.auto_geom_ids).lower()

    override_eles = dagmc_ele.find('material_overrides').findall('cell_override')
    assert len(override_eles) == 4

    for i, override_ele in enumerate(override_eles):
        cell_id = override_ele.get('id')
        assert dag_univ.material_overrides[int(cell_id)][0].id == int(override_ele.find('material_ids').text)

    model.export_to_model_xml()

    xml_model = openmc.Model.from_model_xml()

    for univ in xml_model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            xml_dagmc_univ = univ
            break

    assert xml_dagmc_univ._material_overrides.keys() == dag_univ._material_overrides.keys()

    for xml_mats, model_mats in zip(xml_dagmc_univ._material_overrides.values(), dag_univ._material_overrides.values()):
        assert all([xml_mat.id == orig_mat.id for xml_mat, orig_mat in zip(xml_mats, model_mats)])

def test_dagmc_vacuum(request, run_in_tmpdir):
    tmpdir = Path(os.getcwd())
    # Required initial mats.
    mats = {}
    mats["41"] = openmc.Material(name="41")
    mats["41"].add_nuclide("U235", 1)
    mats["41"].set_density("g/cm3", 1)

    # Not a vacuum material
    na_vacuum_str = "not_a_vacuum"
    mats[na_vacuum_str] = openmc.Material(name=na_vacuum_str)
    mats[na_vacuum_str].add_nuclide("U235", 1)
    mats[na_vacuum_str].set_density("g/cm3", 1)

    # Tweaking the h5m file to change the material assignment
    na_vacuum_h5 = tmpdir / f"dagmc_{na_vacuum_str}.h5m"
    shutil.copyfile(Path(request.fspath).parent / "dagmc.h5m",
                    na_vacuum_h5)
    hf = h5py.File(na_vacuum_h5, 'r+')
    new_assignment = 'mat:not_a_vacuum'.encode('utf-8')
    hf['/tstt/tags/NAME']['values'][1] = new_assignment
    hf.close()

    # Set the Model
    daguniv = openmc.DAGMCUniverse(na_vacuum_h5,
                                   auto_geom_ids=True).bounded_universe()
    settings = openmc.Settings()
    settings.batches = 10
    settings.particles = 100
    model = openmc.Model()
    model.materials = openmc.Materials(mats.values())
    model.geometry = openmc.Geometry(root=daguniv)
    model.settings = settings

    # Set the Tally
    model.init_lib()
    model.sync_dagmc_universes()
    for univ in model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            dag_univ = univ
            break
    my_tally = openmc.Tally()
    cell_filter = openmc.CellFilter([x.id for x in dag_univ.get_all_cells().values()])
    my_tally.scores = ['nu-fission']
    my_tally.filters = [cell_filter]
    model.tallies.append(my_tally)

    # Only a confined point source
    src = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.source = src

    # Run and check the results
    sub_dir = f"sub_dir_{na_vacuum_str}"
    model.export_to_model_xml(f'{sub_dir}/model.xml',)
    openmc.run(cwd=sub_dir)
    # If cell are not vacuum, fission site should not spread to other cells
    with openmc.StatePoint(f"{sub_dir}/statepoint.{settings.batches}.h5") as sp:
        tally = list(sp.tallies.values())[0]
        non_zero_tal = [x for x in tally.mean.flatten() if x > 0]
        assert len(non_zero_tal) > 1

    # Vacuum material
    vacuum_str = "vacuum"
    mats[vacuum_str] = openmc.Material(name=vacuum_str)
    mats[vacuum_str].add_nuclide("U235", 0.03)
    mats[vacuum_str].add_nuclide("U238", 0.97)
    mats[vacuum_str].add_nuclide("O16", 2.0)
    mats[vacuum_str].set_density("g/cm3", 10.0)

    # Tweaking the h5m file to change the material assignment
    vacuum_h5 = tmpdir / f"dagmc_{vacuum_str}.h5m"
    shutil.copyfile(Path(request.fspath).parent / "dagmc.h5m", vacuum_h5)
    hf = h5py.File(vacuum_h5, 'r+')
    new_assignment = 'mat:vacuum'.encode('utf-8')
    hf['/tstt/tags/NAME']['values'][1] = new_assignment
    hf.close()

    # Set the Model
    daguniv = openmc.DAGMCUniverse(vacuum_h5,
                                   auto_geom_ids=True).bounded_universe()
    settings = openmc.Settings()
    settings.batches = 10
    settings.particles = 100
    model = openmc.Model()
    model.materials = openmc.Materials(mats.values())
    model.geometry = openmc.Geometry(root=daguniv)
    model.settings = settings

    # Set the Tally
    model.init_lib()
    model.sync_dagmc_universes()
    for univ in model.geometry.get_all_universes().values():
        if isinstance(univ, openmc.DAGMCUniverse):
            dag_univ = univ
            break
    my_tally = openmc.Tally()
    cell_filter = openmc.CellFilter([x.id for x in dag_univ.get_all_cells().values()])
    my_tally.scores = ['nu-fission']
    my_tally.filters = [cell_filter]
    model.tallies.append(my_tally)

    # Only a confined point source
    src = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.source = src

    # Run and check the results
    sub_dir = f"sub_dir_{vacuum_str}"
    model.export_to_model_xml(f'{sub_dir}/model.xml',)
    openmc.run(cwd=sub_dir)
    # If cell are not vacuum, fission site should not spread to other cells
    with openmc.StatePoint(f"{sub_dir}/statepoint.{settings.batches}.h5") as sp:
        tally = list(sp.tallies.values())[0]
        non_zero_tal = [x for x in tally.mean.flatten() if x > 0]
        assert len(non_zero_tal) == 1
