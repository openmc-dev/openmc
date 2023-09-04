""" Tests for Batchwise class """

from pathlib import Path
from math import exp

import pytest
import numpy as np

import openmc
import openmc.lib
from openmc.deplete import CoupledOperator
from openmc.deplete import (BatchwiseCellGeometrical, BatchwiseCellTemperature,
    BatchwiseMaterialRefuel)

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"

@pytest.fixture
def model():
    f = openmc.Material(name="f")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)
    f.temperature = 293.15
    w = openmc.Material(name="w")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.depletable = True
    h = openmc.Material(name='h')
    h.set_density('g/cm3', 0.001598)
    h.add_element('He', 2.4044e-4)
    radii = [0.42, 0.45]
    f.volume = np.pi * radii[0] ** 2
    w.volume = np.pi * (radii[1]**2 - radii[0]**2)
    materials = openmc.Materials([f, w, h])
    surf_wh = openmc.ZPlane(z0=0)
    surf_f = openmc.Sphere(r=radii[0])
    surf_w = openmc.Sphere(r=radii[1], boundary_type='vacuum')
    cell_w = openmc.Cell(fill=w, region=-surf_wh)
    cell_h = openmc.Cell(fill=h, region=+surf_wh)
    uni_wh = openmc.Universe(cells=(cell_w, cell_h))
    cell_f = openmc.Cell(name='f',fill=f, region=-surf_f)
    cell_wh = openmc.Cell(name='wh',fill=uni_wh, region=+surf_f & -surf_w)
    geometry = openmc.Geometry([cell_f, cell_wh])
    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50
    return openmc.Model(geometry, materials, settings)

@pytest.mark.parametrize("object_name, attribute", [
    ('wh', 'translation'),
    ('wh', 'rotation'),
    ('f', 'temperature'),
    ('f', 'refuel'),
    ])
def test_attributes(model, object_name, attribute):
    """
    Test classes attributes are set correctly
    """
    op = CoupledOperator(model, CHAIN_PATH)

    bracket = [-1,1]
    bracket_limit = [-10,10]
    axis = 2
    mat_vector = {'U235':0.1, 'U238':0.9}
    if attribute in ('translation','rotation'):
        bw = BatchwiseCellGeometrical(object_name, attribute, op, model, axis, bracket, bracket_limit)
        assert bw.cell_materials == [cell.fill for cell in  \
                model.geometry.get_cells_by_name(object_name)[0].fill.cells.values() \
                if cell.fill.depletable]
        assert bw.attrib_name == attribute
        assert bw.axis == axis

    elif attribute == 'temperature':
        bw = BatchwiseCellTemperature(object_name, op, model, bracket, bracket_limit)
    elif attribute == 'refuel':
        bw = BatchwiseMaterialRefuel(object_name, op, model, mat_vector, bracket, bracket_limit)
        assert bw.mat_vector == mat_vector

    assert bw.bracket == bracket
    assert bw.bracket_limit == bracket_limit
    assert bw.burn_mats == op.burnable_mats
    assert bw.local_mats == op.local_mats

@pytest.mark.parametrize("attribute", [
    ('translation'),
    ('rotation')
    ])
def test_cell_methods(run_in_tmpdir, model, attribute):
    """
    Test cell base class internal method
    """
    bracket = [-1,1]
    bracket_limit = [-10,10]
    axis = 2
    op = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, timestep_units = 'd')
    integrator.add_batchwise('wh', attribute, axis=axis, bracket=bracket,
                             bracket_limit=bracket_limit)

    model.export_to_xml()
    openmc.lib.init()
    integrator.batchwise._set_cell_attrib(0)
    assert integrator.batchwise._get_cell_attrib() == 0
    vol = integrator.batchwise._calculate_volumes()
    for cell in integrator.batchwise.cell_materials:
        assert vol[str(cell.id)] == pytest.approx([mat.volume for mat in model.materials if mat.id == cell.id])
    openmc.lib.finalize()

def test_update_materials_methods(run_in_tmpdir):
    """
    It's the same methods for all classe so one class is enough
    """
    bracket = [-1,1]
    bracket_limit = [-10,10]
    axis = 2
    op = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, timestep_units = 'd')
    integrator.add_batchwise('wh', 'translation', axis=axis, bracket=bracket,
                             bracket_limit=bracket_limit)

    model.export_to_xml()
    openmc.lib.init()
    x = op.number.number
    #Double  number of atoms of U238 in fuel
    mat_index = op.number.index_mat['1']
    nuc_index = op.number.index_nuc['U238']
    vol = op.number.get_mat_volume('1')
    op.number.number[mat_index][nuc_idx] *= 2
    integrator.batchwise._update_volumes(op.number.number)
    
    vol = integrator.batchwise._calculate_volumes()
