""" Tests for Batchwise class """

from pathlib import Path

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
    f = openmc.Material(name="fuel")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)
    f.temperature = 293.15

    w = openmc.Material(name="water")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.temperature = 273.15
    w.depletable = True

    h = openmc.Material(name='helium')
    h.add_element('He', 1)
    h.set_density('g/cm3', 0.001598)

    radii = [0.42, 0.45]
    height = 0.5

    f.volume = np.pi * radii[0] ** 2 * height
    w.volume = np.pi * (radii[1]**2 - radii[0]**2) * height/2

    materials = openmc.Materials([f, w, h])

    surf_interface = openmc.ZPlane(z0=0)
    surf_top = openmc.ZPlane(z0=height/2)
    surf_bot = openmc.ZPlane(z0=-height/2)
    surf_in = openmc.Sphere(r=radii[0])
    surf_out = openmc.Sphere(r=radii[1], boundary_type='vacuum')

    cell_water = openmc.Cell(fill=w, region=-surf_interface)
    cell_helium = openmc.Cell(fill=h, region=+surf_interface)
    universe = openmc.Universe(cells=(cell_water, cell_helium))
    cell_fuel = openmc.Cell(name='fuel_cell', fill=f,
                            region=-surf_in & -surf_top & +surf_bot)
    cell_universe = openmc.Cell(name='universe_cell',fill=universe,
                            region=+surf_in & -surf_out & -surf_top & +surf_bot)
    geometry = openmc.Geometry([cell_fuel, cell_universe])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)

@pytest.mark.parametrize("case_name, obj, attribute, bracket, limit, axis, vec", [
    ('cell translation','universe_cell', 'translation', [-1,1], [-10,10], 2, None),
    ('cell rotation', 'universe_cell', 'rotation', [-1,1], [-10,10], 2, None),
    ('single cell temperature', 'fuel_cell', 'temperature', [-1,1], [-10,10], 2, None ),
    ('multi cell', 'universe_cell', 'temperature',  [-1,1], [-10,10], 2, None ),
    ('material refuel', 'fuel', 'refuel', [-1,1], [-10,10], None, {'U235':0.1, 'U238':0.9}),
    ('invalid_1', 'universe_cell', 'refuel', [-1,1], [-10,10], None, {'U235':0.1, 'U238':0.9}),
    ('invalid_2', 'fuel', 'temperature',  [-1,1], [-10,10], 2, None ),
    ('invalid_3', 'fuel', 'translation', [-1,1], [-10,10], 2, None),
    ])
def test_attributes(case_name, model, obj, attribute, bracket, limit, axis, vec):
    """
    Test classes attributes are set correctly
    """
    op = CoupledOperator(model, CHAIN_PATH)

    int = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, timestep_units = 'd')

    kwargs = {'bracket': bracket, 'bracket_limit':limit, 'tol': 0.1,}

    if vec is not None:
        kwargs['mat_vector']=vec

    if axis is not None:
        kwargs['axis'] = axis

    if case_name == "invalid_1":
        with pytest.raises(ValueError) as e:
            int.add_batchwise(obj, attribute, **kwargs)
        assert str(e.value) == 'Unable to set "Material name" to "universe_cell" '\
                               'since it is not in "[\'fuel\', \'water\']"'
    elif case_name == "invalid_2":
        with pytest.raises(ValueError) as e:
            int.add_batchwise(obj, attribute, **kwargs)
        assert str(e.value) == 'Unable to set "Cell name exists" to "fuel" since '\
                   'it is not in "[\'fuel_cell\', \'universe_cell\', \'\', \'\']"'

    elif case_name == "invalid_3":
        with pytest.raises(ValueError) as e:
            int.add_batchwise(obj, attribute, **kwargs)
        assert str(e.value) == 'Unable to set "Cell name exists" to "fuel" since '\
                   'it is not in "[\'fuel_cell\', \'universe_cell\', \'\', \'\']"'
    else:
        int.add_batchwise(obj, attribute, **kwargs)
        if attribute in ('translation','rotation'):
            assert int.batchwise.cell_materials == [cell.fill for cell in  \
                    model.geometry.get_cells_by_name(obj)[0].fill.cells.values() \
                    if cell.fill.depletable]
            assert int.batchwise.axis == axis

        elif attribute == 'refuel':
            assert int.batchwise.mat_vector == vec

        assert int.batchwise.attrib_name == attribute
        assert int.batchwise.bracket == bracket
        assert int.batchwise.bracket_limit == limit
        assert int.batchwise.burn_mats == op.burnable_mats
        assert int.batchwise.local_mats == op.local_mats

@pytest.mark.parametrize("obj, attribute, value_to_set", [
    ('universe_cell', 'translation', 0),
    ('universe_cell', 'rotation',    0)
    ])
def test_cell_methods(run_in_tmpdir, model, obj, attribute, value_to_set):
    """
    Test cell base class internal method
    """
    kwargs = {'bracket':[-1,1], 'bracket_limit':[-10,10], 'axis':2, 'tol':0.1}

    op = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, timestep_units = 'd')
    integrator.add_batchwise(obj, attribute, **kwargs)

    model.export_to_xml()
    openmc.lib.init()
    integrator.batchwise._set_cell_attrib(value_to_set)
    assert integrator.batchwise._get_cell_attrib() == value_to_set

    vol = integrator.batchwise._calculate_volumes()
    for cell in integrator.batchwise.cell_materials:
        assert vol[str(cell.id)] == pytest.approx([
                                    mat.volume for mat in model.materials \
                                    if mat.id == cell.id][0], rel=tolerance)

    openmc.lib.finalize()

@pytest.mark.parametrize("nuclide, atoms_to_add", [
    ('U238', 1.0e22),
    ('Xe135', 1.0e21)
    ])
def test_internal_methods(run_in_tmpdir, model, nuclide, atoms_to_add):
    """
    Method to update volume in AtomNumber after depletion step.
    Method inheritated by all derived classes so one check is enough.
    """

    kwargs = {'bracket':[-1,1], 'bracket_limit':[-10,10], 'mat_vector':{}}

    op = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(op, [1,1], 0.0,
                                                    timestep_units = 'd')
    integrator.add_batchwise('fuel', 'refuel', **kwargs)

    model.export_to_xml()
    openmc.lib.init()

    #Increase  number of atoms of U238 in fuel by fix amount and check the
    # volume increase at constant-density
    #extract fuel material from model materials
    mat = integrator.batchwise.material
    mat_index = op.number.index_mat[str(mat.id)]
    nuc_index = op.number.index_nuc[nuclide]
    vol = op.number.get_mat_volume(str(mat.id))
    op.number.number[mat_index][nuc_index] += atoms_to_add
    integrator.batchwise._update_volumes()

    vol_to_compare = vol + (atoms_to_add * openmc.data.atomic_mass(nuclide) /\
                            openmc.data.AVOGADRO / mat.density)

    assert op.number.get_mat_volume(str(mat.id)) == pytest.approx(vol_to_compare)

    x = [i[:op.number.n_nuc_burn] for i in op.number.number]
    integrator.batchwise._update_materials(x)
    nuc_index_lib = openmc.lib.materials[mat.id].nuclides.index(nuclide)
    dens_to_compare = 1.0e-24 * op.number.get_atom_density(str(mat.id), nuclide)

    assert openmc.lib.materials[mat.id].densities[nuc_index_lib] == pytest.approx(dens_to_compare)

    volumes = {str(mat.id): vol + 1}
    new_x = integrator.batchwise._update_x_and_set_volumes(x, volumes)
    dens_to_compare = 1.0e24 *  volumes[str(mat.id)] *\
                      openmc.lib.materials[mat.id].densities[nuc_index_lib]
    assert new_x[mat_index][nuc_index] == pytest.approx(dens_to_compare)

    # assert volume in AtomNumber is set correctly
    assert op.number.get_mat_volume(str(mat.id)) == volumes[str(mat.id)]

    openmc.lib.finalize()
