""" Tests for ReactivityController class """

from pathlib import Path

import pytest
import numpy as np

import openmc
import openmc.lib
from openmc.deplete import CoupledOperator
from openmc.deplete import (
    CellReactivityController,
    MaterialReactivityController
)

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
    w.temperature = 293.15
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

@pytest.fixture
def operator(model):
    return CoupledOperator(model, CHAIN_PATH)

@pytest.fixture
def integrator(operator):
    return openmc.deplete.PredictorIntegrator(
            operator, [1,1], 0.0, timestep_units = 'd')

@pytest.fixture
def mix_mat():
    mix_mat = openmc.Material()
    mix_mat.add_element("U", 1, percent_type="ao", enrichment=90)
    mix_mat.add_element("O", 2)
    mix_mat.set_density("g/cc", 10.4)
    return mix_mat

@pytest.mark.parametrize("case_name, obj, attribute, bracket, limit, axis", [
    ('cell translation','universe_cell', 'translation', [-1,1], [-10,10], 2),
    ('cell rotation', 'universe_cell', 'rotation', [-1,1], [-10,10], 2),
    ('material mix', 'fuel', None, [0,5], [-10,10], None),
    ])
def test_attributes(case_name, model, operator, integrator, obj, attribute,
                    bracket, limit, axis, mix_mat):
    """
    Test classes attributes are set correctly
    """
    if attribute:
        integrator.reactivity_control = CellReactivityController(
                cell=obj,
                operator=operator,
                attribute=attribute,
                bracket=bracket,
                bracket_limit=limit,
                axis=axis
        )
        assert integrator.reactivity_control.depletable_cells == [cell for cell in  \
                model.geometry.get_cells_by_name(obj)[0].fill.cells.values() \
                if cell.fill.depletable]
        assert integrator.reactivity_control.axis == axis

    else:
        integrator.reactivity_control = MaterialReactivityController(
                material=obj,
                operator=operator,
                material_to_mix=mix_mat,
                bracket=bracket,
                bracket_limit=limit
        )
        assert integrator.reactivity_control.material_to_mix == mix_mat

    assert integrator.reactivity_control.bracket == bracket
    assert integrator.reactivity_control.bracket_limit == limit
    assert integrator.reactivity_control.burn_mats == operator.burnable_mats
    assert integrator.reactivity_control.local_mats == operator.local_mats

@pytest.mark.parametrize("obj, attribute, value_to_set", [
    ('universe_cell', 'translation', 0),
    ('universe_cell', 'rotation',    0)
    ])
def test_cell_methods(run_in_tmpdir, model, operator, integrator, obj, attribute,
                      value_to_set):
    """
    Test cell base class internal method
    """
    integrator.reactivity_control = CellReactivityController(
                    cell=obj,
                    operator=operator,
                    attribute=attribute,
                    bracket=[-1,1],
                    bracket_limit=[-10,10],
                    axis=2
            )

    model.export_to_xml()
    openmc.lib.init()
    integrator.reactivity_control._set_lib_cell()
    integrator.reactivity_control._set_cell_attrib(value_to_set)
    assert integrator.reactivity_control._get_cell_attrib() == value_to_set

    vol = integrator.reactivity_control._adjust_volumes()
    integrator.reactivity_control._update_x_and_set_volumes(operator.number.number, vol)

    for cell in integrator.reactivity_control.depletable_cells:
        mat_id = str(cell.fill.id)
        index_mat = operator.number.index_mat[mat_id]
        assert vol[mat_id] == operator.number.volume[index_mat]

    openmc.lib.finalize()

@pytest.mark.parametrize("units, value, dilute", [
    ('cc', 1.0, True),
    ('cm3', 1.0, False),
    ('grams', 1.0, True),
    ('grams', 1.0, False),
    ('atoms', 1.0, True),
    ('atoms', 1.0, False),
    ])
def test_internal_methods(run_in_tmpdir, model, operator, integrator, mix_mat,
                          units, value, dilute):
    """
    Method to update volume in AtomNumber after depletion step.
    Method inheritated by all derived classes so one check is enough.
    """

    integrator.reactivity_control = MaterialReactivityController(
            material='fuel',
            operator=operator,
            material_to_mix=mix_mat,
            bracket=[-1,1],
            bracket_limit=[-10,10],
            units=units,
            dilute=dilute
    )

    model.export_to_xml()
    openmc.lib.init()

    # check material volume are adjusted correctly
    mat = integrator.reactivity_control.material
    mat_index = operator.number.index_mat[str(mat.id)]
    nuc_index = operator.number.index_nuc['U235']
    vol_before_mix = operator.number.get_mat_volume(str(mat.id))

    # set mix_mat volume
    integrator.reactivity_control._set_mix_material_volume(value)
    mix_mat_vol = mix_mat.volume
    #extract nuclide in atoms
    mix_nuc_atoms = mix_mat.get_nuclide_atoms()['U235']
    operator.number.number[mat_index][nuc_index] += mix_nuc_atoms
    #update volume
    vol_after_mix = integrator.reactivity_control._adjust_volumes(mix_mat_vol)

    if dilute:
        assert vol_before_mix + mix_mat_vol == vol_after_mix[str(mat.id)]
    else:
        assert vol_before_mix == vol_after_mix[str(mat.id)]

    #check nuclide densities get assigned correctly in memory
    x = [i[:operator.number.n_nuc_burn] for i in operator.number.number]
    integrator.reactivity_control._update_materials(x)
    nuc_index_lib = openmc.lib.materials[mat.id].nuclides.index('U235')
    dens_to_compare = 1.0e-24 * operator.number.get_atom_density(str(mat.id), 'U235')

    assert openmc.lib.materials[mat.id].densities[nuc_index_lib] == pytest.approx(dens_to_compare)

    # check volumes get assigned correctly in AtomNumber
    new_x = integrator.reactivity_control._update_x_and_set_volumes(x, vol_after_mix)
    dens_to_compare = 1.0e24 *  vol_after_mix[str(mat.id)] *\
                      openmc.lib.materials[mat.id].densities[nuc_index_lib]
    assert new_x[mat_index][nuc_index] == pytest.approx(dens_to_compare)

    # assert volume in AtomNumber is set correctly
    assert operator.number.get_mat_volume(str(mat.id)) == vol_after_mix[str(mat.id)]

    openmc.lib.finalize()
