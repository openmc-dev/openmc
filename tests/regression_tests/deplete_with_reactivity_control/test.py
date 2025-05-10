""" Tests for ReactivityController class """

from pathlib import Path
import shutil
import sys

import pytest
import numpy as np

import openmc
import openmc.lib
from openmc.deplete import CoupledOperator
from openmc.deplete import (
    CellReactivityController,
    MaterialReactivityController
)

from tests.regression_tests import config

@pytest.fixture
def model():
    f = openmc.Material(name='f')
    f.set_density('g/cm3', 10.29769)
    f.add_element('U', 1., enrichment=2.4)
    f.add_element('O', 2.)

    h = openmc.Material(name='h')
    h.set_density('g/cm3', 0.001598)
    h.add_element('He', 2.4044e-4)

    w = openmc.Material(name='w')
    w.set_density('g/cm3', 0.740582)
    w.add_element('H', 2)
    w.add_element('O', 1)

    # Define overall material
    materials = openmc.Materials([f, h, w])

    # Define surfaces
    radii = [0.5, 0.8, 1]
    height = 80
    surf_in = openmc.ZCylinder(r=radii[0])
    surf_mid = openmc.ZCylinder(r=radii[1])
    surf_out = openmc.ZCylinder(r=radii[2], boundary_type='reflective')
    surf_top = openmc.ZPlane(z0=height/2, boundary_type='vacuum')
    surf_bot = openmc.ZPlane(z0=-height/2, boundary_type='vacuum')

    surf_trans = openmc.ZPlane(z0=0)
    surf_rot1 = openmc.XPlane(x0=0)
    surf_rot2 = openmc.YPlane(y0=0)

    # Define cells
    cell_f = openmc.Cell(name='fuel_cell', fill=f,
            region=-surf_in  & -surf_top & +surf_bot)
    cell_g = openmc.Cell(fill=h,
            region = +surf_in & -surf_mid & -surf_top & +surf_bot & +surf_rot2)

    # Define unbounded cells for rotation universe
    cell_w = openmc.Cell(fill=w, region = -surf_rot1)
    cell_h = openmc.Cell(fill=h, region = +surf_rot1)
    universe_rot = openmc.Universe(cells=(cell_w, cell_h))
    cell_rot = openmc.Cell(name="rot_cell", fill=universe_rot,
            region = +surf_in & -surf_mid & -surf_top & +surf_bot & -surf_rot2)

    # Define unbounded cells for translation universe
    cell_w = openmc.Cell(fill=w, region=+surf_in  & -surf_trans )
    cell_h = openmc.Cell(fill=h,  region=+surf_in & +surf_trans)
    universe_trans = openmc.Universe(cells=(cell_w, cell_h))
    cell_trans = openmc.Cell(name="trans_cell", fill=universe_trans,
                          region=+surf_mid & -surf_out & -surf_top & +surf_bot)

    # Define overall geometry
    geometry = openmc.Geometry([cell_f, cell_g, cell_rot, cell_trans])

    # Set material volume for depletion fuel.
    f.volume = np.pi * radii[0]**2 * height

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)

@pytest.fixture
def mix_mat():
    mix_mat = openmc.Material()
    mix_mat.add_element("U", 1, percent_type="ao", enrichment=90)
    mix_mat.add_element("O", 2)
    mix_mat.set_density("g/cc", 10.4)
    return mix_mat

@pytest.mark.parametrize("obj, attribute, bracket, bracket_limit, axis, ref_result", [
    ('trans_cell', 'translation', [-5,5], [-40,40], 2, 'depletion_with_translation'),
    ('rot_cell', 'rotation', [-5,5], [-90,90], 2, 'depletion_with_rotation'),
    ('f', None, [-1,2], [-100,100], None, 'depletion_with_refuel')
    ])
def test_reactivity_control(run_in_tmpdir, model, obj, attribute, bracket,
                    bracket_limit, axis, ref_result, mix_mat):

    chain_file = Path(__file__).parents[2] / 'chain_simple.xml'
    op = CoupledOperator(model, chain_file)

    integrator = openmc.deplete.PredictorIntegrator(
        op, [1], 174., timestep_units = 'd')

    if attribute:
        integrator.reactivity_control = CellReactivityController(
                cell=obj,
                operator=op,
                attribute=attribute,
                bracket=bracket,
                bracket_limit=bracket_limit,
                axis=axis
        )
    else:
        integrator.reactivity_control = MaterialReactivityController(
                material=obj,
                operator=op,
                material_to_mix=mix_mat,
                bracket=bracket,
                bracket_limit=bracket_limit,
        )
    integrator.integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'
    path_reference = Path(__file__).with_name(f'ref_{ref_result}.h5')

    # If updating results, do so and return
    if config['update']:
        shutil.copyfile(str(path_test), str(path_reference))
        return

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)
    res_ref = openmc.deplete.Results(path_reference)

    # Use high tolerance here
    assert res_test[0].reac_cont == pytest.approx(res_ref[0].reac_cont, rel=2)
