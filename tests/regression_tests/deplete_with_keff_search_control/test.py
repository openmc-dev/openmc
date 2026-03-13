""" Tests for KeffSearchControl class """

from hmac import new
from pathlib import Path
import shutil
import sys

import pytest
import numpy as np

import openmc
import openmc.lib
from openmc.deplete import CoupledOperator

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


def translate_cell(position):
    cell_trans = [cell for cell in openmc.lib.cells.values() if cell.name == 'trans_cell'][0]
    openmc.lib.cells[cell_trans.id].translation = [0, 0, position]


def rotate_cell(angle):
    cell_rot = [cell for cell in openmc.lib.cells.values() if cell.name == 'rot_cell'][0]
    openmc.lib.cells[cell_rot.id].rotation = [0, 0, angle]


def set_u235_density(u235_density):
    fuel = [material for material in openmc.lib.materials.values()
            if material.name == 'f'][0]
    nuclides = openmc.lib.materials[fuel.id].nuclides
    densities = openmc.lib.materials[fuel.id].densities
    nuc_idx = nuclides.index('U235')
    densities[nuc_idx] = u235_density
    openmc.lib.materials[fuel.id].set_densities(nuclides, densities)


@pytest.mark.parametrize("function, x0, x1, bracket, ref_result", [
    (translate_cell, -11, -5, (-15, 0), 'depletion_with_translation'),
    (rotate_cell, -80, -50, (-90, 0), 'depletion_with_rotation'),
    (set_u235_density, 2e-4, 1e-3, (1e-4, 2e-3), 'depletion_with_refuel')
])
def test_keff_search_control(run_in_tmpdir, model, function, x0, x1, bracket, ref_result):
    chain_file = Path(__file__).parents[2] / 'chain_simple.xml'
    model.settings.verbosity = 1
    op = CoupledOperator(model, chain_file)

    integrator = openmc.deplete.PredictorIntegrator(
        op, [1], 174., timestep_units = 'd')
    integrator.add_keff_search_control(
        function=function,
        x0=x0,
        x1=x1,
        bracket=bracket,
        output=True,
        k_tol=0.1,
        sigma_final=5e-2)

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
    assert res_test[0].keff_search_root == pytest.approx(res_ref[0].keff_search_root, rel=2)
