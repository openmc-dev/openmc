""" TransferRates depletion test suite """

from pathlib import Path
import shutil
import sys

import numpy as np
import pytest
import openmc
import openmc.deplete
from openmc.deplete import CoupledOperator

from tests.regression_tests import config, assert_reaction_rates_equal, \
    assert_atoms_equal, assert_same_mats


@pytest.fixture
def model():
    openmc.reset_auto_ids()
    f = openmc.Material(name="f")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)

    w = openmc.Material(name="w")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.depletable = True

    radii = [0.42, 0.45]
    f.volume = np.pi * radii[0] ** 2
    w.volume = np.pi * (radii[1]**2 - radii[0]**2)

    materials = openmc.Materials([f, w])

    surf_f = openmc.Sphere(r=radii[0])
    surf_w = openmc.Sphere(r=radii[1], boundary_type='reflective')
    cell_f = openmc.Cell(fill=f, region=-surf_f)
    cell_w = openmc.Cell(fill=w, region=+surf_f & -surf_w)
    geometry = openmc.Geometry([cell_f, cell_w])

    settings = openmc.Settings()
    settings.particles = 100
    settings.inactive = 0
    settings.batches = 10

    return openmc.Model(geometry, materials, settings)

@pytest.mark.parametrize("rate, dest_mat, power, ref_result", [
    (1e-5, None, 0.0, 'no_depletion_only_removal'),
    (-1e-5, None, 0.0, 'no_depletion_only_feed'),
    (1e-5, None, 174.0, 'depletion_with_removal'),
    (-1e-5, None, 174.0, 'depletion_with_feed'),
    (-1e-5, 'w', 0.0, 'no_depletion_with_transfer'),
    (1e-5, 'w', 174.0, 'depletion_with_transfer'),
    ])
def test_transfer_rates(run_in_tmpdir, model, rate, dest_mat, power, ref_result):
    """Tests transfer_rates depletion class with transfer rates"""

    chain_file = Path(__file__).parents[2] / 'chain_simple.xml'

    transfer_elements = ['Xe']

    op = CoupledOperator(model, chain_file)
    op.round_number = True
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1], power, timestep_units = 'd')
    integrator.add_transfer_rate('f', transfer_elements, rate,
                                 destination_material=dest_mat)
    integrator.integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'
    path_reference = Path(__file__).with_name(f'ref_{ref_result}.h5')

    # If updating results, do so and return
    if config['update']:
        shutil.copyfile(str(path_test), str(path_reference))
        return

    # Load the reference/test results
    res_ref = openmc.deplete.Results(path_reference)
    res_test = openmc.deplete.Results(path_test)

    assert_same_mats(res_ref, res_test)
    assert_atoms_equal(res_ref, res_test, 1e-6)
    assert_reaction_rates_equal(res_ref, res_test)
