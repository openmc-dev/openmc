""" Transport-free depletion test suite """

from math import floor
import shutil
from pathlib import Path

from difflib import unified_diff
import numpy as np
import pytest
import openmc
from openmc.data import JOULE_PER_EV
import openmc.deplete
from openmc.deplete import FluxDepletionOperator


from tests.regression_tests import config


@pytest.fixture(scope="module")
def vol_nuc():
    fuel = openmc.Material(name="uo2")
    fuel.add_element("U", 1, percent_type="ao", enrichment=4.25)
    fuel.add_element("O", 2)
    fuel.set_density("g/cc", 10.4)

    fuel.volume = np.pi * 0.42 ** 2

    nuclides = {}
    for nuc, dens in fuel.get_nuclide_atom_densities().items():
        nuclides[nuc] = dens * 1e24

    # Load geometry from example
    return (fuel.volume, nuclides)


@pytest.mark.parametrize("multiproc", [True, False])
def test_no_transport(run_in_tmpdir, vol_nuc, multiproc):
    """Transport free system test suite.

    Runs an OpenMC transport-free depletion calculation and verifies
    that the outputs match a reference file.

    """

    # Create operator
    micro_xs_file = Path(__file__).parents[2] / 'micro_xs_simple.csv'
    micro_xs = FluxDepletionOperator.create_micro_xs_from_csv(micro_xs_file)
    chain_file = Path(__file__).parents[2] / 'chain_simple.xml'
    flux = 1164719970082145.0 # flux from pincell example
    op = FluxDepletionOperator(vol_nuc[0], vol_nuc[1], micro_xs, flux, chain_file)

    # Power and timesteps
    dt = [30] # single step
    power = 174 # W/cm

    # Perform simulation using the predictor algorithm
    openmc.deplete.pool.USE_MULTIPROCESSING = multiproc
    openmc.deplete.PredictorIntegrator(op, dt, power, timestep_units='d').integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'
    path_reference = Path(__file__).with_name('test_reference.h5')

    # If updating results, do so and return
    if config['update']:
        shutil.copyfile(str(path_test), str(path_reference))
        return

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)
    res_ref = openmc.deplete.Results(path_reference)

    # Assert same mats
    for mat in res_ref[0].mat_to_ind:
        assert mat in res_test[0].mat_to_ind, \
            "Material {} not in new results.".format(mat)
    for nuc in res_ref[0].nuc_to_ind:
        assert nuc in res_test[0].nuc_to_ind, \
            "Nuclide {} not in new results.".format(nuc)

    for mat in res_test[0].mat_to_ind:
        assert mat in res_ref[0].mat_to_ind, \
            "Material {} not in old results.".format(mat)
    for nuc in res_test[0].nuc_to_ind:
        assert nuc in res_ref[0].nuc_to_ind, \
            "Nuclide {} not in old results.".format(nuc)

    tol = 1.0e-6
    for mat in res_test[0].mat_to_ind:
        for nuc in res_test[0].nuc_to_ind:
            _, y_test = res_test.get_atoms(mat, nuc)
            _, y_old = res_ref.get_atoms(mat, nuc)

            # Test each point
            correct = True
            for i, ref in enumerate(y_old):
                if ref != y_test[i]:
                    if ref != 0.0:
                        correct = np.abs(y_test[i] - ref) / ref <= tol
                    else:
                        correct = False

            assert correct, "Discrepancy in mat {} and nuc {}\n{}\n{}".format(
                mat, nuc, y_old, y_test)

