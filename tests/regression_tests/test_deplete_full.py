""" Full system test suite. """

from math import floor
import shutil
from pathlib import Path

import numpy as np
import openmc
from openmc.data import JOULE_PER_EV
import openmc.deplete
from openmc.deplete import results
from openmc.deplete import utilities

from tests.regression_tests import config
from .example_geometry import generate_problem


def test_full(run_in_tmpdir):
    """Full system test suite.

    Runs an entire OpenMC simulation with depletion coupling and verifies
    that the outputs match a reference file.  Sensitive to changes in
    OpenMC.

    This test runs a complete OpenMC simulation and tests the outputs.
    It will take a while.
    """

    n_rings = 2
    n_wedges = 4

    # Load geometry from example
    geometry, lower_left, upper_right = generate_problem(n_rings, n_wedges)

    # Create dt vector for 3 steps with 15 day timesteps
    dt1 = 15.*24*60*60  # 15 days
    dt2 = 1.5*30*24*60*60  # 1.5 months
    N = floor(dt2/dt1)
    dt = np.full(N, dt1)

    # Depletion settings
    settings = openmc.deplete.OpenMCSettings()
    settings.chain_file = str(Path(__file__).parents[2] / 'chains' /
                              'chain_simple.xml')
    settings.power = 2.337e15*4*JOULE_PER_EV*1e6  # MeV/second cm from CASMO
    settings.dt_vec = dt
    settings.round_number = True

    # Add OpenMC-specific settings
    settings.particles = 100
    settings.batches = 100
    settings.inactive = 40
    space = openmc.stats.Box(lower_left, upper_right)
    settings.source = openmc.Source(space=space)
    settings.seed = 1
    settings.verbosity = 3

    op = openmc.deplete.OpenMCOperator(geometry, settings)

    # Perform simulation using the predictor algorithm
    openmc.deplete.integrator.predictor(op)

    # Get path to test and reference results
    path_test = settings.output_dir / 'depletion_results.h5'
    path_reference = Path(__file__).with_name('test_reference.h5')

    # If updating results, do so and return
    if config['update']:
        shutil.copyfile(str(path_test), str(path_reference))
        return

    # Load the reference/test results
    res_test = results.read_results(path_test)
    res_ref = results.read_results(path_reference)

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
            _, y_test = utilities.evaluate_single_nuclide(res_test, mat, nuc)
            _, y_old = utilities.evaluate_single_nuclide(res_ref, mat, nuc)

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
