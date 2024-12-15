""" Full system test suite. """

from math import floor
import shutil
from pathlib import Path
from collections import defaultdict

from difflib import unified_diff
import numpy as np
import pytest
import openmc
from openmc.data import JOULE_PER_EV
import openmc.deplete

from tests.regression_tests import config, assert_atoms_equal
from .example_geometry import generate_problem


@pytest.fixture(scope="module")
def problem():
    n_rings = 2
    n_wedges = 4

    # Load geometry from example
    return generate_problem(n_rings, n_wedges)


@pytest.mark.parametrize("multiproc", [True, False])
def test_full(run_in_tmpdir, problem, multiproc):
    """Full system test suite.

    Runs an entire OpenMC simulation with depletion coupling and verifies
    that the outputs match a reference file.  Sensitive to changes in
    OpenMC.

    This test runs a complete OpenMC simulation and tests the outputs.
    It will take a while.

    """

    geometry, lower_left, upper_right = problem

    # OpenMC-specific settings
    settings = openmc.Settings()
    settings.particles = 100
    settings.batches = 10
    settings.inactive = 0
    space = openmc.stats.Box(lower_left, upper_right)
    settings.source = openmc.IndependentSource(space=space)
    settings.seed = 1
    settings.verbosity = 1

    model = openmc.Model(geometry=geometry, settings=settings)

    # Create operator
    chain_file = Path(__file__).parents[2] / "chain_simple.xml"
    op = openmc.deplete.CoupledOperator(model, chain_file)
    op.round_number = True

    # Power and timesteps
    dt1 = 15.0 * 24 * 60 * 60  # 15 days
    dt2 = 1.5 * 30 * 24 * 60 * 60  # 1.5 months
    N = floor(dt2 / dt1)
    dt = np.full(N, dt1)
    power = 2.337e15 * 4 * JOULE_PER_EV * 1e6  # MeV/second cm from CASMO

    # Perform simulation using the predictor algorithm
    openmc.deplete.pool.USE_MULTIPROCESSING = multiproc
    openmc.deplete.PredictorIntegrator(op, dt, power).integrate()

    # Get path to test and reference results
    path_test = op.output_dir / "depletion_results.h5"
    path_reference = Path(__file__).with_name("test_reference.h5")

    # If updating results, do so and return
    if config["update"]:
        shutil.copyfile(str(path_test), str(path_reference))
        return

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)
    res_ref = openmc.deplete.Results(path_reference)

    # Assert same mats
    for mat in res_ref[0].index_mat:
        assert mat in res_test[0].index_mat, "Material {} not in new results.".format(
            mat
        )
    for nuc in res_ref[0].index_nuc:
        assert nuc in res_test[0].index_nuc, "Nuclide {} not in new results.".format(
            nuc
        )

    for mat in res_test[0].index_mat:
        assert mat in res_ref[0].index_mat, "Material {} not in old results.".format(
            mat
        )
    for nuc in res_test[0].index_nuc:
        assert nuc in res_ref[0].index_nuc, "Nuclide {} not in old results.".format(nuc)

    assert_atoms_equal(res_ref, res_test, tol=1e-6)

    # Compare statepoint files with depletion results

    t_test, k_test = res_test.get_keff()
    t_ref, k_ref = res_ref.get_keff()
    k_state = np.empty_like(k_ref)

    n_tallies = np.empty(N + 1, dtype=int)

    # Get statepoint files for all BOS points and EOL
    runtimes = defaultdict(list)
    for n in range(N + 1):
        statepoint = openmc.StatePoint(f"openmc_simulation_n{n}.h5")
        for measure, time in statepoint.runtime.items():
            runtimes[measure].append(time)
        k_n = statepoint.keff
        k_state[n] = [k_n.nominal_value, k_n.std_dev]
        n_tallies[n] = len(statepoint.tallies)
    # Look for exact match pulling from statepoint and depletion_results
    assert np.all(k_state == k_test)
    assert np.allclose(k_test, k_ref)

    # Check that no additional tallies are loaded from the files
    assert np.all(n_tallies == 0)

    # Convert values in runtimes to arrays
    runtimes = {k: np.array(v) for k, v in runtimes.items()}

    # Check that runtimes are qualitatively correct
    assert runtimes["reading cross sections"][0] != 0
    assert runtimes["total initialization"][0] != 0
    assert np.all(runtimes["reading cross sections"][1:] == 0)
    assert np.all(runtimes["total initialization"][1:] == 0)
    assert np.all(runtimes["inactive batches"] == 0)
    del runtimes["reading cross sections"]
    del runtimes["total initialization"]
    del runtimes["inactive batches"]
    for measure, times in runtimes.items():
        assert np.all(times != 0)


def test_depletion_results_to_material(run_in_tmpdir, problem):
    """Checks openmc.Materials objects can be created from depletion results"""
    # Load the reference/test results
    path_reference = Path(__file__).with_name("test_reference.h5")
    res_ref = openmc.deplete.Results(path_reference)

    # Firstly need to export materials.xml file for the initial simulation state
    geometry, lower_left, upper_right = problem
    materials = openmc.Materials()
    for mat in geometry.root_universe.get_all_materials().values():
        materials.append(mat)
    materials.export_to_xml()

    # Export last step of depletion to its own openmc.Materials object,
    # using only nuclides available in the current nuclear data library
    last_step_materials = res_ref.export_to_materials(-1)

    # Export final depletion  step materials to XML
    output_xml_file = "last_step_materials.xml"
    last_step_materials.export_to_xml(path=output_xml_file)
    with open(output_xml_file, "r") as result_file:
        result_file_lines = result_file.readlines()

    # If updating results, do so and return. We write out the last-step
    # depleted materials as an XML, and save the list of lines to diff.
    reference_file = Path(__file__).with_name("last_step_reference_materials.xml")
    if config["update"]:
        with open(reference_file, "w") as ref_file:
            ref_file.writelines(result_file_lines)
        return

    # Check text of final depletion point material XML matches reference
    with open(reference_file) as ref_file:
        reference_lines = ref_file.readlines()
    diff_vs_expected = unified_diff(reference_lines, result_file_lines)

    # Check all lines match, printing errors along the way
    success = True
    for line in diff_vs_expected:
        success = False
        print(line.rstrip())
    assert success
