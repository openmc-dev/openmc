""" Transport-free depletion test suite """

from pathlib import Path

import numpy as np
import pytest
import openmc
import openmc.deplete
from openmc.deplete import IndependentOperator, MicroXS


@pytest.fixture(scope="module")
def fuel():
    fuel = openmc.Material(name="uo2")
    fuel.add_element("U", 1, percent_type="ao", enrichment=4.25)
    fuel.add_element("O", 2)
    fuel.set_density("g/cc", 10.4)
    fuel.depletable = True

    fuel.volume = np.pi * 0.42 ** 2

    return fuel

@pytest.fixture(scope="module")
def micro_xs():
    micro_xs_file = Path(__file__).parents[2] / 'micro_xs_simple.csv'
    return MicroXS.from_csv(micro_xs_file)


@pytest.fixture(scope="module")
def chain_file():
    return Path(__file__).parents[2] / 'chain_simple.xml'


@pytest.mark.parametrize("multiproc, from_nuclides, normalization_mode, power, flux", [
    (True, True, 'source-rate', None, 1164719970082145.0),
    (False, True, 'source-rate', None, 1164719970082145.0),
    (True, True, 'fission-q', 174, None),
    (False, True, 'fission-q', 174, None),
    (True, False, 'source-rate', None, 1164719970082145.0),
    (False, False, 'source-rate', None, 1164719970082145.0),
    (True, False, 'fission-q', 174, None),
    (False, False, 'fission-q', 174, None)])
def test_against_self(run_in_tmpdir,
                      fuel,
                      micro_xs,
                      chain_file,
                      multiproc,
                      from_nuclides,
                      normalization_mode,
                      power,
                      flux):
    """Transport free system test suite.

    Runs an OpenMC transport-free depletion calculation and verifies
    that the outputs match a reference file.

    """
    # Create operator
    op = _create_operator(from_nuclides,
                          fuel,
                          micro_xs,
                          chain_file,
                          normalization_mode)

    # Power and timesteps
    dt = [360]  # single step

    # Perform simulation using the predictor algorithm
    openmc.deplete.pool.USE_MULTIPROCESSING = multiproc
    openmc.deplete.PredictorIntegrator(op,
                                       dt,
                                       power=power,
                                       source_rates=flux,
                                       timestep_units='s').integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'
    if flux is not None:
        ref_path = 'test_reference_source_rate.h5'
    else:
        ref_path = 'test_reference_fission_q.h5'
    path_reference = Path(__file__).with_name(ref_path)

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)
    res_ref = openmc.deplete.Results(path_reference)

    # Assert same mats
    _assert_same_mats(res_test, res_ref)

    tol = 1.0e-14
    _assert_atoms_equal(res_test, res_ref, tol)
    _assert_reaction_rates_equal(res_test, res_ref, tol)


@pytest.mark.parametrize("multiproc, dt, time_units, time_type, atom_tol, rx_tol ", [
    (True, 360, 's', 'minutes', 2.0e-3, 3.0e-2),
    (False, 360, 's', 'minutes', 2.0e-3, 3.0e-2),
    (True, 4, 'h', 'hours', 2.0e-3, 6.0e-2),
    (False, 4, 'h', 'hours', 2.0e-3, 6.0e-2),
    (True, 5, 'd', 'days', 2.0e-3, 5.0e-2),
    (False, 5, 'd', 'days', 2.0e-3, 5.0e-2),
    (True, 100, 'd', 'months', 4.0e-3, 9.0e-2),
    (False, 100, 'd', 'months', 4.0e-3, 9.0e-2)])
def test_against_coupled(run_in_tmpdir,
                         fuel,
                         micro_xs,
                         chain_file,
                         multiproc,
                         dt,
                         time_units,
                         time_type,
                         atom_tol,
                         rx_tol):
    # Create operator
    op = _create_operator(False, fuel, micro_xs, chain_file, 'fission-q')

    # Power and timesteps
    dt = [dt]  # single step

    # Perform simulation using the predictor algorithm
    openmc.deplete.pool.USE_MULTIPROCESSING = multiproc
    openmc.deplete.PredictorIntegrator(
        op, dt, power=174, timestep_units=time_units).integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'

    ref_path = f'test_reference_coupled_{time_type}.h5'
    path_reference = Path(__file__).with_name(ref_path)

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)
    res_ref = openmc.deplete.Results(path_reference)

    # Assert same mats
    _assert_same_mats(res_test, res_ref)

    _assert_atoms_equal(res_test, res_ref, atom_tol)
    _assert_reaction_rates_equal(res_test, res_ref, rx_tol)


def _create_operator(from_nuclides,
                     fuel,
                     micro_xs,
                     chain_file,
                     normalization_mode):
    if from_nuclides:
        nuclides = {}
        for nuc, dens in fuel.get_nuclide_atom_densities().items():
            nuclides[nuc] = dens

        op = IndependentOperator.from_nuclides(fuel.volume,
                                               nuclides,
                                               micro_xs,
                                               chain_file,
                                               normalization_mode=normalization_mode)

    else:
        op = IndependentOperator(openmc.Materials([fuel]),
                                 micro_xs,
                                 chain_file,
                                 normalization_mode=normalization_mode)

    return op


def _assert_same_mats(res_ref, res_test):
    for mat in res_ref[0].index_mat:
        assert mat in res_test[0].index_mat, \
            f"Material {mat} not in new results."
    for nuc in res_ref[0].index_nuc:
        assert nuc in res_test[0].index_nuc, \
            f"Nuclide {nuc} not in new results."

    for mat in res_test[0].index_mat:
        assert mat in res_ref[0].index_mat, \
            f"Material {mat} not in old results."
    for nuc in res_test[0].index_nuc:
        assert nuc in res_ref[0].index_nuc, \
            f"Nuclide {nuc} not in old results."


def _assert_atoms_equal(res_ref, res_test, tol):
    for mat in res_test[0].index_mat:
        for nuc in res_test[0].index_nuc:
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


def _assert_reaction_rates_equal(res_ref, res_test, tol):
    for reactions in res_test[0].rates:
        for mat in reactions.index_mat:
            for nuc in reactions.index_nuc:
                for rx in reactions.index_rx:
                    y_test = res_test.get_reaction_rate(mat, nuc, rx)[1] / \
                        res_test.get_atoms(mat, nuc)[1]
                    y_old = res_ref.get_reaction_rate(mat, nuc, rx)[1] / \
                        res_ref.get_atoms(mat, nuc)[1]

                    # Test each point
                    correct = True
                    for i, ref in enumerate(y_old):
                        if ref != y_test[i]:
                            if ref != 0.0:
                                correct = np.abs(y_test[i] - ref) / ref <= tol
                            else:
                                if y_test[i] != 0.0:
                                    correct = False

                    assert correct, "Discrepancy in mat {}, nuc {}, and rx {}\n{}\n{}".format(
                        mat, nuc, rx, y_old, y_test)
