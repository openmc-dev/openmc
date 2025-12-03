""" Transport-free depletion test suite """

from pathlib import Path
import shutil

import numpy as np
import pytest
import openmc
import openmc.deplete
from openmc.deplete import IndependentOperator, MicroXS

from tests.regression_tests import config, assert_atoms_equal, \
    assert_reaction_rates_equal, assert_same_mats


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


neutron_per_cm2_sec = 1164719970082145.0


@pytest.mark.parametrize("multiproc, from_nuclides, normalization_mode, power, source_rate", [
    (True, True, 'source-rate', None, 1.0),
    (False, True, 'source-rate', None, 1.0),
    (True, True, 'fission-q', 174, None),
    (False, True, 'fission-q', 174, None),
    (True, False, 'source-rate', None, 1.0),
    (False, False, 'source-rate', None, 1.0),
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
                      source_rate):
    """Transport free system test suite.

    Runs an OpenMC transport-free depletion calculation and verifies
    that the outputs match a reference file.

    """
    # Create operator
    flux = neutron_per_cm2_sec * fuel.volume
    op = _create_operator(from_nuclides,
                          fuel,
                          flux,
                          micro_xs,
                          chain_file,
                          normalization_mode)

    # Power and timesteps
    dt = [360]  # single step

    # Perform simulation using the predictor algorithm
    if config['mpi'] and multiproc:
        pytest.skip("Multiprocessing depletion is disabled when MPI is enabled.")
    openmc.deplete.pool.USE_MULTIPROCESSING = multiproc
    openmc.deplete.PredictorIntegrator(op,
                                       dt,
                                       power=power,
                                       source_rates=source_rate,
                                       timestep_units='s').integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'
    if power is None:
        ref_path = 'test_reference_source_rate.h5'
    else:
        ref_path = 'test_reference_fission_q.h5'
    path_reference = Path(__file__).with_name(ref_path)

    # If updating results, do so and return
    if config['update']:
        shutil.copyfile(str(path_test), str(path_reference))
        return

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)
    res_ref = openmc.deplete.Results(path_reference)

    # Assert same mats
    assert_same_mats(res_ref, res_test)

    tol = 1.0e-14
    assert_atoms_equal(res_ref, res_test, tol)
    assert_reaction_rates_equal(res_ref, res_test, tol)


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
    flux = neutron_per_cm2_sec * fuel.volume
    op = _create_operator(False, fuel, flux, micro_xs, chain_file, 'fission-q')

    # Power and timesteps
    dt = [dt]  # single step

    # Perform simulation using the predictor algorithm
    if config['mpi'] and multiproc:
        pytest.skip("Multiprocessing depletion is disabled when MPI is enabled.")
    openmc.deplete.pool.USE_MULTIPROCESSING = multiproc
    openmc.deplete.PredictorIntegrator(
        op, dt, power=174, timestep_units=time_units).integrate()

    # Get path to test and reference results
    path_test = op.output_dir / 'depletion_results.h5'

    ref_path = f'test_reference_coupled_{time_type}.h5'
    path_reference = Path(__file__).with_name(ref_path)

    # If updating results, do so and return
    if config['update']:
        shutil.copyfile(str(path_test), str(path_reference))
        return

    # Load the reference/test results
    res_test = openmc.deplete.Results(path_test)
    res_ref = openmc.deplete.Results(path_reference)

    # Assert same mats
    assert_same_mats(res_test, res_ref)

    assert_atoms_equal(res_ref, res_test, atom_tol)
    assert_reaction_rates_equal(res_ref, res_test, rx_tol)


def _create_operator(from_nuclides,
                     fuel,
                     flux,
                     micro_xs,
                     chain_file,
                     normalization_mode):
    if from_nuclides:
        nuclides = {}
        for nuc, dens in fuel.get_nuclide_atom_densities().items():
            nuclides[nuc] = dens

        openmc.reset_auto_ids()
        op = IndependentOperator.from_nuclides(fuel.volume,
                                               nuclides,
                                               flux,
                                               micro_xs,
                                               chain_file,
                                               normalization_mode=normalization_mode)

    else:
        op = IndependentOperator(openmc.Materials([fuel]),
                                 [flux],
                                 [micro_xs],
                                 chain_file,
                                 normalization_mode=normalization_mode)

    return op
