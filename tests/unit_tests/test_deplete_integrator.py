"""Tests for saving results

It is worth noting that openmc.deplete.integrate is extremely complex, to the
point I am unsure if it can be reasonably unit-tested.  For the time being, it
will be left unimplemented and testing will be done via regression.

"""

import copy
from random import uniform
from unittest.mock import MagicMock

import numpy as np
from uncertainties import ufloat
import pytest

from openmc.deplete import (
    ReactionRates, Results, ResultsList, comm, OperatorResult,
    PredictorIntegrator, CECMIntegrator, CF4Integrator, CELIIntegrator,
    EPCRK4Integrator, LEQIIntegrator, SICELIIntegrator, SILEQIIntegrator)

from tests import dummy_operator


INTEGRATORS = [
    PredictorIntegrator,
    CECMIntegrator,
    CF4Integrator,
    CELIIntegrator,
    EPCRK4Integrator,
    LEQIIntegrator,
    SICELIIntegrator,
    SILEQIIntegrator
]


def test_results_save(run_in_tmpdir):
    """Test data save module"""

    stages = 3

    np.random.seed(comm.rank)

    # Mock geometry
    op = MagicMock()

    # Avoid DummyOperator thinking it's doing a restart calculation
    op.prev_res = None

    vol_dict = {}
    full_burn_list = []

    for i in range(comm.size):
        vol_dict[str(2*i)] = 1.2
        vol_dict[str(2*i + 1)] = 1.2
        full_burn_list.append(str(2*i))
        full_burn_list.append(str(2*i + 1))

    burn_list = full_burn_list[2*comm.rank: 2*comm.rank + 2]
    nuc_list = ["na", "nb"]

    op.get_results_info.return_value = (
        vol_dict, nuc_list, burn_list, full_burn_list)

    # Construct x
    x1 = []
    x2 = []

    for i in range(stages):
        x1.append([np.random.rand(2), np.random.rand(2)])
        x2.append([np.random.rand(2), np.random.rand(2)])

    # Construct r
    r1 = ReactionRates(burn_list, ["na", "nb"], ["ra", "rb"])
    r1[:] = np.random.rand(2, 2, 2)

    rate1 = []
    rate2 = []

    for i in range(stages):
        rate1.append(copy.deepcopy(r1))
        r1[:] = np.random.rand(2, 2, 2)
        rate2.append(copy.deepcopy(r1))
        r1[:] = np.random.rand(2, 2, 2)

    # Create global terms
    # Col 0: eig, Col 1: uncertainty
    eigvl1 = np.random.rand(stages, 2)
    eigvl2 = np.random.rand(stages, 2)

    eigvl1 = comm.bcast(eigvl1, root=0)
    eigvl2 = comm.bcast(eigvl2, root=0)

    t1 = [0.0, 1.0]
    t2 = [1.0, 2.0]

    op_result1 = [OperatorResult(ufloat(*k), rates)
                  for k, rates in zip(eigvl1, rate1)]
    op_result2 = [OperatorResult(ufloat(*k), rates)
                  for k, rates in zip(eigvl2, rate2)]
    Results.save(op, x1, op_result1, t1, 0, 0)
    Results.save(op, x2, op_result2, t2, 0, 1)

    # Load the files
    res = ResultsList.from_hdf5("depletion_results.h5")

    for i in range(stages):
        for mat_i, mat in enumerate(burn_list):
            for nuc_i, nuc in enumerate(nuc_list):
                assert res[0][i, mat, nuc] == x1[i][mat_i][nuc_i]
                assert res[1][i, mat, nuc] == x2[i][mat_i][nuc_i]
        np.testing.assert_array_equal(res[0].rates[i], rate1[i])
        np.testing.assert_array_equal(res[1].rates[i], rate2[i])

    np.testing.assert_array_equal(res[0].k, eigvl1)
    np.testing.assert_array_equal(res[0].time, t1)

    np.testing.assert_array_equal(res[1].k, eigvl2)
    np.testing.assert_array_equal(res[1].time, t2)


def test_bad_integrator_inputs():
    """Test failure modes for Integrator inputs"""

    op = MagicMock()
    op.prev_res = None
    op.chain = None
    op.heavy_metal = 1.0
    timesteps = [1]

    # No power nor power density given
    with pytest.raises(ValueError, match="Either power or power density"):
        PredictorIntegrator(op, timesteps)

    # Length of power != length time
    with pytest.raises(ValueError, match="number of powers"):
        PredictorIntegrator(op, timesteps, power=[1, 2])

    # Length of power density != length time
    with pytest.raises(ValueError, match="number of powers"):
        PredictorIntegrator(op, timesteps, power_density=[1, 2])

    # SI integrator with bad steps
    with pytest.raises(TypeError, match="n_steps"):
        SICELIIntegrator(op, timesteps, [1], n_steps=2.5)

    with pytest.raises(ValueError, match="n_steps"):
        SICELIIntegrator(op, timesteps, [1], n_steps=0)


@pytest.mark.parametrize("scheme", dummy_operator.SCHEMES)
def test_integrator(run_in_tmpdir, scheme):
    """Test the integrators against their expected values"""

    bundle = dummy_operator.SCHEMES[scheme]
    operator = dummy_operator.DummyOperator()
    bundle.solver(operator, [0.75, 0.75], 1.0).integrate()

    # get expected results

    res = ResultsList.from_hdf5(
        operator.output_dir / "depletion_results.h5")

    t1, y1 = res.get_atoms("1", "1")
    t2, y2 = res.get_atoms("1", "2")

    assert (t1 == [0.0, 0.75, 1.5]).all()
    assert y1 == pytest.approx(bundle.atoms_1)
    assert (t2 == [0.0, 0.75, 1.5]).all()
    assert y2 == pytest.approx(bundle.atoms_2)

    # test structure of depletion time dataset
    dep_time = res.get_depletion_time()
    assert dep_time.shape == (2, )
    assert all(dep_time > 0)


@pytest.mark.parametrize("integrator", INTEGRATORS)
def test_timesteps(integrator):
    # Crate fake operator
    op = MagicMock()
    op.prev_res = None
    op.chain = None

    # Set heavy metal mass and power randomly
    op.heavy_metal = uniform(0, 10000)
    power = uniform(0, 1e6)

    # Reference timesteps in seconds
    day = 86400.0
    ref_timesteps = [1*day, 2*day, 5*day, 10*day]

    # Case 1, timesteps in seconds
    timesteps = ref_timesteps
    x = integrator(op, timesteps, power, timestep_units='s')
    assert np.allclose(x.timesteps, ref_timesteps)

    # Case 2, timesteps in minutes
    minute = 60
    timesteps = [t / minute for t in ref_timesteps]
    x = integrator(op, timesteps, power, timestep_units='min')
    assert np.allclose(x.timesteps, ref_timesteps)

    # Case 3, timesteps in hours
    hour = 60*60
    timesteps = [t / hour for t in ref_timesteps]
    x = integrator(op, timesteps, power, timestep_units='h')
    assert np.allclose(x.timesteps, ref_timesteps)

    # Case 4, timesteps in days
    timesteps = [t / day for t in ref_timesteps]
    x = integrator(op, timesteps, power, timestep_units='d')
    assert np.allclose(x.timesteps, ref_timesteps)

    # Case 5, timesteps in MWd/kg
    kilograms = op.heavy_metal / 1000.0
    days = [t/day for t in ref_timesteps]
    megawatts = power / 1000000.0
    burnup = [t * megawatts / kilograms for t in days]
    x = integrator(op, burnup, power, timestep_units='MWd/kg')
    assert np.allclose(x.timesteps, ref_timesteps)

    # Case 6, mixed units
    burnup_per_day = (1e-6*power) / kilograms
    timesteps = [(burnup_per_day, 'MWd/kg'), (2*day, 's'), (5, 'd'),
                 (10*burnup_per_day, 'MWd/kg')]
    x = integrator(op, timesteps, power)
    assert np.allclose(x.timesteps, ref_timesteps)

    # Bad units should raise an exception
    with pytest.raises(ValueError, match="unit"):
        integrator(op, ref_timesteps, power, timestep_units='🐨')
    with pytest.raises(ValueError, match="unit"):
        integrator(op, [(800.0, 'gorillas')], power)
