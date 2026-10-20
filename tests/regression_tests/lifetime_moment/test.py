"""Analytic test of the lifetime moment filter in an infinite medium.

A neutron born at time zero in an infinite, homogeneous, one-group medium with
speed v, total cross section Sigma_t, and absorption cross section Sigma_a
collides at the rate v Sigma_t and is still alive at time tau with probability
exp(-v Sigma_a tau). A collision estimator scores Sigma_x / Sigma_t at every
collision, so the expected score Sigma_x weighted by tau^n per source neutron
is

    M_n(x) = (Sigma_x / Sigma_t) int_0^inf tau^n v Sigma_t exp(-v Sigma_a tau)
             dtau = Sigma_x n! / (Sigma_a^(n+1) v^n),

which gives M_n(nu-fission) = nu Sigma_f n! / (Sigma_a^(n+1) v^n) and
M_n(absorption) = n! / (Sigma_a v)^n. Analog absorption and fission estimators
score at the absorption time, which is exponentially distributed with rate
v Sigma_a, and have the same expectations; because every history ends with one
absorption at unit weight, their zeroth moments are exact. In the fixed-source
calculation the source neutrons are emitted uniformly over 1 ms, much longer
than their mean lifetime 1 / (v Sigma_a) = 40 us; the moments of the time since
birth are unaffected, but moments of the absolute time would not match the
closed forms. In an eigenvalue calculation the same expressions hold per source
neutron because fission neutrons are banked rather than tracked, with
k = nu Sigma_f / Sigma_a and the prompt and delayed parts in the proportions
1 - beta and beta.

The tallies are checked against these closed forms to within four standard
deviations instead of being hashed and compared with a results_true.dat file.
Reference results of hash-based tests are only portable when they are generated
with a strict floating-point build (-DOPENMC_ENABLE_STRICT_FP=on), whereas the
closed forms check the filter weights quantitatively with any build. The
eigenvalue test also checks that analog fission tallies with an outgoing energy
filter, which are scored by a separate routine for each banked fission neutron
(with and without a delayed group filter), agree with the same tallies without
that filter to round-off.

"""
from math import factorial

import numpy as np
import openmc
import pytest

from tests.regression_tests import config

# Some estimates are exact in this problem, and the sample standard deviation
# of a constant can evaluate to the square root of a slightly negative number
pytestmark = pytest.mark.filterwarnings(
    'ignore:invalid value encountered in sqrt:RuntimeWarning')

# One-group data. The large delayed neutron fraction keeps the statistical
# uncertainty of the delayed tallies small.
INVERSE_VELOCITY = 1.0e-5  # s/cm
SIGMA_T = 1.0  # 1/cm
SIGMA_A = 0.25  # 1/cm
SIGMA_F = 0.12  # 1/cm
NU = 2.5
BETA = 0.25
DECAY_RATE = 0.1  # 1/s
ORDER = 2

NU_SIGMA_F = NU * SIGMA_F
K_INF = NU_SIGMA_F / SIGMA_A


def lifetime_moments(sigma):
    """Closed-form moments M_n, n = 0, ..., ORDER, of a reaction rate with
    macroscopic cross section sigma per source neutron."""
    v = 1.0 / INVERSE_VELOCITY
    return np.array([sigma * factorial(n) / (SIGMA_A**(n + 1) * v**n)
                     for n in range(ORDER + 1)])


def make_model(library):
    """Reflective box of a homogeneous one-group medium with one delayed
    group."""
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    xsdata = openmc.XSdata('medium', groups, num_delayed_groups=1)
    xsdata.order = 0
    xsdata.set_total([SIGMA_T])
    xsdata.set_absorption([SIGMA_A])
    xsdata.set_fission([SIGMA_F])
    xsdata.set_nu_fission([NU_SIGMA_F])
    xsdata.set_scatter_matrix([[[SIGMA_T - SIGMA_A]]])
    xsdata.set_chi([1.0])
    xsdata.set_beta([BETA])
    xsdata.set_decay_rate([DECAY_RATE])
    xsdata.set_inverse_velocity([INVERSE_VELOCITY])
    mgxs_library = openmc.MGXSLibrary(groups, num_delayed_groups=1)
    mgxs_library.add_xsdata(xsdata)
    mgxs_library.export_to_hdf5(library)

    medium = openmc.Material()
    medium.set_density('macro', 1.0)
    medium.add_macroscopic('medium')

    model = openmc.Model()
    model.materials = openmc.Materials([medium])
    model.materials.cross_sections = str(library)
    box = openmc.model.RectangularParallelepiped(
        -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='reflective')
    model.geometry = openmc.Geometry([openmc.Cell(fill=medium, region=-box)])

    model.settings.energy_mode = 'multi-group'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))
    return model


def run_model(model, path):
    kwargs = {'cwd': path, 'openmc_exec': config['exe'],
              'event_based': config['event']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    return model.run(**kwargs)


def check_moments(tally, score, sigma, exact_zeroth=False, rel_err=0.02):
    """Check the moments of a score against the closed form.

    The moments must agree to within four standard deviations, and each
    standard deviation must be small enough for the check to be meaningful.
    With exact_zeroth, the zeroth moment has no statistical uncertainty and
    must agree to round-off.

    """
    mean = tally.get_values(scores=[score]).ravel()
    std_dev = tally.get_values(scores=[score], value='std_dev').ravel()
    expected = lifetime_moments(sigma)
    assert mean.size == ORDER + 1
    first = 0
    if exact_zeroth:
        assert mean[0] == pytest.approx(expected[0], rel=1.0e-12)
        first = 1
    mean, std_dev, expected = mean[first:], std_dev[first:], expected[first:]
    assert np.all(np.abs(mean - expected) <= 4.0 * std_dev), \
        f'{score}: {mean} differs from {expected} (std. dev. {std_dev})'
    assert np.all(std_dev <= rel_err * expected)


def test_lifetime_moment_fixed_source(tmp_path):
    model = make_model(tmp_path / 'one_group.h5')
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 50
    model.settings.particles = 10000
    model.settings.create_fission_neutrons = False
    model.settings.source[0].time = openmc.stats.Uniform(0.0, 1.0e-3)

    # Default estimator, which the filter turns into a collision estimator
    collision = openmc.Tally()
    collision.filters = [openmc.LifetimeMomentFilter(ORDER)]
    collision.scores = ['nu-fission', 'absorption', 'fission', 'flux']

    analog = openmc.Tally()
    analog.filters = [openmc.LifetimeMomentFilter(ORDER)]
    analog.scores = ['absorption', 'fission']
    analog.estimator = 'analog'
    model.tallies = [collision, analog]

    with openmc.StatePoint(run_model(model, tmp_path)) as sp:
        t_collision = sp.tallies[collision.id]
        t_analog = sp.tallies[analog.id]

    # The filter and the selected estimators are read back from the statepoint
    assert t_collision.estimator == 'collision'
    assert t_analog.estimator == 'analog'
    lifetime_filter = t_collision.find_filter(openmc.LifetimeMomentFilter)
    assert lifetime_filter.order == ORDER
    assert lifetime_filter.bins == [f'tau^{n}' for n in range(ORDER + 1)]
    df = t_collision.get_pandas_dataframe()
    assert list(df['lifetimemoment'].unique()) == lifetime_filter.bins

    check_moments(t_collision, 'nu-fission', NU_SIGMA_F)
    check_moments(t_collision, 'absorption', SIGMA_A)
    check_moments(t_collision, 'fission', SIGMA_F)
    check_moments(t_collision, 'flux', 1.0)
    check_moments(t_analog, 'absorption', SIGMA_A, exact_zeroth=True)
    check_moments(t_analog, 'fission', SIGMA_F, exact_zeroth=True)


def test_lifetime_moment_eigenvalue(tmp_path):
    model = make_model(tmp_path / 'one_group.h5')
    model.settings.run_mode = 'eigenvalue'
    model.settings.batches = 60
    model.settings.inactive = 10
    model.settings.particles = 10000

    group_filter = openmc.DelayedGroupFilter([1])
    energyout_filter = openmc.EnergyoutFilter([0.0, 20.0e6])

    def lifetime_tally(filters, scores, estimator=None):
        tally = openmc.Tally()
        tally.filters = filters + [openmc.LifetimeMomentFilter(ORDER)]
        tally.scores = scores
        tally.estimator = estimator
        return tally

    fission_scores = ['nu-fission', 'prompt-nu-fission']
    collision = lifetime_tally([], fission_scores + ['absorption'])
    collision_dg = lifetime_tally([group_filter], ['delayed-nu-fission'])
    analog = lifetime_tally([], fission_scores, 'analog')
    analog_dg = lifetime_tally([group_filter], ['delayed-nu-fission'],
                               'analog')
    # The outgoing energy filter requires an analog estimator
    eout = lifetime_tally([energyout_filter], fission_scores)
    eout_dg = lifetime_tally([energyout_filter, group_filter],
                             ['delayed-nu-fission'])
    tallies = [collision, collision_dg, analog, analog_dg, eout, eout_dg]
    model.tallies = tallies

    with openmc.StatePoint(run_model(model, tmp_path)) as sp:
        k_collision = sp.global_tallies[0]
        results = {t.id: sp.tallies[t.id] for t in tallies}

    # The absorption estimate of k is exact here, so check the collision one
    assert k_collision['name'] == b'k-collision'
    assert abs(k_collision['mean'] - K_INF) <= 4.0 * k_collision['std_dev']

    assert results[collision.id].estimator == 'collision'
    assert results[collision_dg.id].estimator == 'collision'
    for t in (analog, analog_dg, eout, eout_dg):
        assert results[t.id].estimator == 'analog'

    for t in (collision, analog, eout):
        check_moments(results[t.id], 'nu-fission', NU_SIGMA_F)
        check_moments(results[t.id], 'prompt-nu-fission',
                      (1.0 - BETA) * NU_SIGMA_F)
    for t in (collision_dg, analog_dg, eout_dg):
        check_moments(results[t.id], 'delayed-nu-fission', BETA * NU_SIGMA_F)
    check_moments(results[collision.id], 'absorption', SIGMA_A)

    # Scoring each banked fission neutron separately gives the same analog
    # estimate as scoring the banked weight at once
    np.testing.assert_allclose(results[eout.id].mean, results[analog.id].mean,
                               rtol=1.0e-10)
    np.testing.assert_allclose(results[eout_dg.id].mean,
                               results[analog_dg.id].mean, rtol=1.0e-10)
