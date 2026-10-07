"""Analog fission tallies with a neutron time cutoff.

A delayed fission neutron that would be emitted after the time cutoff is not
banked, since it would be killed at birth, but it was still produced by the
fission. Analog fission tallies must count it, as the track-length and
collision estimators do, and must keep the record of the fission neutrons of a
collision consistent with the number of neutrons they loop over.

"""

from pathlib import Path

import numpy as np
import openmc
import pytest
from openmc.examples import slab_mg


T_CUTOFF = 1.0e-2


def create_library(all_delayed=False):
    """Two-group cross sections with two delayed groups.

    The delayed neutron fractions are unphysically large so that a short run
    produces many delayed neutrons, nearly all of which are emitted after the
    time cutoff. With all_delayed, every fission neutron is delayed.

    """
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 20.0e6])
    fission = np.array([0.002817, 0.097])
    nu_fission = 2.5 * fission
    capture = np.array([0.008708, 0.02518])
    scatter = np.array(
        [[[0.31980, 0.06694], [0.004555, -0.0003972]],
         [[0.00000, 0.00000], [0.424100, 0.05439000]]])
    if all_delayed:
        beta = np.array([1/3, 2/3])
        prompt_nu_fission = np.zeros(2)
    else:
        beta = np.array([0.1, 0.2])
        prompt_nu_fission = (1.0 - beta.sum()) * nu_fission

    xsdata = openmc.XSdata('mat_1', groups, num_delayed_groups=2)
    xsdata.order = 1
    xsdata.set_total([0.33588, 0.54628])
    xsdata.set_absorption(capture + fission)
    xsdata.set_scatter_matrix(scatter)
    xsdata.set_fission(fission)
    xsdata.set_prompt_nu_fission(prompt_nu_fission)
    xsdata.set_delayed_nu_fission(np.outer(beta, nu_fission))
    xsdata.set_chi_prompt([1.0, 0.0])
    xsdata.set_chi_delayed([[0.8, 0.2], [0.6, 0.4]])
    xsdata.set_decay_rate([0.0127, 0.0317])
    xsdata.set_inverse_velocity([1.0e-9, 4.5e-6])

    library = openmc.MGXSLibrary(groups, num_delayed_groups=2)
    library.add_xsdata(xsdata)
    library.export_to_hdf5('2g.h5')
    return Path('2g.h5').resolve()


def mg_model(all_delayed=False):
    # The library path is absolute because the model may run in another
    # directory
    model = slab_mg(mgxslib_name=str(create_library(all_delayed)))
    model.settings.seed = 1
    return model


def ce_model():
    model = openmc.Model()
    material = openmc.Material()
    material.add_nuclide('U235', 1.0)
    material.set_density('g/cm3', 16.0)
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=material)
    model.geometry = openmc.Geometry([cell])
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(*cell.bounding_box),
        constraints={'fissionable': True})
    model.settings.seed = 1
    return model


def add_tallies(model, energyout_edges, delayed_groups):
    scores = ['nu-fission', 'prompt-nu-fission', 'delayed-nu-fission']
    for name, filters, tally_scores in [
        ('analog', [], scores),
        ('energyout', [openmc.EnergyoutFilter(energyout_edges)], scores),
        ('delayed', [openmc.DelayedGroupFilter(delayed_groups)],
         ['delayed-nu-fission', 'decay-rate']),
    ]:
        tally = openmc.Tally(name=name)
        tally.filters = filters
        tally.scores = tally_scores
        tally.estimator = 'analog'
        model.tallies.append(tally)
    tally = openmc.Tally(name='flux')
    tally.scores = ['flux']
    tally.estimator = 'tracklength'
    model.tallies.append(tally)


def get_mean(sp, name, score):
    tally = sp.get_tally(name=name)
    return tally.get_values(scores=[score]).ravel()


def assert_equal(a, b):
    """Equal to round-off in the sums of the tally results."""
    a = np.asarray(a)
    b = np.asarray(b)
    assert np.abs(a - b).max() <= 1e-10 * np.abs(b).max()


@pytest.mark.parametrize('event_based', [False, True], ids=['history', 'event'])
@pytest.mark.parametrize('energy_mode', ['mg', 'ce'])
def test_time_cutoff_eigenvalue(run_in_tmpdir, energy_mode, event_based):
    """The first generation with and without a time cutoff.

    The source particles are transported the same way whether or not the
    cutoff is set, because no history lasts as long as the cutoff, so analog
    tallies of the fission neutrons produced must be the same. Only the
    fission bank differs: with the cutoff, delayed neutrons emitted after it
    are not banked.

    """
    if energy_mode == 'mg':
        model = mg_model()
        model.settings.particles = 1000
        energyout_edges = [0.0, 0.625, 20.0e6]
        delayed_groups = [1, 2]
    else:
        model = ce_model()
        model.settings.particles = 10000
        # Fission neutron energies are resampled below the maximum energy of
        # the data, so these bins cover every fission neutron
        energyout_edges = [0.0, 1.0e6, 1.0e8]
        delayed_groups = [1, 2, 3, 4, 5, 6]
    model.settings.batches = 1
    model.settings.inactive = 0
    model.settings.event_based = event_based
    add_tallies(model, energyout_edges, delayed_groups)

    results = {}
    for run in ('no-cutoff', 'cutoff'):
        if run == 'cutoff':
            model.settings.cutoff = {'time_neutron': T_CUTOFF}
        sp_path = model.run(cwd=run)
        with openmc.StatePoint(sp_path) as sp:
            res = {'k': sp.k_generation, 'source': sp.source}
            for name in ('flux', 'analog', 'energyout', 'delayed'):
                for score in sp.get_tally(name=name).scores:
                    res[name, score] = get_mean(sp, name, score)
            results[run] = res
    a = results['cutoff']
    b = results['no-cutoff']

    # The transport is the same in both runs
    assert_equal(a['k'], b['k'])
    assert_equal(a['flux', 'flux'], b['flux', 'flux'])

    # Delayed neutrons emitted after the cutoff are produced in both runs, but
    # banked only without the cutoff
    source = b['source']
    assert np.any((source['delayed_group'] > 0) & (source['time'] > T_CUTOFF))
    assert np.all(a['source']['time'] <= T_CUTOFF)

    for key in a:
        if key[0] in ('analog', 'energyout', 'delayed'):
            assert_equal(a[key], b[key])

    for res in (a, b):
        # The tally with an EnergyoutFilter loops over the fission neutrons of
        # each collision, so it must sum to the tally without the filter
        for score in ('nu-fission', 'prompt-nu-fission', 'delayed-nu-fission'):
            assert_equal(res['energyout', score].sum(), res['analog', score])

        if energy_mode == 'mg':
            # The decay rate loops over the fission neutrons of each collision,
            # and the delayed-nu-fission score counts them
            decay_rate = np.array([0.0127, 0.0317])
            assert_equal(res['delayed', 'decay-rate'],
                         decay_rate * res['delayed', 'delayed-nu-fission'])

    if energy_mode == 'mg':
        assert np.all(a['delayed', 'delayed-nu-fission'] > 0.0)
    else:
        assert a['analog', 'delayed-nu-fission'].sum() > 0.0


@pytest.mark.parametrize('event_based', [False, True], ids=['history', 'event'])
def test_time_cutoff_all_delayed(run_in_tmpdir, event_based):
    """No fission neutron is prompt, so no prompt-nu-fission may be scored."""
    model = mg_model(all_delayed=True)
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000
    model.settings.event_based = event_based
    model.settings.cutoff = {'time_neutron': T_CUTOFF}
    scores = ['nu-fission', 'prompt-nu-fission', 'delayed-nu-fission']
    for name, filters in [
        ('analog', []),
        ('energyout', [openmc.EnergyoutFilter([0.0, 0.625, 20.0e6])]),
    ]:
        tally = openmc.Tally(name=name)
        tally.filters = filters
        tally.scores = scores
        tally.estimator = 'analog'
        model.tallies.append(tally)

    sp_path = model.run()
    with openmc.StatePoint(sp_path) as sp:
        for name in ('analog', 'energyout'):
            nu_fission = get_mean(sp, name, 'nu-fission')
            prompt = get_mean(sp, name, 'prompt-nu-fission')
            delayed = get_mean(sp, name, 'delayed-nu-fission')
            assert np.all(prompt == 0.0)
            assert_equal(delayed, nu_fission)
            assert np.all(delayed > 0.0)
        for score in scores:
            assert_equal(get_mean(sp, 'energyout', score).sum(),
                         get_mean(sp, 'analog', score))


# Decay rates so small that every delayed neutron is emitted long after the
# time cutoff
SLOW_DECAY_RATE = np.array([1.0e-6, 2.0e-6])


def create_exact_library(nu, prompt_fraction):
    """Two-group cross sections in which a collision produces nu neutrons.

    nu-fission / total equals nu in both groups and every number is a binary
    fraction, so a collision of a particle of unit weight samples exactly nu
    fission neutrons while keff is 1, i.e., in a fixed-source run and in the
    first generation of an eigenvalue run.

    """
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 20.0e6])
    nu_fission = np.full(2, 0.5 * nu)
    xsdata = openmc.XSdata('mat_1', groups, num_delayed_groups=2)
    xsdata.order = 0
    xsdata.set_total([0.5, 0.5])
    xsdata.set_absorption([0.375, 0.375])
    xsdata.set_scatter_matrix([[[0.0625], [0.0625]], [[0.0], [0.125]]])
    xsdata.set_fission([0.25, 0.25])
    xsdata.set_prompt_nu_fission(prompt_fraction * nu_fission)
    xsdata.set_delayed_nu_fission(
        np.outer([0.5, 0.5], (1.0 - prompt_fraction) * nu_fission))
    xsdata.set_chi_prompt([1.0, 0.0])
    xsdata.set_chi_delayed([[0.75, 0.25], [0.5, 0.5]])
    xsdata.set_decay_rate(SLOW_DECAY_RATE)
    xsdata.set_inverse_velocity([1.0e-9, 4.5e-6])

    library = openmc.MGXSLibrary(groups, num_delayed_groups=2)
    library.add_xsdata(xsdata)
    library.export_to_hdf5('exact.h5')
    return Path('exact.h5').resolve()


@pytest.mark.parametrize('event_based', [False, True], ids=['history', 'event'])
def test_time_cutoff_one_neutron_per_collision(run_in_tmpdir, event_based):
    """Every collision produces one delayed neutron, which is never banked.

    Each collision of a fixed-source run produces exactly one fission neutron,
    which is delayed and emitted long after the time cutoff. The analog
    estimates must then equal the collision estimates collision by collision,
    and no fission neutron may enter the simulation or change the particle
    that produced it.

    """
    model = slab_mg(mgxslib_name=str(create_exact_library(1, 0.0)))
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 2
    model.settings.inactive = 0
    model.settings.particles = 1000
    model.settings.seed = 1
    model.settings.event_based = event_based
    model.settings.cutoff = {'time_neutron': T_CUTOFF}
    model.settings.collision_track = {'max_collisions': 100000}

    scores = ['nu-fission', 'prompt-nu-fission', 'delayed-nu-fission']
    for name, filters, tally_scores, estimator in [
        ('analog', [], scores, 'analog'),
        ('collision', [], ['nu-fission', 'delayed-nu-fission'], 'collision'),
        ('energyout', [openmc.EnergyoutFilter([0.0, 0.625, 20.0e6])], scores,
         'analog'),
        ('delayed', [openmc.DelayedGroupFilter([1, 2])],
         ['delayed-nu-fission', 'decay-rate'], 'analog'),
        ('flux', [], ['flux'], 'tracklength'),
        ('production', [openmc.ParticleProductionFilter('neutron')],
         ['events'], 'analog'),
    ]:
        tally = openmc.Tally(name=name)
        tally.filters = filters
        tally.scores = tally_scores
        tally.estimator = estimator
        model.tallies.append(tally)

    sp_path = model.run()
    with openmc.StatePoint(sp_path) as sp:
        nu_fission = get_mean(sp, 'analog', 'nu-fission')
        assert np.all(nu_fission > 0.0)

        # One neutron per collision, all of them delayed
        assert_equal(nu_fission, get_mean(sp, 'collision', 'nu-fission'))
        assert_equal(get_mean(sp, 'collision', 'delayed-nu-fission'),
                     nu_fission)
        assert_equal(get_mean(sp, 'analog', 'delayed-nu-fission'), nu_fission)
        assert np.all(get_mean(sp, 'analog', 'prompt-nu-fission') == 0.0)
        assert np.all(get_mean(sp, 'energyout', 'prompt-nu-fission') == 0.0)
        for score in scores:
            assert_equal(get_mean(sp, 'energyout', score).sum(),
                         get_mean(sp, 'analog', score))
        delayed = get_mean(sp, 'delayed', 'delayed-nu-fission')
        assert_equal(delayed.sum(), nu_fission)
        assert_equal(get_mean(sp, 'delayed', 'decay-rate'),
                     SLOW_DECAY_RATE * delayed)

        # No fission neutron is transported
        assert np.all(get_mean(sp, 'production', 'events') == 0.0)
        assert np.all(get_mean(sp, 'flux', 'flux') > 0.0)

    # The particle that produced the fission neutrons is unchanged by them: it
    # keeps the delayed group of its source site and takes no progeny
    collisions = openmc.read_collision_track_hdf5('collision_track.h5')
    assert collisions.size > 0
    assert np.all(collisions['delayed_group'] == 0)
    assert np.all(collisions['progeny_id'] == 0)


@pytest.mark.parametrize('event_based', [False, True], ids=['history', 'event'])
def test_time_cutoff_fixed_source_ce(run_in_tmpdir, event_based):
    """Analog and track-length estimates agree in a subcritical sphere.

    Only fission neutrons emitted before the cutoff may enter the simulation:
    one emitted after it would be moved back to the cutoff and score a
    negative track length, so the flux must be positive. The analog estimates
    of the fission neutrons produced must agree with the track-length
    estimates, which count every neutron produced whatever the cutoff.

    """
    model = openmc.Model()
    material = openmc.Material()
    material.add_nuclide('U235', 1.0)
    material.set_density('g/cm3', 16.0)
    sphere = openmc.Sphere(r=8.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=material)
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(), energy=openmc.stats.delta_function(2.0e6))
    model.settings.batches = 10
    model.settings.particles = 1000
    model.settings.seed = 1
    model.settings.event_based = event_based
    model.settings.cutoff = {'time_neutron': T_CUTOFF}

    scores = ['nu-fission', 'prompt-nu-fission', 'delayed-nu-fission']
    for name, filters, tally_scores, estimator in [
        ('analog', [], scores, 'analog'),
        ('tracklength', [], scores + ['flux'], 'tracklength'),
        ('energyout', [openmc.EnergyoutFilter([0.0, 1.0e6, 1.0e8])], scores,
         'analog'),
    ]:
        tally = openmc.Tally(name=name)
        tally.filters = filters
        tally.scores = tally_scores
        tally.estimator = estimator
        model.tallies.append(tally)

    sp_path = model.run()
    with openmc.StatePoint(sp_path) as sp:
        assert np.all(get_mean(sp, 'tracklength', 'flux') > 0.0)
        for score in scores:
            analog = get_mean(sp, 'analog', score)
            tracklength = get_mean(sp, 'tracklength', score)
            # The relative standard deviation of the analog estimate is below
            # 6% for delayed neutrons and below 3% otherwise
            assert np.all(np.abs(analog / tracklength - 1.0) < 0.5)
            assert_equal(get_mean(sp, 'energyout', score).sum(), analog)


@pytest.mark.parametrize('event_based', [False, True], ids=['history', 'event'])
@pytest.mark.parametrize('energy_mode', ['mg', 'ce'])
def test_time_cutoff_full_fission_bank(run_in_tmpdir, energy_mode,
                                       event_based):
    """The fission bank fills up in the first generation.

    The fission bank holds three sites per source particle. The site that
    fails to enter the full bank and the later sites of its collision are
    neither banked nor counted, so analog tallies count exactly the sites in
    the bank, plus the delayed neutrons that the time cutoff rejected before
    that site. With the multigroup library and the cutoff, every delayed
    neutron is rejected, so the bank holds prompt neutrons only.

    """
    if energy_mode == 'mg':
        # Every collision produces eight neutrons and the bank size, 1500, is
        # not a multiple of eight, so the bank fills in the middle of a
        # collision
        model = slab_mg(mgxslib_name=str(create_exact_library(8, 0.5)))
        runs = ('no-cutoff', 'cutoff')
    else:
        # k-infinity is about 4, so the fission sites of the first generation
        # do not fit in the bank
        model = openmc.Model()
        material = openmc.Material()
        material.add_nuclide('Cf249', 1.0)
        material.set_density('g/cm3', 15.0)
        box = openmc.model.RectangularParallelepiped(
            -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='reflective')
        cell = openmc.Cell(region=-box, fill=material)
        model.geometry = openmc.Geometry([cell])
        model.settings.source = openmc.IndependentSource(
            space=openmc.stats.Box(*cell.bounding_box))
        # A delayed neutron of the evaluated data can be emitted before the
        # cutoff and banked, so the count is exact only without a cutoff
        runs = ('no-cutoff',)
    model.settings.batches = 1
    model.settings.inactive = 0
    model.settings.particles = 500
    model.settings.seed = 1
    model.settings.event_based = event_based
    scores = ['nu-fission', 'prompt-nu-fission', 'delayed-nu-fission']
    tally = openmc.Tally(name='analog')
    tally.scores = scores
    tally.estimator = 'analog'
    model.tallies.append(tally)

    for run in runs:
        if run == 'cutoff':
            model.settings.cutoff = {'time_neutron': T_CUTOFF}
        sp_path = model.run(cwd=run)
        with openmc.StatePoint(sp_path) as sp:
            nu_fission, prompt, delayed = (
                get_mean(sp, 'analog', score) for score in scores)
        assert_equal(prompt + delayed, nu_fission)
        if run == 'no-cutoff':
            # The tallies are normalized by the number of source particles
            assert_equal(nu_fission, [3.0])
        else:
            assert_equal(prompt, [3.0])
            assert np.all(delayed > 0.0)


@pytest.mark.parametrize('event_based', [False, True], ids=['history', 'event'])
@pytest.mark.parametrize('energy_mode', ['mg', 'ce'])
def test_time_cutoff_ufs(run_in_tmpdir, energy_mode, event_based):
    """Fission neutrons whose weight is set by uniform fission site weighting.

    From the second generation on, a fission neutron has the weight 1/w of
    its collision's UFS weight w, so the analog tallies that loop over the
    fission neutrons of a collision must sum to the tallies that use the
    total weight of the collision's fission neutrons.

    """
    mesh = openmc.RegularMesh()
    if energy_mode == 'mg':
        model = mg_model()
        model.settings.particles = 1000
        mesh.lower_left = (0.0, -1.0e6, -1.0e6)
        mesh.upper_right = (929.45, 1.0e6, 1.0e6)
        mesh.dimension = (4, 1, 1)
        energyout_edges = [0.0, 0.625, 20.0e6]
        delayed_groups = [1, 2]
    else:
        model = ce_model()
        model.settings.particles = 2000
        mesh.lower_left = (-10.0, -10.0, -10.0)
        mesh.upper_right = (10.0, 10.0, 10.0)
        mesh.dimension = (3, 3, 3)
        energyout_edges = [0.0, 1.0e6, 1.0e8]
        delayed_groups = [1, 2, 3, 4, 5, 6]
    model.settings.batches = 3
    model.settings.inactive = 1
    model.settings.event_based = event_based
    model.settings.cutoff = {'time_neutron': T_CUTOFF}
    model.settings.ufs_mesh = mesh
    add_tallies(model, energyout_edges, delayed_groups)

    sp_path = model.run()
    with openmc.StatePoint(sp_path) as sp:
        # The weights are not all one
        assert np.any(sp.source['wgt'] != 1.0)
        for score in ('nu-fission', 'prompt-nu-fission', 'delayed-nu-fission'):
            assert_equal(get_mean(sp, 'energyout', score).sum(),
                         get_mean(sp, 'analog', score))
        delayed = get_mean(sp, 'delayed', 'delayed-nu-fission')
        assert delayed.sum() > 0.0
        assert_equal(delayed.sum(),
                     get_mean(sp, 'analog', 'delayed-nu-fission'))
        if energy_mode == 'mg':
            decay_rate = np.array([0.0127, 0.0317])
            assert_equal(get_mean(sp, 'delayed', 'decay-rate'),
                         decay_rate * delayed)
