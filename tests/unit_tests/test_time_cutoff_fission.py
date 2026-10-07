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
