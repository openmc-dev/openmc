"""Filter weights on delayed-nu-fission tallies with an outgoing energy filter.

With an EnergyoutFilter, delayed-nu-fission is scored once per banked fission
neutron, and a DelayedGroupFilter sends each delayed neutron to its group bin.
The weight of any other filter on the tally, here a SpatialLegendreFilter,
must be applied exactly once. The same delayed neutrons are also tallied
without the EnergyoutFilter and without the DelayedGroupFilter, which are
scored along other paths, so summing the first tally over either of those
filters must reproduce the corresponding tally to round-off.

"""

import os

import numpy as np
import openmc
from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def create_library():
    # Two-group cross sections from the mg_tallies test with two delayed
    # groups added. The delayed neutron fractions are unphysically large so
    # that a short run banks many delayed neutrons, and delayed neutrons are
    # born in both groups so that every outgoing energy bin is scored.
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 20.0e6])
    beta = np.array([0.1, 0.2])
    fission = np.array([0.002817, 0.097])
    nu_fission = 2.5 * fission
    capture = np.array([0.008708, 0.02518])
    scatter = np.array(
        [[[0.31980, 0.06694], [0.004555, -0.0003972]],
         [[0.00000, 0.00000], [0.424100, 0.05439000]]])

    xsdata = openmc.XSdata('mat_1', groups, num_delayed_groups=2)
    xsdata.order = 1
    xsdata.set_total([0.33588, 0.54628])
    xsdata.set_absorption(capture + fission)
    xsdata.set_scatter_matrix(scatter)
    xsdata.set_fission(fission)
    xsdata.set_prompt_nu_fission((1.0 - beta.sum()) * nu_fission)
    xsdata.set_delayed_nu_fission(np.outer(beta, nu_fission))
    xsdata.set_chi_prompt([1.0, 0.0])
    xsdata.set_chi_delayed([[0.8, 0.2], [0.6, 0.4]])
    xsdata.set_decay_rate([0.0127, 0.0317])

    library = openmc.MGXSLibrary(groups, num_delayed_groups=2)
    library.add_xsdata(xsdata)
    library.export_to_hdf5('2g.h5')


class MGDelayedEoutHarness(PyAPITestHarness):
    def _compare_results(self):
        # The tallies must agree with each other whatever the reference
        # results are, so check that first
        with openmc.StatePoint(self._sp_name) as sp:
            eout_dg = sp.get_tally(name='energyout-delayedgroup')
            dg = sp.get_tally(name='delayedgroup')
            eout = sp.get_tally(name='energyout')
            eout_dg = eout_dg.get_reshaped_data(value='mean')[..., 0, 0]
            dg = dg.get_reshaped_data(value='mean')[..., 0, 0]
            eout = eout.get_reshaped_data(value='mean')[..., 0, 0]

        # Every bin must be scored for the comparisons to mean anything
        assert np.all(dg > 0.0)
        assert np.all(eout > 0.0)

        # Indices of eout_dg are (energyout, delayed group, Legendre order)
        np.testing.assert_allclose(eout_dg.sum(axis=0), dg, rtol=1e-10)
        np.testing.assert_allclose(eout_dg.sum(axis=1), eout, rtol=1e-10)

        super()._compare_results()

    def _cleanup(self):
        super()._cleanup()
        f = '2g.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mg_delayed_eout_weights():
    create_library()
    model = slab_mg()

    # The slab spans 0 <= x <= 929.45 cm. With the filter defined over twice
    # that width, the P1 weight is x/929.45, which lies between 0 and 1 in the
    # slab, so the P1 results are sums of non-negative terms that can be
    # compared with a relative tolerance.
    legendre = openmc.SpatialLegendreFilter(1, 'x', -929.45, 929.45)
    energyout = openmc.EnergyoutFilter([0.0, 0.625, 20.0e6])
    delayed_groups = openmc.DelayedGroupFilter([1, 2])

    for name, filters in [
        ('energyout-delayedgroup', [energyout, delayed_groups, legendre]),
        ('delayedgroup', [delayed_groups, legendre]),
        ('energyout', [energyout, legendre]),
    ]:
        tally = openmc.Tally(name=name)
        tally.filters = filters
        tally.scores = ['delayed-nu-fission']
        # Required: the SpatialLegendreFilter would otherwise give the
        # delayedgroup tally a collision estimator, and the checks rely on all
        # three tallies scoring the same analog events
        tally.estimator = 'analog'
        model.tallies.append(tally)

    harness = MGDelayedEoutHarness('statepoint.10.h5', model)
    harness.main()
