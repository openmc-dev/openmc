#!/usr/bin/env python

import glob
import os
import pandas as pd
try:
    from StringIO import StringIO
except:
    from io import StringIO
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from openmc import Filter, Mesh, Tally, TalliesFile, Summary, StatePoint
#from openmc.statepoint import StatePoint
from openmc.source import Source
from openmc.stats import Box

class DiffTallyTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Build default materials/geometry
        self._input_set.build_default_materials_and_geometry()

        # Set settings explicitly
        self._input_set.settings.batches = 5
        self._input_set.settings.inactive = 0
        self._input_set.settings.particles = 400
        self._input_set.settings.source = Source(space=Box(
            [-160, -160, -183], [160, 160, 183]))
        self._input_set.settings.output = {'summary':True}

        self._input_set.tallies = TalliesFile()

        filt_mats = Filter(type='material', bins=(1, 3))
        filt_eout = Filter(type='energyout', bins=(0.0, 1.0, 20.0))

        def add_derivs(tally_list):
            assert len(tally_list) == 4
            # We want density derivatives for both water and fuel to get
            # coverage for both fissile and non-fissile materials.
            tally_list[0].diff_variable = 'density'
            tally_list[0].diff_material = 3
            tally_list[1].diff_variable = 'density'
            tally_list[1].diff_material = 1

            # O-16 is a good nuclide to test against because it is present
            # in both water and fuel.  Some routines need to recognize that they
            # have the perturbed nuclide but not the perturbed material.
            tally_list[2].diff_variable = 'nuclide_density'
            tally_list[2].diff_material = 1
            tally_list[2].diff_nuclide = 'O-16'

            # A fissile nuclide, just for good measure.
            tally_list[3].diff_variable = 'nuclide_density'
            tally_list[3].diff_material = 1
            tally_list[3].diff_nuclide = 'U-235'

        # Cover the flux score.
        tallies = [Tally() for i in range(4)]
        for t in tallies: t.add_score('flux')
        for t in tallies: t.add_filter(filt_mats)
        add_derivs(tallies)
        for t in tallies: self._input_set.tallies.add_tally(t)

        # Cover supported scores with a collision estimator.
        tallies = [Tally() for i in range(4)]
        for t in tallies: t.add_score('total')
        for t in tallies: t.add_score('absorption')
        for t in tallies: t.add_score('fission')
        for t in tallies: t.add_score('nu-fission')
        for t in tallies: t.add_filter(filt_mats)
        for t in tallies: t.add_nuclide('total')
        for t in tallies: t.add_nuclide('U-235')
        add_derivs(tallies)
        for t in tallies: self._input_set.tallies.add_tally(t)

        # Cover an analog estimator.
        tallies = [Tally() for i in range(4)]
        for t in tallies: t.add_score('absorption')
        for t in tallies: t.add_filter(filt_mats)
        for t in tallies: t.estimator = 'analog'
        add_derivs(tallies)
        for t in tallies: self._input_set.tallies.add_tally(t)

        # And the special fission with energyout filter.
        tallies = [Tally() for i in range(4)]
        for t in tallies: t.add_score('nu-fission')
        for t in tallies: t.add_filter(filt_mats)
        for t in tallies: t.add_filter(filt_eout)
        add_derivs(tallies)
        for t in tallies: self._input_set.tallies.add_tally(t)

        self._input_set.export()

    def _get_results(self):
        #return super(DiffTallyTestHarness, self)._get_results(hash_output=True)
        # Read the statepoint and summary files.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = StatePoint(statepoint)
        su = Summary('summary.h5')
        sp.link_with_summary(su)

        # Extract the tally data as a Pandas DataFrame.
        df = pd.DataFrame()
        for t in sp.tallies.values():
            df = df.append(t.get_pandas_dataframe(), ignore_index=True)

        # Extract the relevant data as a CSV string.
        out = StringIO()
        cols = ('d_material', 'd_nuclide', 'd_variable', 'score', 'mean',
                'std. dev.')
        df.to_csv(out, columns=cols, index=False, float_format='%.2e')

        return out.getvalue()

    def _cleanup(self):
        super(DiffTallyTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = DiffTallyTestHarness('statepoint.5.*', True)
    harness.main()
