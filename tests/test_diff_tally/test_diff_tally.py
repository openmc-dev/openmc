#!/usr/bin/env python

import glob
import os
try:
    from StringIO import StringIO
except:
    from io import StringIO
import sys

import pandas as pd

sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from openmc import Filter, Mesh, Tally, Tallies, Summary, StatePoint, \
                   TallyDerivative
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
        self._input_set.settings.output = {'summary': True}

        self._input_set.tallies = Tallies()

        filt_mats = Filter(type='material', bins=(1, 3))
        filt_eout = Filter(type='energyout', bins=(0.0, 1.0, 20.0))

        # We want density derivatives for both water and fuel to get coverage
        # for both fissile and non-fissile materials.
        d1 = TallyDerivative(derivative_id=1)
        d1.variable = 'density'
        d1.material = 3
        d2 = TallyDerivative(derivative_id=2)
        d2.variable = 'density'
        d2.material = 1

        # O-16 is a good nuclide to test against because it is present in both
        # water and fuel.  Some routines need to recognize that they have the
        # perturbed nuclide but not the perturbed material.
        d3 = TallyDerivative(derivative_id=3)
        d3.variable = 'nuclide_density'
        d3.material = 1
        d3.nuclide = 'O16'

        # A fissile nuclide, just for good measure.
        d4 = TallyDerivative(derivative_id=4)
        d4.variable = 'nuclide_density'
        d4.material = 1
        d4.nuclide = 'U235'

        derivs = [d1, d2, d3, d4]

        # Cover the flux score.
        for i in range(4):
            t = Tally()
            t.add_score('flux')
            t.add_filter(filt_mats)
            t.derivative = derivs[i]
            self._input_set.tallies.add_tally(t)

        # Cover supported scores with a collision estimator.
        for i in range(4):
            t = Tally()
            t.add_score('total')
            t.add_score('absorption')
            t.add_score('fission')
            t.add_score('nu-fission')
            t.add_filter(filt_mats)
            t.add_nuclide('total')
            t.add_nuclide('U235')
            t.derivative = derivs[i]
            self._input_set.tallies.add_tally(t)

        # Cover an analog estimator.
        for i in range(4):
            t = Tally()
            t.add_score('absorption')
            t.add_filter(filt_mats)
            t.estimator = 'analog'
            t.derivative = derivs[i]
            self._input_set.tallies.add_tally(t)

        # And the special fission with energyout filter.
        for i in range(4):
            t = Tally()
            t.add_score('nu-fission')
            t.add_filter(filt_mats)
            t.add_filter(filt_eout)
            t.derivative = derivs[i]
            self._input_set.tallies.add_tally(t)

        self._input_set.export()

    def _get_results(self):
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
