#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc

class FilterPolarTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt1 = openmc.Filter(type='polar',
                              bins=(0.0, 0.6283, 1.2566, 1.8850, 2.5132,
                                    3.1416))
        tally1 = openmc.Tally(tally_id=1)
        tally1.add_filter(filt1)
        tally1.add_score('flux')
        tally1.estimator = 'tracklength'

        tally2 = openmc.Tally(tally_id=2)
        tally2.add_filter(filt1)
        tally2.add_score('flux')
        tally2.estimator = 'analog'

        filt3 = openmc.Filter(type='polar', bins=(5,))
        tally3 = openmc.Tally(tally_id=3)
        tally3.add_filter(filt3)
        tally3.add_score('flux')
        tally3.estimator = 'tracklength'

        mesh = openmc.Mesh(mesh_id=1)
        mesh.lower_left  = [-182.07, -182.07]
        mesh.upper_right = [182.07,  182.07]
        mesh.dimension = [2, 2]
        filt_mesh = openmc.Filter(type='mesh', bins=(1,))
        tally4 = openmc.Tally(tally_id=4)
        tally4.add_filter(filt3)
        tally4.add_filter(filt_mesh)
        tally4.add_score('flux')
        tally4.estimator = 'tracklength'


        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(tally1)
        self._input_set.tallies.add_tally(tally2)
        self._input_set.tallies.add_tally(tally3)
        self._input_set.tallies.add_tally(tally4)
        self._input_set.tallies.add_mesh(mesh)

        super(FilterPolarTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(FilterPolarTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterPolarTestHarness('statepoint.10.*', True)
    harness.main()
