#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc


class FilterMuTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        filt1 = openmc.Filter(type='mu',
                              bins=(-1.0, -0.5, 0.0, 0.5, 1.0))
        tally1 = openmc.Tally(tally_id=1)
        tally1.add_filter(filt1)
        tally1.add_score('scatter')
        tally1.add_score('nu-scatter')

        filt2 = openmc.Filter(type='mu', bins=(5,))
        tally2 = openmc.Tally(tally_id=2)
        tally2.add_filter(filt2)
        tally2.add_score('scatter')
        tally2.add_score('nu-scatter')

        mesh = openmc.Mesh(mesh_id=1)
        mesh.lower_left  = [-182.07, -182.07]
        mesh.upper_right = [182.07,  182.07]
        mesh.dimension = [2, 2]
        filt_mesh = openmc.Filter(type='mesh', bins=(1,))
        tally3 = openmc.Tally(tally_id=3)
        tally3.add_filter(filt2)
        tally3.add_filter(filt_mesh)
        tally3.add_score('scatter')
        tally3.add_score('nu-scatter')


        self._input_set.tallies = openmc.TalliesFile()
        self._input_set.tallies.add_tally(tally1)
        self._input_set.tallies.add_tally(tally2)
        self._input_set.tallies.add_tally(tally3)
        self._input_set.tallies.add_mesh(mesh)

        super(FilterMuTestHarness, self)._build_inputs()

    def _cleanup(self):
        super(FilterMuTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = FilterMuTestHarness('statepoint.10.*', True)
    harness.main()
