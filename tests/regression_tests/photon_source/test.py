from math import pi

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class SourceTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        mat = openmc.Material()
        mat.set_density('g/cm3', 0.998207)
        mat.add_element('H', 0.111894)
        mat.add_element('O', 0.888106)
        materials = openmc.Materials([mat])
        materials.export_to_xml()

        sphere = openmc.Sphere(R=1.0e9, boundary_type='reflective')
        inside_sphere = openmc.Cell()
        inside_sphere.region = -sphere
        inside_sphere.fill = mat
        root = openmc.Universe()
        root.add_cell(inside_sphere)
        geometry = openmc.Geometry(root)
        geometry.export_to_xml()
 
        source = openmc.Source()
        source.space = openmc.stats.Point((0, 0, 0))
        source.angle = openmc.stats.Isotropic()
        source.energy = openmc.stats.Discrete([10.0e6], [1.0])
        source.particle = 'photon'
 
        settings = openmc.Settings()
        settings.particles = 10000
        settings.batches = 1
        settings.photon_transport = True
        settings.electron_treatment = 'ttb'
        settings.cutoff = {'energy_photon' : 1000.0}
        settings.run_mode = 'fixed source'
        settings.source = source
        settings.export_to_xml()
 
        tally = openmc.Tally()
        tally.scores = ['flux']
        tallies = openmc.Tallies([tally])
        tallies.export_to_xml()

    def _get_results(self):
        sp = openmc.StatePoint(self._sp_name)
        outstr = ''
        t = sp.get_tally()
        outstr += 'tally {}:\n'.format(t.id)
        outstr += 'sum = {:12.6E}\n'.format(t.sum[0, 0, 0])
        outstr += 'sum_sq = {:12.6E}\n'.format(t.sum_sq[0, 0, 0])

        return outstr


def test_source():
    harness = SourceTestHarness('statepoint.1.h5')
    harness.main()
