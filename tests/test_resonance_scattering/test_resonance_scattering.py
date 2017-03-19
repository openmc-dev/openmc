#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class ResonanceScatteringTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Materials
        mat = openmc.Material(material_id=1)
        mat.set_density('g/cc', 1.0)
        mat.add_nuclide('U238', 1.0)
        mat.add_nuclide('U235', 0.02)
        mat.add_nuclide('Pu239', 0.02)
        mat.add_nuclide('H1', 20.0)

        mats_file = openmc.Materials([mat])
        mats_file.export_to_xml()

        # Geometry
        dumb_surface = openmc.XPlane(x0=100, boundary_type='reflective')
        c1 = openmc.Cell(cell_id=1, fill=mat, region=-dumb_surface)
        root_univ = openmc.Universe(universe_id=0, cells=[c1])
        geometry = openmc.Geometry(root_univ)
        geometry.export_to_xml()

        # Settings
        res_scatt_dbrc = openmc.ResonanceScattering('U238', 'DBRC', 1.0, 210.0)
        res_scatt_wcm = openmc.ResonanceScattering('U235', 'WCM', 1.0, 210.0)
        res_scatt_ares = openmc.ResonanceScattering('Pu239', 'ARES', 1.0, 210.0)

        sets_file = openmc.Settings()
        sets_file.batches = 10
        sets_file.inactive = 5
        sets_file.particles = 1000
        sets_file.source = openmc.source.Source(
             space=openmc.stats.Box([-4, -4, -4], [4, 4, 4]))
        sets_file.resonance_scattering = [res_scatt_dbrc, res_scatt_wcm,
             res_scatt_ares]
        sets_file.export_to_xml()


if __name__ == '__main__':
    harness = ResonanceScatteringTestHarness('statepoint.10.h5')
    harness.main()
