#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class ResonanceScatteringTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Nuclides
        u238 = openmc.Nuclide('U238')
        u235 = openmc.Nuclide('U235')
        pu239 = openmc.Nuclide('Pu239')
        h1 = openmc.Nuclide('H1')

        # Materials
        mat = openmc.Material(material_id=1)
        mat.set_density('g/cc', 1.0)
        mat.add_nuclide(u238, 1.0)
        mat.add_nuclide(u235, 0.02)
        mat.add_nuclide(pu239, 0.02)
        mat.add_nuclide(h1, 20.0)

        mats_file = openmc.Materials([mat])
        mats_file.export_to_xml()

        # Geometry
        dumb_surface = openmc.XPlane(x0=100)
        dumb_surface.boundary_type = 'reflective'

        c1 = openmc.Cell(cell_id=1)
        c1.fill = mat
        c1.region = -dumb_surface

        root_univ = openmc.Universe(universe_id=0)
        root_univ.add_cell(c1)

        geometry = openmc.Geometry()
        geometry.root_universe = root_univ
        geometry.export_to_xml()

        # Settings
        res_scatt_dbrc = openmc.ResonanceScattering(u238, 'DBRC', 1.0, 210.0)
        res_scatt_wcm = openmc.ResonanceScattering(u235, 'WCM', 1.0, 210.0)
        res_scatt_ares = openmc.ResonanceScattering(pu239, 'ARES', 1.0, 210.0)

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
    harness = ResonanceScatteringTestHarness('statepoint.10.*')
    harness.main()
