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
        mat.add_nuclide('U-238', 1.0)
        mat.add_nuclide('U-235', 0.02)
        mat.add_nuclide('Pu-239', 0.02)
        mat.add_nuclide('H-1', 20.0)

        mats_file = openmc.MaterialsFile()
        mats_file.default_xs = '71c'
        mats_file.add_material(mat)
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
        geo_file = openmc.GeometryFile()
        geo_file.geometry = geometry
        geo_file.export_to_xml()

        # Settings
        nuclide = openmc.Nuclide('U-238', '71c')
        nuclide.zaid = 92238
        res_scatt_dbrc = openmc.ResonanceScattering()
        res_scatt_dbrc.nuclide = nuclide
        res_scatt_dbrc.nuclide_0K = nuclide # This is a bad idea! Just for tests
        res_scatt_dbrc.method = 'DBRC'
        res_scatt_dbrc.E_min = 1e-6
        res_scatt_dbrc.E_max = 210e-6

        nuclide = openmc.Nuclide('U-235', '71c')
        nuclide.zaid = 92235
        res_scatt_wcm = openmc.ResonanceScattering()
        res_scatt_wcm.nuclide = nuclide
        res_scatt_wcm.nuclide_0K = nuclide
        res_scatt_wcm.method = 'WCM'
        res_scatt_wcm.E_min = 1e-6
        res_scatt_wcm.E_max = 210e-6

        nuclide = openmc.Nuclide('Pu-239', '71c')
        nuclide.zaid = 94239
        res_scatt_ares = openmc.ResonanceScattering()
        res_scatt_ares.nuclide = nuclide
        res_scatt_ares.nuclide_0K = nuclide
        res_scatt_ares.method = 'ARES'
        res_scatt_ares.E_min = 1e-6
        res_scatt_ares.E_max = 210e-6

        sets_file = openmc.SettingsFile()
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
