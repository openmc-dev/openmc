from math import pi

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class SourceTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        mat = openmc.Material()
        mat.set_density('g/cm3', 2.6989)
        mat.add_nuclide('Al27', 1.0)
        materials = openmc.Materials([mat])
        materials.export_to_xml()

        cyl = openmc.XCylinder(boundary_type='vacuum', R=1.0e-6)
        x_plane_left = openmc.XPlane(boundary_type='vacuum', x0=-1.0)
        x_plane_center = openmc.XPlane(boundary_type='transmission', x0=1.0)
        x_plane_right = openmc.XPlane(boundary_type='vacuum', x0=1.0e9)

        inner_cyl_left = openmc.Cell()
        inner_cyl_right = openmc.Cell()
        outer_cyl = openmc.Cell()

        inner_cyl_left.region = -cyl & +x_plane_left & -x_plane_center
        inner_cyl_right.region = -cyl & +x_plane_center & -x_plane_right
        outer_cyl.region = ~(-cyl & +x_plane_left & -x_plane_right)
        inner_cyl_right.fill = mat
        geometry = openmc.Geometry([inner_cyl_left, inner_cyl_right, outer_cyl])
        geometry.export_to_xml()

        source = openmc.Source()
        source.space = openmc.stats.Point((0,0,0))
        source.angle = openmc.stats.Monodirectional()
        source.energy = openmc.stats.Discrete([14.0], [1.0])
        source.particle = 'neutron'

        settings = openmc.Settings()
        settings.particles = 10000
        settings.run_mode = 'fixed source'
        settings.batches = 1
        settings.photon_transport = True
        settings.electron_treatment = 'ttb'
        settings.cutoff = {'energy_photon' : 1000.0}
        settings.source = source
        settings.export_to_xml()

        surface_filter = openmc.SurfaceFilter(cyl)
        particle_filter = openmc.ParticleFilter('photon')
        tally = openmc.Tally()
        tally.filters = [surface_filter, particle_filter]
        tally.scores = ['current']
        tallies = openmc.Tallies([tally])
        tallies.export_to_xml()

    def _get_results(self):
        with openmc.StatePoint(self._sp_name) as sp:
            outstr = ''
            t = sp.get_tally()
            outstr += 'tally {}:\n'.format(t.id)
            outstr += 'sum = {:12.6E}\n'.format(t.sum[0, 0, 0])
            outstr += 'sum_sq = {:12.6E}\n'.format(t.sum_sq[0, 0, 0])

            return outstr


def test_source():
    harness = SourceTestHarness('statepoint.1.h5')
    harness.main()
