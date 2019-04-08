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

        cyl = openmc.XCylinder(boundary_type='vacuum', r=1.0)
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
        source.energy = openmc.stats.Discrete([14.0e6], [1.0])
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
        current_tally = openmc.Tally()
        current_tally.filters = [surface_filter, particle_filter]
        current_tally.scores = ['current']
        total_tally_tracklength = openmc.Tally()
        total_tally_tracklength.filters = [particle_filter]
        total_tally_tracklength.scores = ['total']
        total_tally_tracklength.estimator = 'tracklength'
        total_tally_collision = openmc.Tally()
        total_tally_collision.filters = [particle_filter]
        total_tally_collision.scores = ['total']
        total_tally_collision.estimator = 'collision'
        total_tally_analog = openmc.Tally()
        total_tally_analog.filters = [particle_filter]
        total_tally_analog.scores = ['total']
        total_tally_analog.estimator = 'analog'
        tallies = openmc.Tallies([current_tally, total_tally_tracklength,
                                  total_tally_collision, total_tally_analog])
        tallies.export_to_xml()

    def _get_results(self):
        with openmc.StatePoint(self._sp_name) as sp:
            outstr = ''
            for i, tally_ind in enumerate(sp.tallies):
                tally = sp.tallies[tally_ind]
                results = np.zeros((tally.sum.size * 2, ))
                results[0::2] = tally.sum.ravel()
                results[1::2] = tally.sum_sq.ravel()
                results = ['{0:12.6E}'.format(x) for x in results]

                outstr += 'tally {}:\n'.format(i + 1)
                outstr += '\n'.join(results) + '\n'
            return outstr


def test_photon_production():
    harness = SourceTestHarness('statepoint.1.h5')
    harness.main()
