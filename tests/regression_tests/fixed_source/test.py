import numpy as np

import openmc
import openmc.stats

from tests.testing_harness import PyAPITestHarness


class FixedSourceTestHarness(PyAPITestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        outstr = ''
        with openmc.StatePoint(self._sp_name) as sp:
            # Write out tally data.
            for i, tally_ind in enumerate(sp.tallies):
                tally = sp.tallies[tally_ind]
                results = np.zeros((tally.sum.size*2, ))
                results[0::2] = tally.sum.ravel()
                results[1::2] = tally.sum_sq.ravel()
                results = ['{0:12.6E}'.format(x) for x in results]

                outstr += 'tally ' + str(i + 1) + ':\n'
                outstr += '\n'.join(results) + '\n'

            gt = sp.global_tallies
            outstr += 'leakage:\n'
            outstr += '{0:12.6E}'.format(gt[gt['name'] == b'leakage'][0]['sum']) + '\n'
            outstr += '{0:12.6E}'.format(gt[gt['name'] == b'leakage'][0]['sum_sq']) + '\n'

        return outstr


def test_fixed_source():
    mat = openmc.Material()
    mat.add_nuclide('O16', 1.0)
    mat.add_nuclide('U238', 0.0001)
    mat.set_density('g/cc', 7.5)

    surf = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-surf)

    model = openmc.model.Model()
    model.geometry.root_universe = openmc.Universe(cells=[cell])
    model.materials.append(mat)

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 100
    model.settings.temperature = {'default': 294}
    model.settings.source = openmc.Source(space=openmc.stats.Point(),
                                          strength=10.0)

    tally = openmc.Tally()
    tally.scores = ['flux']
    model.tallies.append(tally)

    harness = FixedSourceTestHarness('statepoint.10.h5', model)
    harness.main()
