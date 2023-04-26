import openmc

from tests.testing_harness import PyAPITestHarness


class TimeCutoffTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set time cutoff
        time_cutoff = 1e-7

        # Two spheres in each other
        s1 = openmc.Sphere(r=100, surface_id=1)
        s2 = openmc.Sphere(r=200, surface_id=2, boundary_type='vacuum')
        inner_sphere = openmc.Cell(cell_id=1, name='inner sphere')
        inner_sphere.region = -s1
        outer_sphere = openmc.Cell(cell_id=2, name='outer sphere')
        outer_sphere.region = +s1 & -s2
        root = openmc.Universe(universe_id=0, name='root universe')
        root.add_cell(inner_sphere)
        root.add_cell(outer_sphere)
        self._model.geometry = openmc.Geometry(root)

        # Set the running parameters
        settings_file = openmc.Settings()
        settings_file.run_mode = 'fixed source'
        settings_file.batches = 10
        settings_file.particles = 100
        settings_file.cutoff = {'time_neutron': time_cutoff}
        settings_file.source = openmc.source.Source(space=openmc.stats.Point(),
                                                    energy=openmc.stats.Discrete([1e4],[1]))
        self._model.settings = settings_file

        # Tally flux under time cutoff
        tallies = openmc.Tallies()
        tally = openmc.Tally(1)
        tally.scores = ['flux']
        cell_filter = openmc.CellFilter(outer_sphere)
        tally.filters = [cell_filter]
        tallies.append(tally)
        self._model.tallies = tallies

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Write out tally data.
        outstr = ''
        t = sp.get_tally()
        outstr += 'tally {}:\n'.format(t.id)
        outstr += 'sum = {:12.6E}\n'.format(t.sum[0, 0, 0])
        outstr += 'sum_sq = {:12.6E}\n'.format(t.sum_sq[0, 0, 0])

        return outstr


def test_time_cutoff():
    harness = TimeCutoffTestHarness('statepoint.10.h5', model=openmc.Model())
    harness.main()
