import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class TimeCutoffTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        NaI = openmc.Material(material_id=1)
        NaI.set_density('g/cc', 3.7)
        NaI.add_element('Na', 1.0)
        NaI.add_element('I', 1.0)

        self._model.materials = openmc.Materials([NaI])

        # Two spheres in each other
        s1 = openmc.Sphere(r=1, surface_id=1)
        s2 = openmc.Sphere(r=2, surface_id=2, boundary_type='vacuum')
        inner_sphere = openmc.Cell(cell_id=1, name='inner sphere')
        inner_sphere.region = -s1
        inner_sphere.fill = NaI
        outer_sphere = openmc.Cell(cell_id=2, name='outer sphere')
        outer_sphere.region = +s1 & -s2
        root = openmc.Universe(universe_id=0, name='root universe')
        root.add_cell(inner_sphere)
        root.add_cell(outer_sphere)
        self._model.geometry = openmc.Geometry(root)

        # Set the running parameters
        settings_file = openmc.Settings()
        settings_file.run_mode = 'fixed source'
        settings_file.batches = 5
        settings_file.particles = 100
        settings_file.photon_transport = True
        settings_file.source = openmc.source.Source(space=openmc.stats.Point(),
                                                    energy=openmc.stats.Discrete([1e6],[1]),
                                                    particle='photon')
        self._model.settings = settings_file

        # Tally flux under time cutoff
        tallies = openmc.Tallies()
        tally = openmc.Tally(name="pht tally")
        tally.scores = ['pulse-height']
        cell_filter = openmc.CellFilter(inner_sphere)
        energy_filter = openmc.EnergyFilter(np.linspace(0, 1_000_000, 101))
        tally.filters = [cell_filter, energy_filter]
        tallies.append(tally)
        self._model.tallies = tallies

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        
        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)
        tally = sp.get_tally(name='pht tally')
        
        # Write out tally data.
        outstr = ''
        outstr += '\n'.join(map('{:.8e}'.format, tally.mean.flatten())) + '\n'
        outstr += '\n'.join(map('{:.8e}'.format, tally.std_dev.flatten())) + '\n'
        
        return outstr


def test_time_cutoff():
    harness = TimeCutoffTestHarness('statepoint.5.h5', model=openmc.Model())
    harness.main()