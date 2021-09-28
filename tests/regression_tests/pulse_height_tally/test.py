import numpy as np
import openmc
import pandas as pd

from tests.testing_harness import PyAPITestHarness


class PulseHeightTallyTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        nai = openmc.Material(material_id=1, name='sodium iodide detector material')
        nai.set_density('g/cc', 3.67)
        nai.add_element('Na', 1.0)
        nai.add_element('I', 1.0)

        air = openmc.Material(material_id=2, name='air')
        air.set_density('g/cc', 0.001225)
        air.add_element('N', 4)
        air.add_nuclide('O16', 1)

        materials_file = openmc.Materials([nai, air])
        materials_file.export_to_xml()

        cylinder = openmc.ZCylinder(x0=0, y0=0, r=2.54, boundary_type='transmission',
            surface_id=1, name='cylinder')
        bottom = openmc.ZPlane(z0=0, boundary_type='transmission',
            surface_id=2, name='bottom')
        top = openmc.ZPlane(z0=5, boundary_type='transmission',
            surface_id=3, name='top')
        sphere = openmc.Sphere(x0=0, y0=0, z0=0, r=10, boundary_type='vacuum',
            surface_id=4, name='outer sphere')

        detector = openmc.Cell(cell_id=1, name='detector')
        detector.region = -cylinder & +bottom & -top
        detector.fill = nai

        surrounding = openmc.Cell(cell_id=2, name='surrounding')
        surrounding.region = ~(-cylinder & +bottom & -top) & -sphere
        surrounding.fill = air

        detector_cell = openmc.Universe(universe_id=1, name='detector cell')
        detector_cell.add_cells([detector, surrounding])

        # Instantiate root Cell and Universe
        root_cell = openmc.Cell(cell_id=3, name='root cell')
        root_cell.region = -sphere
        root_cell.fill = detector_cell
        root_univ = openmc.Universe(universe_id=2, name='root universe')
        root_univ.add_cell(root_cell)

        geometry = openmc.Geometry(root_univ)
        geometry.export_to_xml()

        # Instantiate a Settings object, set all runtime parameters
        settings_file = openmc.Settings()
        settings_file.seed = 1
        settings_file.batches = 1
        settings_file.inactive = 0
        settings_file.particles = 1000
        settings_file.photon_transport = True

        # create the source, particles start uniformly distributed
        # on a disc in z-directrion
        source = openmc.Source()
        source.space = openmc.stats.Point((0.0, 0.0, -1.0))
        source.particle = 'photon'
        source.angle = openmc.stats.Monodirectional(np.array([0.0, 0.0, 1.0]))
        source.energy = openmc.stats.Discrete([0.662e6], [1.0])
        settings_file.source = source

        settings_file.run_mode = 'fixed source'

        settings_file.export_to_xml()

        # Tallies file
        tallies_file = openmc.Tallies()

        cell_filter = openmc.CellFilter(detector)
        energy_filter = openmc.EnergyFilter(np.linspace(0, 1e6, 51))

        tally = openmc.Tally(name='pulse-height tally')
        tally.filters = [cell_filter, energy_filter]
        tally.scores = ['pulse-height']
        tallies_file.append(tally)

        tallies_file.export_to_xml()


    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Extract the tally data as a Pandas DataFrame.
        df = pd.DataFrame()
        for t in sp.tallies.values():
            df = df.append(t.get_pandas_dataframe(), ignore_index=True)

        # Extract the relevant data as a CSV string.
        return df['mean'].to_csv(None, columns=['mean'], index=False, float_format='%.3e')
        return outstr


def test_surface_tally():
    harness = PulseHeightTallyTestHarness('statepoint.1.h5')
    harness.main()
