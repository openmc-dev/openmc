import numpy as np
import openmc
import pandas as pd

from tests.testing_harness import PyAPITestHarness


class SurfaceTallyTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Instantiate some Materials and register the appropriate Nuclides
        uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment')
        uo2.set_density('g/cc', 10.0)
        uo2.add_nuclide('U238', 1.0)
        uo2.add_nuclide('U235', 0.02)
        uo2.add_nuclide('O16', 2.0)

        borated_water = openmc.Material(name='Borated water')
        borated_water.set_density('g/cm3', 1)
        borated_water.add_nuclide('B10', 10e-5)
        borated_water.add_nuclide('H1', 2.0)
        borated_water.add_nuclide('O16', 1.0)

        # Instantiate a Materials collection and export to XML
        materials_file = openmc.Materials([uo2, borated_water])
        materials_file.export_to_xml()

        # Instantiate ZCylinder surfaces
        fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, r=1,
            name='Fuel OR')
        left = openmc.XPlane(surface_id=2, x0=-2, name='left')
        right = openmc.XPlane(surface_id=3, x0=2, name='right')
        bottom = openmc.YPlane(y0=-2, name='bottom')
        top = openmc.YPlane(y0=2, name='top')

        left.boundary_type = 'vacuum'
        right.boundary_type = 'reflective'
        top.boundary_type = 'reflective'
        bottom.boundary_type = 'reflective'

        # Instantiate Cells
        fuel = openmc.Cell(name='fuel')
        water = openmc.Cell(name='water')

        # Use surface half-spaces to define regions
        fuel.region = -fuel_or
        water.region = +fuel_or & -right & +bottom & -top

        # Register Materials with Cells
        fuel.fill = uo2
        water.fill = borated_water

        # Instantiate pin cell Universe
        pin_cell = openmc.Universe(name='pin cell')
        pin_cell.add_cells([fuel, water])

        # Instantiate root Cell and Universe
        root_cell = openmc.Cell(name='root cell')
        root_cell.region = +left & -right & +bottom & -top
        root_cell.fill = pin_cell
        root_univ = openmc.Universe(universe_id=0, name='root universe')
        root_univ.add_cell(root_cell)

        # Instantiate a Geometry, register the root Universe
        geometry = openmc.Geometry(root_univ)
        geometry.export_to_xml()

        # Instantiate a Settings object, set all runtime parameters
        settings_file = openmc.Settings()
        settings_file.batches = 10
        settings_file.inactive = 0
        settings_file.particles = 1000
        #settings_file.output = {'tallies': True}

        # Create an initial uniform spatial source distribution
        bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],\
             only_fissionable=True)
        settings_file.source = openmc.source.Source(space=uniform_dist)
        settings_file.export_to_xml()

        # Tallies file
        tallies_file = openmc.Tallies()

        # Create partial current tallies from fuel to water
        # Filters
        two_groups = [0., 4e6, 20e6]
        energy_filter = openmc.EnergyFilter(two_groups)
        polar_filter = openmc.PolarFilter([0, np.pi / 4, np.pi])
        azimuthal_filter = openmc.AzimuthalFilter([0, np.pi / 4, np.pi])
        surface_filter = openmc.SurfaceFilter([1])
        cell_from_filter = openmc.CellFromFilter(fuel)
        cell_filter = openmc.CellFilter(water)

        # Use Cell to cell filters for partial current
        cell_to_cell_tally = openmc.Tally(name=str('fuel_to_water_1'))
        cell_to_cell_tally.filters = [cell_from_filter, cell_filter, \
             energy_filter, polar_filter, azimuthal_filter]
        cell_to_cell_tally.scores = ['current']
        tallies_file.append(cell_to_cell_tally)

        # Use a Cell from + surface filters for partial current
        cell_to_cell_tally = openmc.Tally(name=str('fuel_to_water_2'))
        cell_to_cell_tally.filters = [cell_from_filter, surface_filter, \
             energy_filter, polar_filter, azimuthal_filter]
        cell_to_cell_tally.scores = ['current']
        tallies_file.append(cell_to_cell_tally)

        # Create partial current tallies from water to fuel
        # Filters
        cell_from_filter = openmc.CellFromFilter(water)
        cell_filter = openmc.CellFilter(fuel)

        # Cell to cell filters for partial current
        cell_to_cell_tally = openmc.Tally(name=str('water_to_fuel_1'))
        cell_to_cell_tally.filters = [cell_from_filter, cell_filter, \
             energy_filter, polar_filter, azimuthal_filter]
        cell_to_cell_tally.scores = ['current']
        tallies_file.append(cell_to_cell_tally)

        # Cell from + surface filters for partial current
        cell_to_cell_tally = openmc.Tally(name=str('water_to_fuel_2'))
        cell_to_cell_tally.filters = [cell_from_filter, surface_filter, \
             energy_filter, polar_filter, azimuthal_filter]
        cell_to_cell_tally.scores = ['current']
        tallies_file.append(cell_to_cell_tally)

        # Create a net current tally on inner surface using a surface filter
        surface_filter = openmc.SurfaceFilter([1])
        surf_tally1 = openmc.Tally(name='net_cylinder')
        surf_tally1.filters = [surface_filter, energy_filter, polar_filter, \
             azimuthal_filter]
        surf_tally1.scores = ['current']
        tallies_file.append(surf_tally1)

        # Create a net current tally on left surface using a surface filter
        # This surface has a vacuum boundary condition, so leakage is tallied
        surface_filter = openmc.SurfaceFilter([2])
        surf_tally2 = openmc.Tally(name='leakage_left')
        surf_tally2.filters = [surface_filter, energy_filter, polar_filter, \
            azimuthal_filter]
        surf_tally2.scores = ['current']
        tallies_file.append(surf_tally2)

        # Create a net current tally on right surface using a surface filter
        # This surface has a reflective boundary condition, so the net current
        # should be zero.
        surface_filter = openmc.SurfaceFilter([3])
        surf_tally3 = openmc.Tally(name='net_right')
        surf_tally3.filters = [surface_filter, energy_filter]
        surf_tally3.scores = ['current']
        tallies_file.append(surf_tally3)

        surface_filter = openmc.SurfaceFilter([3])
        surf_tally3 = openmc.Tally(name='net_right')
        surf_tally3.filters = [surface_filter, energy_filter]
        surf_tally3.scores = ['current']
        tallies_file.append(surf_tally3)

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
        cols = ('mean', 'std. dev.')
        return df.to_csv(None, columns=cols, index=False, float_format='%.7e')
        return outstr


def test_surface_tally():
    harness = SurfaceTallyTestHarness('statepoint.10.h5')
    harness.main()
