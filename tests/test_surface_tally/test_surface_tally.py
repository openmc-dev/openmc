#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import numpy as np
import openmc


class CreateSurfaceTallyTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Instantiate some Materials and register the appropriate Nuclides
        uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment')
        uo2.set_density('g/cm3', 10.29769)
        uo2.add_element('U', 1., enrichment=2.4)
        uo2.add_element('O', 2.)

        borated_water = openmc.Material(name='Borated water')
        borated_water.set_density('g/cm3', 0.740582)
        borated_water.add_element('B', 4.0e-5)
        borated_water.add_element('H', 5.0e-2)
        borated_water.add_element('O', 2.4e-2)
        borated_water.add_s_alpha_beta('c_H_in_H2O')

        # Instantiate a Materials collection and export to XML
        materials_file = openmc.Materials([uo2, borated_water])

        materials_file.export_to_xml()

        # Instantiate ZCylinder surfaces
        fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, R=0.4, name='Fuel OR')
        left    = openmc.XPlane(surface_id=2, x0=-0.62992, name='left')
        right   = openmc.XPlane(surface_id=3, x0=0.62992, name='right')
        bottom  = openmc.YPlane(y0=-0.62992, name='bottom')
        top     = openmc.YPlane(y0=0.62992, name='top')

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
        fuel.fill  = uo2
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

        # Instantiate a Geometry, register the root Universe, and export to XML
        geometry = openmc.Geometry(root_univ)
        geometry.export_to_xml()

        # Instantiate a Settings object, set all runtime parameters, and export to XML
        settings_file = openmc.Settings()
        settings_file.batches = 10
        settings_file.inactive = 5
        settings_file.particles = 1000

        # Create an initial uniform spatial source distribution over fissionable zones
        bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
        settings_file.source = openmc.source.Source(space=uniform_dist)
        settings_file.export_to_xml()

        # Tallies file
        tallies_file = openmc.Tallies()

        # Cell to cell tallies
        # These filters are same for all tallies
        two_groups = np.array([0., 4, 20.]) * 1e6
        energy_filter     = openmc.EnergyFilter(two_groups)
        polar_filter      = openmc.PolarFilter([0, np.pi / 4, np.pi])
        azimuthal_filter  = openmc.AzimuthalFilter([0, np.pi / 4, np.pi])


        cell_to_cell_tallies = []
        tally_index = 0 

        for cell1 in pin_cell.get_all_cells():
            for cell2 in pin_cell.get_all_cells():

                if cell1 != cell2 and abs(abs(cell1-cell2)-1) < 0.1:  # no need for cell1 to cell1

                    # Cell to cell filters for partial current
                    cell_from_filter = openmc.CellFromFilter(cell1)
                    cell_to_filter = openmc.CellFilter(cell2)

                    cell_to_cell_tallies.append(openmc.Tally(tally_id=tally_index, name=str(cell1)+'-'+str(cell2)))
                    cell_to_cell_tallies[tally_index].filters   = [cell_from_filter, cell_to_filter, energy_filter, polar_filter, azimuthal_filter]
                    cell_to_cell_tallies[tally_index].scores    = ['current']
                    tallies_file.append(cell_to_cell_tallies[tally_index])
                    tally_index += 1

                    # Cell from + surface filters for partial current
                    surface_filter = openmc.SurfaceFilter([1])
                    cell_to_cell_tallies.append(openmc.Tally(tally_id=tally_index, name=str(cell1)+'-surface1'))
                    cell_to_cell_tallies[tally_index].filters   = [cell_from_filter, surface_filter, energy_filter, polar_filter, azimuthal_filter]
                    cell_to_cell_tallies[tally_index].scores    = ['current']
                    tallies_file.append(cell_to_cell_tallies[tally_index])
                    tally_index += 1

        # Surface filter on inner surface, for net current
        surface_filter = openmc.SurfaceFilter([1])
        surf_tally1 = openmc.Tally(tally_id=tally_index, name='net_cylinder')
        surf_tally1.filters   = [surface_filter, energy_filter, polar_filter, azimuthal_filter]
        surf_tally1.scores    = ['current']
        tallies_file.append(surf_tally1)
        tally_index += 1

        # Surface filter on left surface, vacuum BC, for net current = leakage
        surface_filter = openmc.SurfaceFilter([2])
        surf_tally2 = openmc.Tally(tally_id=tally_index, name='leakage_left')
        surf_tally2.filters   = [surface_filter, energy_filter, polar_filter, azimuthal_filter]
        surf_tally2.scores    = ['current']
        tallies_file.append(surf_tally2)
        tally_index += 1

        # Surface filter on right surface, reflective, for net current = 0
        surface_filter = openmc.SurfaceFilter([3])
        surf_tally3 = openmc.Tally(tally_id=tally_index, name='net_right')
        surf_tally3.filters   = [surface_filter, energy_filter, polar_filter, azimuthal_filter]
        surf_tally3.scores    = ['current']
        tallies_file.append(surf_tally3)
        tally_index += 1

        tallies_file.export_to_xml()

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


if __name__ == '__main__':
    harness = CreateSurfaceTallyTestHarness('statepoint.10.h5')
    harness.main()
