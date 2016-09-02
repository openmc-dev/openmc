#!/usr/bin/env python
import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
import openmc
from openmc.stats import Box
from openmc.source import Source

class MultipoleTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        ####################
        # Materials
        ####################

        moderator = openmc.Material(material_id=1)
        moderator.set_density('g/cc', 1.0)
        moderator.add_nuclide('H1', 2.0)
        moderator.add_nuclide('O16', 1.0)

        dense_fuel = openmc.Material(material_id=2)
        dense_fuel.set_density('g/cc', 4.5)
        dense_fuel.add_nuclide('U235', 1.0)

        mats_file = openmc.Materials([moderator, dense_fuel])
        mats_file.export_to_xml()

        ####################
        # Geometry
        ####################

        c1 = openmc.Cell(cell_id=1, fill=moderator)
        mod_univ = openmc.Universe(universe_id=1, cells=(c1,))

        r0 = openmc.ZCylinder(R=0.3)
        c11 = openmc.Cell(cell_id=11, fill=dense_fuel, region=-r0)
        c11.temperature = [500, 0, 700, 800]
        c12 = openmc.Cell(cell_id=12, fill=moderator, region=+r0)
        fuel_univ = openmc.Universe(universe_id=11, cells=(c11, c12))

        lat = openmc.RectLattice(lattice_id=101)
        lat.dimension = [2, 2]
        lat.lower_left = [-2.0, -2.0]
        lat.pitch = [2.0, 2.0]
        lat.universes = [[fuel_univ]*2]*2
        lat.outer = mod_univ

        x0 = openmc.XPlane(x0=-3.0)
        x1 = openmc.XPlane(x0=3.0)
        y0 = openmc.YPlane(y0=-3.0)
        y1 = openmc.YPlane(y0=3.0)
        for s in [x0, x1, y0, y1]:
            s.boundary_type = 'reflective'
        c101 = openmc.Cell(cell_id=101, fill=lat, region=+x0 & -x1 & +y0 & -y1)
        root_univ = openmc.Universe(universe_id=0, cells=(c101,))

        geometry = openmc.Geometry(root_univ)
        geometry.export_to_xml()

        ####################
        # Settings
        ####################

        sets_file = openmc.Settings()
        sets_file.batches = 5
        sets_file.inactive = 0
        sets_file.particles = 1000
        sets_file.source = Source(space=Box([-1, -1, -1], [1, 1, 1]))
        sets_file.output = {'summary': True}
        sets_file.temperature = {'method': 'multipole'}
        sets_file.export_to_xml()

        ####################
        # Plots
        ####################

        plots_file = openmc.Plots()

        plot = openmc.Plot(plot_id=1)
        plot.basis = 'xy'
        plot.color = 'cell'
        plot.filename = 'cellplot'
        plot.origin = (0, 0, 0)
        plot.width = (7, 7)
        plot.pixels = (400, 400)
        plots_file.append(plot)

        plot = openmc.Plot(plot_id=2)
        plot.basis = 'xy'
        plot.color = 'mat'
        plot.filename = 'matplot'
        plot.origin = (0, 0, 0)
        plot.width = (7, 7)
        plot.pixels = (400, 400)
        plots_file.append(plot)

        plots_file.export_to_xml()

    def execute_test(self):
        if not 'OPENMC_MULTIPOLE_LIBRARY' in os.environ:
            raise RuntimeError("The 'OPENMC_MULTIPOLE_LIBRARY' environment "
                 "variable must be specified for this test.")
        else:
            super(MultipoleTestHarness, self).execute_test()

    def _get_results(self):
        outstr = super(MultipoleTestHarness, self)._get_results()
        su = openmc.Summary('summary.h5')
        outstr += str(su.get_cell_by_id(11))
        return outstr

    def _cleanup(self):
        f = os.path.join(os.getcwd(), 'plots.xml')
        if os.path.exists(f):
            os.remove(f)
        super(MultipoleTestHarness, self)._cleanup()


if __name__ == '__main__':
    harness = MultipoleTestHarness('statepoint.5.*')
    harness.main()
