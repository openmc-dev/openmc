import openmc

from tests.testing_harness import TestHarness, PyAPITestHarness


class DistribmatTestHarness(PyAPITestHarness):
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

        light_fuel = openmc.Material(material_id=3)
        light_fuel.set_density('g/cc', 2.0)
        light_fuel.add_nuclide('U235', 1.0)

        mats_file = openmc.Materials([moderator, dense_fuel, light_fuel])
        mats_file.export_to_xml()

        ####################
        # Geometry
        ####################

        c1 = openmc.Cell(cell_id=1, fill=moderator)
        mod_univ = openmc.Universe(universe_id=1, cells=[c1])

        r0 = openmc.ZCylinder(r=0.3)
        c11 = openmc.Cell(cell_id=11, region=-r0)
        c11.fill = [dense_fuel, None, light_fuel, dense_fuel]
        c12 = openmc.Cell(cell_id=12, region=+r0, fill=moderator)
        fuel_univ = openmc.Universe(universe_id=11, cells=[c11, c12])

        lat = openmc.RectLattice(lattice_id=101)
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
        c101 = openmc.Cell(cell_id=101, fill=lat)
        c101.region = +x0 & -x1 & +y0 & -y1
        root_univ = openmc.Universe(universe_id=0, cells=[c101])

        geometry = openmc.Geometry(root_univ)
        geometry.export_to_xml()

        ####################
        # Settings
        ####################

        sets_file = openmc.Settings()
        sets_file.batches = 5
        sets_file.inactive = 0
        sets_file.particles = 1000
        sets_file.source = openmc.Source(space=openmc.stats.Box(
            [-1, -1, -1], [1, 1, 1]))
        sets_file.export_to_xml()

        ####################
        # Plots
        ####################

        plot1 = openmc.Plot(plot_id=1)
        plot1.basis = 'xy'
        plot1.color_by = 'cell'
        plot1.filename = 'cellplot'
        plot1.origin = (0, 0, 0)
        plot1.width = (7, 7)
        plot1.pixels = (400, 400)

        plot2 = openmc.Plot(plot_id=2)
        plot2.basis = 'xy'
        plot2.color_by = 'material'
        plot2.filename = 'matplot'
        plot2.origin = (0, 0, 0)
        plot2.width = (7, 7)
        plot2.pixels = (400, 400)

        plots = openmc.Plots([plot1, plot2])
        plots.export_to_xml()

    def _get_results(self):
        outstr = super()._get_results()
        su = openmc.Summary('summary.h5')
        outstr += str(su.geometry.get_all_cells()[11])
        return outstr


def test_distribmat():
    harness = DistribmatTestHarness('statepoint.5.h5')
    harness.main()
