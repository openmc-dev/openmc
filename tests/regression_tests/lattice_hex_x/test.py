from tests.testing_harness import PyAPITestHarness
import openmc
import numpy as np


class HexLatticeOXTestHarness(PyAPITestHarness):

    def _build_inputs(self):
        materials = openmc.Materials()

        fuel_mat = openmc.Material(material_id=1, name="UO2")
        fuel_mat.set_density('sum')
        fuel_mat.add_nuclide('U235', 0.87370e-03)
        fuel_mat.add_nuclide('U238', 1.87440e-02)
        fuel_mat.add_nuclide('O16', 3.92350e-02)
        materials.append(fuel_mat)

        coolant = openmc.Material(material_id=2, name="borated H2O")
        coolant.set_density('sum')
        coolant.add_nuclide('H1', 0.06694)
        coolant.add_nuclide('O16', 0.03347)
        coolant.add_nuclide('B10', 6.6262e-6)
        coolant.add_nuclide('B11', 2.6839e-5)
        materials.append(coolant)

        absorber = openmc.Material(material_id=3, name="pellet B4C")
        absorber.set_density('sum')
        absorber.add_nuclide('C0', 0.01966)
        absorber.add_nuclide('B11', 4.7344e-6)
        absorber.add_nuclide('B10', 1.9177e-5)
        materials.append(absorber)

        zirc = openmc.Material(material_id=4, name="Zirc4")
        zirc.set_density('sum')
        zirc.add_element('Zr', 4.23e-2)
        materials.append(zirc)

        materials.export_to_xml()

        # Geometry #

        pin_rad = 0.7 # cm
        assembly_pitch = 1.235 # cm
        hexagonal_pitch = 23.6 # cm
        length = 10.0 # cm

        # Fuel pin surfaces

        cylfuelin = openmc.ZCylinder(surface_id=1, r=0.386)
        cylfuelout = openmc.ZCylinder(surface_id=2, r=0.4582)

        # Fuel cells

        infcell = openmc.Cell(cell_id=1)
        infcell.region = -cylfuelin
        infcell.fill = fuel_mat

        clfcell = openmc.Cell(cell_id=2)
        clfcell.region = -cylfuelout & +cylfuelin
        clfcell.fill = zirc

        outfcell = openmc.Cell(cell_id=3)
        outfcell.region = +cylfuelout
        outfcell.fill = coolant

        # Fuel universe

        fuel_ch_univ = openmc.Universe(universe_id=1, name="Fuel channel",
                                       cells=[infcell, clfcell, outfcell])

        # Central tube surfaces

        cyltubein = openmc.ZCylinder(surface_id=3, r=0.45)
        cyltubeout = openmc.ZCylinder(surface_id=4, r=0.5177)

        # Central tube cells

        inctcell = openmc.Cell(cell_id=4)
        inctcell.region = -cyltubein
        inctcell.fill = coolant

        clctcell = openmc.Cell(cell_id=5)
        clctcell.region = -cyltubeout & +cyltubein
        clctcell.fill = zirc

        outctcell = openmc.Cell(cell_id=6)
        outctcell.region = +cyltubeout
        outctcell.fill = coolant

        # Central tubel universe

        tube_ch_univ = openmc.Universe(universe_id=2,
                                       name="Central tube channel",
                                       cells=[inctcell, clctcell, outctcell])

        # Absorber tube surfaces

        cylabsin = openmc.ZCylinder(surface_id=5, r=0.35)
        cylabsout = openmc.ZCylinder(surface_id=6, r=0.41)
        cylabsclin = openmc.ZCylinder(surface_id=7, r=0.545)
        cylabsclout = openmc.ZCylinder(surface_id=8, r=0.6323)

        # Absorber tube cells

        inabscell = openmc.Cell(cell_id=7)
        inabscell.region = -cylabsin
        inabscell.fill = absorber

        clabscell = openmc.Cell(cell_id=8)
        clabscell.region = -cylabsout & +cylabsin
        clabscell.fill = zirc

        interabscell = openmc.Cell(cell_id=9)
        interabscell.region = -cylabsclin & +cylabsout
        interabscell.fill = coolant

        clatcell = openmc.Cell(cell_id=10)
        clatcell.region = -cylabsclout & +cylabsclin
        clatcell.fill = zirc

        outabscell = openmc.Cell(cell_id=11)
        outabscell.region = +cylabsclout
        outabscell.fill = coolant

        # Absorber tube universe

        abs_ch_univ = openmc.Universe(universe_id=3,
                                      name="Central tube channel",
                                      cells=[inabscell, clabscell,
                                             interabscell,
                                             clatcell, outabscell])
        # Assembly surfaces

        edge_length = (1./np.sqrt(3.0)) * hexagonal_pitch
        fuel_bottom = openmc.ZPlane(surface_id=9, z0=0.0,
                                    boundary_type='reflective')
        fuel_top = openmc.ZPlane(surface_id=10, z0=length,
                                 boundary_type='reflective')

        # a hex surface for the core to go inside of

        hexprism = openmc.model.hexagonal_prism(edge_length=edge_length,
                                                origin=(0.0,  0.0),
                                                boundary_type='reflective',
                                                orientation='x')
        region = hexprism & +fuel_bottom & -fuel_top

        inf_mat = openmc.Cell(cell_id=12)
        inf_mat.fill = coolant
        inf_mat_univ = openmc.Universe(universe_id=4, cells=[inf_mat])

        # Fill lattice by channels

        nring = 11
        universes = []
        for ring in range(nring - 1, -1, -1):
            arr = []
            arr.append(fuel_ch_univ)
            for cell in range(ring * 6 - 1):
                arr.append(fuel_ch_univ)
            universes.append(arr)
        universes[-1] = [tube_ch_univ]
        channels = [(7, 2), (7, 5), (7, 8), (7, 11), (7, 14), (7, 17), (5, 0),
                    (4, 3), (5, 5), (4, 9), (5, 10), (4, 15), (5, 15),
                    (4, 21), (5, 20), (4, 27), (5, 25), (4, 33)]
        for i, j in channels:
            universes[i][j] = abs_ch_univ
        lattice = openmc.HexLattice(name="regular fuel assembly")
        lattice.orientation = "x"
        lattice.center = (0., 0., length/2.0)
        lattice.pitch = (assembly_pitch, length/2.0)
        lattice.universes = 2*[universes]
        lattice.outer = inf_mat_univ

        assembly_cell = openmc.Cell(cell_id=13,
                                    name="container assembly cell")
        assembly_cell.region = region
        assembly_cell.fill = lattice

        root_univ = openmc.Universe(universe_id=5, name="root universe",
                                    cells=[assembly_cell])

        geom = openmc.Geometry(root_univ)
        geom.export_to_xml()

        # Settings #

        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'

        source = openmc.Source()
        ll = [-edge_length, -edge_length, 0.0]
        ur = [edge_length, edge_length, 10.0]
        source.space = openmc.stats.Box(ll, ur)
        source.strength = 1.0
        settings.source = source
        settings.batches = 10
        settings.inactive = 5
        settings.particles = 1000
        settings.seed = 22
        settings.export_to_xml()


def test_lattice_hex_ox_surf():
    harness = HexLatticeOXTestHarness('statepoint.10.h5')
    harness.main()
