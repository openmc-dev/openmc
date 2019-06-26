from math import sqrt

import openmc

from tests.testing_harness import PyAPITestHarness


class HexLatticeCoincidentTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        materials = openmc.Materials()

        fuel_mat = openmc.Material()
        fuel_mat.add_nuclide('U235', 4.9817E-03, 'ao')
        materials.append(fuel_mat)

        matrix = openmc.Material()
        matrix.set_density('atom/b-cm', 1.7742E-02)
        matrix.add_element('C', 1.0, 'ao')
        matrix.add_s_alpha_beta('c_Graphite')
        materials.append(matrix)

        lead = openmc.Material(name="Lead")
        lead.set_density('g/cm3', 10.32)
        lead.add_nuclide('Pb204', 0.014, 'ao')
        lead.add_nuclide('Pb206', 0.241, 'ao')
        lead.add_nuclide('Pb207', 0.221, 'ao')
        lead.add_nuclide('Pb208', 0.524, 'ao')
        materials.append(lead)

        coolant = openmc.Material()
        coolant.set_density('atom/b-cm', 5.4464E-04)
        coolant.add_nuclide('He4', 1.0, 'ao')
        materials.append(coolant)

        zirc = openmc.Material(name="Zirc4")
        zirc.add_nuclide('Zr90', 2.217E-02, 'ao')
        zirc.add_nuclide('Zr91', 4.781E-03, 'ao')
        zirc.add_nuclide('Zr92', 7.228E-03, 'ao')
        zirc.add_nuclide('Zr94', 7.169E-03, 'ao')
        zirc.add_nuclide('Zr96', 1.131E-03, 'ao')
        materials.append(zirc)

        materials.export_to_xml()

        ### Geometry ###
        pin_rad = 0.7 # cm
        assembly_pitch = 1.4 # cm

        cool_rad = 0.293 # cm
        zirc_clad_thickness = 0.057 # cm
        zirc_ir = cool_rad # cm
        zirc_or = cool_rad + zirc_clad_thickness # cm
        lead_thickness = 0.002 # cm
        lead_ir = zirc_or # cm
        lead_or = zirc_or + lead_thickness # cm

        cyl = openmc.ZCylinder(x0=0., y0=0., r=pin_rad)
        fuel_btm = openmc.ZPlane(z0=0.0, boundary_type = 'reflective')
        fuel_top = openmc.ZPlane(z0=10.0, boundary_type = 'reflective')
        region = -cyl & +fuel_btm & -fuel_top

        container = openmc.Cell(region=region)
        container.fill = fuel_mat

        fuel_outside = openmc.Cell()
        fuel_outside.region = +cyl
        fuel_outside.fill = matrix

        fuel_ch_univ = openmc.Universe(cells=[container, fuel_outside])

        # Coolant Channel
        cool_outer = openmc.ZCylinder(x0=0.0, y0=0.0, r=cool_rad)
        zirc_outer = openmc.ZCylinder(x0=0.0, y0=0.0, r=zirc_or)
        lead_outer = openmc.ZCylinder(x0=0.0, y0=0.0, r=lead_or)

        coolant_ch = openmc.Cell(name="coolant")
        coolant_ch.region = -cool_outer & +fuel_btm & -fuel_top
        coolant_ch.fill = coolant

        zirc_shell = openmc.Cell(name="zirconium_shell")
        zirc_shell.region = +cool_outer & -zirc_outer & +fuel_btm & -fuel_top
        zirc_shell.fill = zirc

        lead_shell = openmc.Cell(name="lead_shell")
        lead_shell.region = +zirc_outer & -lead_outer & +fuel_btm & -fuel_top
        lead_shell.fill = lead

        coolant_matrix = openmc.Cell(name="matrix coolant surround")
        coolant_matrix.region = +lead_outer & +fuel_btm & -fuel_top
        coolant_matrix.fill = matrix

        coolant_channel = [coolant_ch, zirc_shell, lead_shell, coolant_matrix]

        coolant_univ = openmc.Universe(name="coolant universe")
        coolant_univ.add_cells(coolant_channel)

        half_width = assembly_pitch # cm
        edge_length = (2./sqrt(3.0)) * half_width

        inf_mat = openmc.Cell()
        inf_mat.fill = matrix

        inf_mat_univ = openmc.Universe(cells=[inf_mat,])

        # a hex surface for the core to go inside of
        hexprism = openmc.model.hexagonal_prism(edge_length=edge_length,
                                                origin=(0.0, 0.0),
                                                boundary_type = 'reflective',
                                                orientation='x')

        pincell_only_lattice = openmc.HexLattice(name="regular fuel assembly")
        pincell_only_lattice.center = (0., 0.)
        pincell_only_lattice.pitch = (assembly_pitch,)
        pincell_only_lattice.outer = inf_mat_univ

        # setup hex rings
        ring0 = [fuel_ch_univ]
        ring1 = [coolant_univ] * 6
        pincell_only_lattice.universes = [ring1, ring0]

        pincell_only_cell = openmc.Cell(name="container cell")
        pincell_only_cell.region = hexprism & +fuel_btm & -fuel_top
        pincell_only_cell.fill = pincell_only_lattice

        root_univ = openmc.Universe(name="root universe", cells=[pincell_only_cell,])

        geom = openmc.Geometry(root_univ)
        geom.export_to_xml()

        ### Settings ###

        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'

        source = openmc.Source()
        corner_dist = sqrt(2) * pin_rad
        ll = [-corner_dist, -corner_dist, 0.0]
        ur = [corner_dist, corner_dist, 10.0]
        source.space = openmc.stats.Box(ll, ur)
        source.strength = 1.0
        settings.source = source
        settings.output = {'summary' : False}
        settings.batches = 5
        settings.inactive = 2
        settings.particles = 1000
        settings.seed = 22
        settings.export_to_xml()

def test_lattice_hex_coincident_surf():
    harness = HexLatticeCoincidentTestHarness('statepoint.5.h5')
    harness.main()
