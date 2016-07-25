import openmc
from openmc.source import Source
from openmc.stats import Box


class InputSet(object):
    def __init__(self):
        self.settings = openmc.Settings()
        self.materials = openmc.Materials()
        self.geometry = openmc.Geometry()
        self.tallies = None
        self.plots = None

    def export(self):
        self.settings.export_to_xml()
        self.materials.export_to_xml()
        self.geometry.export_to_xml()
        if self.tallies is not None:
            self.tallies.export_to_xml()
        if self.plots is not None:
            self.plots.export_to_xml()

    def build_default_materials_and_geometry(self):
        # Define materials.
        fuel = openmc.Material(name='Fuel', material_id=1)
        fuel.set_density('g/cm3', 10.062)
        fuel.add_nuclide("U-234", 4.9476e-6)
        fuel.add_nuclide("U-235", 4.8218e-4)
        fuel.add_nuclide("U-236", 9.0402e-5)
        fuel.add_nuclide("U-238", 2.1504e-2)
        fuel.add_nuclide("Np-237", 7.3733e-6)
        fuel.add_nuclide("Pu-238", 1.5148e-6)
        fuel.add_nuclide("Pu-239", 1.3955e-4)
        fuel.add_nuclide("Pu-240", 3.4405e-5)
        fuel.add_nuclide("Pu-241", 2.1439e-5)
        fuel.add_nuclide("Pu-242", 3.7422e-6)
        fuel.add_nuclide("Am-241", 4.5041e-7)
        fuel.add_nuclide("Am-242m", 9.2301e-9)
        fuel.add_nuclide("Am-243", 4.7878e-7)
        fuel.add_nuclide("Cm-242", 1.0485e-7)
        fuel.add_nuclide("Cm-243", 1.4268e-9)
        fuel.add_nuclide("Cm-244", 8.8756e-8)
        fuel.add_nuclide("Cm-245", 3.5285e-9)
        fuel.add_nuclide("Mo-95", 2.6497e-5)
        fuel.add_nuclide("Tc-99", 3.2772e-5)
        fuel.add_nuclide("Ru-101", 3.0742e-5)
        fuel.add_nuclide("Ru-103", 2.3505e-6)
        fuel.add_nuclide("Ag-109", 2.0009e-6)
        fuel.add_nuclide("Xe-135", 1.0801e-8)
        fuel.add_nuclide("Cs-133", 3.4612e-5)
        fuel.add_nuclide("Nd-143", 2.6078e-5)
        fuel.add_nuclide("Nd-145", 1.9898e-5)
        fuel.add_nuclide("Sm-147", 1.6128e-6)
        fuel.add_nuclide("Sm-149", 1.1627e-7)
        fuel.add_nuclide("Sm-150", 7.1727e-6)
        fuel.add_nuclide("Sm-151", 5.4947e-7)
        fuel.add_nuclide("Sm-152", 3.0221e-6)
        fuel.add_nuclide("Eu-153", 2.6209e-6)
        fuel.add_nuclide("Gd-155", 1.5369e-9)
        fuel.add_nuclide("O-16", 4.5737e-2)

        clad = openmc.Material(name='Cladding', material_id=2)
        clad.set_density('g/cm3', 5.77)
        clad.add_nuclide("Zr-90", 0.5145)
        clad.add_nuclide("Zr-91", 0.1122)
        clad.add_nuclide("Zr-92", 0.1715)
        clad.add_nuclide("Zr-94", 0.1738)
        clad.add_nuclide("Zr-96", 0.0280)

        cold_water = openmc.Material(name='Cold borated water', material_id=3)
        cold_water.set_density('atom/b-cm', 0.07416)
        cold_water.add_nuclide("H-1", 2.0)
        cold_water.add_nuclide("O-16", 1.0)
        cold_water.add_nuclide("B-10", 6.490e-4)
        cold_water.add_nuclide("B-11", 2.689e-3)
        cold_water.add_s_alpha_beta('HH2O', '71t')

        hot_water = openmc.Material(name='Hot borated water', material_id=4)
        hot_water.set_density('atom/b-cm', 0.06614)
        hot_water.add_nuclide("H-1", 2.0)
        hot_water.add_nuclide("O-16", 1.0)
        hot_water.add_nuclide("B-10", 6.490e-4)
        hot_water.add_nuclide("B-11", 2.689e-3)
        hot_water.add_s_alpha_beta('HH2O', '71t')

        rpv_steel = openmc.Material(name='Reactor pressure vessel steel',
                                    material_id=5)
        rpv_steel.set_density('g/cm3', 7.9)
        rpv_steel.add_nuclide("Fe-54", 0.05437098, 'wo')
        rpv_steel.add_nuclide("Fe-56", 0.88500663, 'wo')
        rpv_steel.add_nuclide("Fe-57", 0.0208008, 'wo')
        rpv_steel.add_nuclide("Fe-58", 0.00282159, 'wo')
        rpv_steel.add_nuclide("Ni-58", 0.0067198, 'wo')
        rpv_steel.add_nuclide("Ni-60", 0.0026776, 'wo')
        rpv_steel.add_nuclide("Ni-61", 0.0001183, 'wo')
        rpv_steel.add_nuclide("Ni-62", 0.0003835, 'wo')
        rpv_steel.add_nuclide("Ni-64", 0.0001008, 'wo')
        rpv_steel.add_nuclide("Mn-55", 0.01, 'wo')
        rpv_steel.add_nuclide("Mo-92", 0.000849, 'wo')
        rpv_steel.add_nuclide("Mo-94", 0.0005418, 'wo')
        rpv_steel.add_nuclide("Mo-95", 0.0009438, 'wo')
        rpv_steel.add_nuclide("Mo-96", 0.0010002, 'wo')
        rpv_steel.add_nuclide("Mo-97", 0.0005796, 'wo')
        rpv_steel.add_nuclide("Mo-98", 0.0014814, 'wo')
        rpv_steel.add_nuclide("Mo-100", 0.0006042, 'wo')
        rpv_steel.add_nuclide("Si-28", 0.00367464, 'wo')
        rpv_steel.add_nuclide("Si-29", 0.00019336, 'wo')
        rpv_steel.add_nuclide("Si-30", 0.000132, 'wo')
        rpv_steel.add_nuclide("Cr-50", 0.00010435, 'wo')
        rpv_steel.add_nuclide("Cr-52", 0.002092475, 'wo')
        rpv_steel.add_nuclide("Cr-53", 0.00024185, 'wo')
        rpv_steel.add_nuclide("Cr-54", 6.1325e-05, 'wo')
        rpv_steel.add_nuclide("C-Nat", 0.0025, 'wo')
        rpv_steel.add_nuclide("Cu-63", 0.0013696, 'wo')
        rpv_steel.add_nuclide("Cu-65", 0.0006304, 'wo')

        lower_rad_ref = openmc.Material(name='Lower radial reflector',
                                        material_id=6)
        lower_rad_ref.set_density('g/cm3', 4.32)
        lower_rad_ref.add_nuclide("H-1", 0.0095661, 'wo')
        lower_rad_ref.add_nuclide("O-16", 0.0759107, 'wo')
        lower_rad_ref.add_nuclide("B-10", 3.08409e-5, 'wo')
        lower_rad_ref.add_nuclide("B-11", 1.40499e-4, 'wo')
        lower_rad_ref.add_nuclide("Fe-54", 0.035620772088, 'wo')
        lower_rad_ref.add_nuclide("Fe-56", 0.579805982228, 'wo')
        lower_rad_ref.add_nuclide("Fe-57", 0.01362750048, 'wo')
        lower_rad_ref.add_nuclide("Fe-58", 0.001848545204, 'wo')
        lower_rad_ref.add_nuclide("Ni-58", 0.055298376566, 'wo')
        lower_rad_ref.add_nuclide("Ni-60", 0.022034425592, 'wo')
        lower_rad_ref.add_nuclide("Ni-61", 0.000973510811, 'wo')
        lower_rad_ref.add_nuclide("Ni-62", 0.003155886695, 'wo')
        lower_rad_ref.add_nuclide("Ni-64", 0.000829500336, 'wo')
        lower_rad_ref.add_nuclide("Mn-55", 0.0182870, 'wo')
        lower_rad_ref.add_nuclide("Si-28", 0.00839976771, 'wo')
        lower_rad_ref.add_nuclide("Si-29", 0.00044199679, 'wo')
        lower_rad_ref.add_nuclide("Si-30", 0.0003017355, 'wo')
        lower_rad_ref.add_nuclide("Cr-50", 0.007251360806, 'wo')
        lower_rad_ref.add_nuclide("Cr-52", 0.145407678031, 'wo')
        lower_rad_ref.add_nuclide("Cr-53", 0.016806340306, 'wo')
        lower_rad_ref.add_nuclide("Cr-54", 0.004261520857, 'wo')
        lower_rad_ref.add_s_alpha_beta('HH2O', '71t')

        upper_rad_ref = openmc.Material(name='Upper radial reflector /'
                                             'Top plate region', material_id=7)
        upper_rad_ref.set_density('g/cm3', 4.28)
        upper_rad_ref.add_nuclide("H-1", 0.0086117, 'wo')
        upper_rad_ref.add_nuclide("O-16", 0.0683369, 'wo')
        upper_rad_ref.add_nuclide("B-10", 2.77638e-5, 'wo')
        upper_rad_ref.add_nuclide("B-11", 1.26481e-4, 'wo')
        upper_rad_ref.add_nuclide("Fe-54", 0.035953677186, 'wo')
        upper_rad_ref.add_nuclide("Fe-56", 0.585224740891, 'wo')
        upper_rad_ref.add_nuclide("Fe-57", 0.01375486056, 'wo')
        upper_rad_ref.add_nuclide("Fe-58", 0.001865821363, 'wo')
        upper_rad_ref.add_nuclide("Ni-58", 0.055815129186, 'wo')
        upper_rad_ref.add_nuclide("Ni-60", 0.022240333032, 'wo')
        upper_rad_ref.add_nuclide("Ni-61", 0.000982608081, 'wo')
        upper_rad_ref.add_nuclide("Ni-62", 0.003185377845, 'wo')
        upper_rad_ref.add_nuclide("Ni-64", 0.000837251856, 'wo')
        upper_rad_ref.add_nuclide("Mn-55", 0.0184579, 'wo')
        upper_rad_ref.add_nuclide("Si-28", 0.00847831314, 'wo')
        upper_rad_ref.add_nuclide("Si-29", 0.00044612986, 'wo')
        upper_rad_ref.add_nuclide("Si-30", 0.000304557, 'wo')
        upper_rad_ref.add_nuclide("Cr-50", 0.00731912987, 'wo')
        upper_rad_ref.add_nuclide("Cr-52", 0.146766614995, 'wo')
        upper_rad_ref.add_nuclide("Cr-53", 0.01696340737, 'wo')
        upper_rad_ref.add_nuclide("Cr-54", 0.004301347765, 'wo')
        upper_rad_ref.add_s_alpha_beta('HH2O', '71t')

        bot_plate = openmc.Material(name='Bottom plate region', material_id=8)
        bot_plate.set_density('g/cm3', 7.184)
        bot_plate.add_nuclide("H-1", 0.0011505, 'wo')
        bot_plate.add_nuclide("O-16", 0.0091296, 'wo')
        bot_plate.add_nuclide("B-10", 3.70915e-6, 'wo')
        bot_plate.add_nuclide("B-11", 1.68974e-5, 'wo')
        bot_plate.add_nuclide("Fe-54", 0.03855611055, 'wo')
        bot_plate.add_nuclide("Fe-56", 0.627585036425, 'wo')
        bot_plate.add_nuclide("Fe-57", 0.014750478, 'wo')
        bot_plate.add_nuclide("Fe-58", 0.002000875025, 'wo')
        bot_plate.add_nuclide("Ni-58", 0.059855207342, 'wo')
        bot_plate.add_nuclide("Ni-60", 0.023850159704, 'wo')
        bot_plate.add_nuclide("Ni-61", 0.001053732407, 'wo')
        bot_plate.add_nuclide("Ni-62", 0.003415945715, 'wo')
        bot_plate.add_nuclide("Ni-64", 0.000897854832, 'wo')
        bot_plate.add_nuclide("Mn-55", 0.0197940, 'wo')
        bot_plate.add_nuclide("Si-28", 0.00909197802, 'wo')
        bot_plate.add_nuclide("Si-29", 0.00047842098, 'wo')
        bot_plate.add_nuclide("Si-30", 0.000326601, 'wo')
        bot_plate.add_nuclide("Cr-50", 0.007848910646, 'wo')
        bot_plate.add_nuclide("Cr-52", 0.157390026871, 'wo')
        bot_plate.add_nuclide("Cr-53", 0.018191270146, 'wo')
        bot_plate.add_nuclide("Cr-54", 0.004612692337, 'wo')
        bot_plate.add_s_alpha_beta('HH2O', '71t')

        bot_nozzle = openmc.Material(name='Bottom nozzle region',
                                     material_id=9)
        bot_nozzle.set_density('g/cm3', 2.53)
        bot_nozzle.add_nuclide("H-1", 0.0245014, 'wo')
        bot_nozzle.add_nuclide("O-16", 0.1944274, 'wo')
        bot_nozzle.add_nuclide("B-10", 7.89917e-5, 'wo')
        bot_nozzle.add_nuclide("B-11", 3.59854e-4, 'wo')
        bot_nozzle.add_nuclide("Fe-54", 0.030411411144, 'wo')
        bot_nozzle.add_nuclide("Fe-56", 0.495012237964, 'wo')
        bot_nozzle.add_nuclide("Fe-57", 0.01163454624, 'wo')
        bot_nozzle.add_nuclide("Fe-58", 0.001578204652, 'wo')
        bot_nozzle.add_nuclide("Ni-58", 0.047211231662, 'wo')
        bot_nozzle.add_nuclide("Ni-60", 0.018811987544, 'wo')
        bot_nozzle.add_nuclide("Ni-61", 0.000831139127, 'wo')
        bot_nozzle.add_nuclide("Ni-62", 0.002694352115, 'wo')
        bot_nozzle.add_nuclide("Ni-64", 0.000708189552, 'wo')
        bot_nozzle.add_nuclide("Mn-55", 0.0156126, 'wo')
        bot_nozzle.add_nuclide("Si-28", 0.007171335558, 'wo')
        bot_nozzle.add_nuclide("Si-29", 0.000377356542, 'wo')
        bot_nozzle.add_nuclide("Si-30", 0.0002576079, 'wo')
        bot_nozzle.add_nuclide("Cr-50", 0.006190885148, 'wo')
        bot_nozzle.add_nuclide("Cr-52", 0.124142524198, 'wo')
        bot_nozzle.add_nuclide("Cr-53", 0.014348496148, 'wo')
        bot_nozzle.add_nuclide("Cr-54", 0.003638294506, 'wo')
        bot_nozzle.add_s_alpha_beta('HH2O', '71t')

        top_nozzle = openmc.Material(name='Top nozzle region', material_id=10)
        top_nozzle.set_density('g/cm3', 1.746)
        top_nozzle.add_nuclide("H-1", 0.0358870, 'wo')
        top_nozzle.add_nuclide("O-16", 0.2847761, 'wo')
        top_nozzle.add_nuclide("B-10", 1.15699e-4, 'wo')
        top_nozzle.add_nuclide("B-11", 5.27075e-4, 'wo')
        top_nozzle.add_nuclide("Fe-54", 0.02644016154, 'wo')
        top_nozzle.add_nuclide("Fe-56", 0.43037146399, 'wo')
        top_nozzle.add_nuclide("Fe-57", 0.0101152584, 'wo')
        top_nozzle.add_nuclide("Fe-58", 0.00137211607, 'wo')
        top_nozzle.add_nuclide("Ni-58", 0.04104621835, 'wo')
        top_nozzle.add_nuclide("Ni-60", 0.0163554502, 'wo')
        top_nozzle.add_nuclide("Ni-61", 0.000722605975, 'wo')
        top_nozzle.add_nuclide("Ni-62", 0.002342513875, 'wo')
        top_nozzle.add_nuclide("Ni-64", 0.0006157116, 'wo')
        top_nozzle.add_nuclide("Mn-55", 0.0135739, 'wo')
        top_nozzle.add_nuclide("Si-28", 0.006234853554, 'wo')
        top_nozzle.add_nuclide("Si-29", 0.000328078746, 'wo')
        top_nozzle.add_nuclide("Si-30", 0.0002239677, 'wo')
        top_nozzle.add_nuclide("Cr-50", 0.005382452306, 'wo')
        top_nozzle.add_nuclide("Cr-52", 0.107931450781, 'wo')
        top_nozzle.add_nuclide("Cr-53", 0.012474806806, 'wo')
        top_nozzle.add_nuclide("Cr-54", 0.003163190107, 'wo')
        top_nozzle.add_s_alpha_beta('HH2O', '71t')

        top_fa = openmc.Material(name='Top of fuel assemblies', material_id=11)
        top_fa.set_density('g/cm3', 3.044)
        top_fa.add_nuclide("H-1", 0.0162913, 'wo')
        top_fa.add_nuclide("O-16", 0.1292776, 'wo')
        top_fa.add_nuclide("B-10", 5.25228e-5, 'wo')
        top_fa.add_nuclide("B-11", 2.39272e-4, 'wo')
        top_fa.add_nuclide("Zr-90", 0.43313403903, 'wo')
        top_fa.add_nuclide("Zr-91", 0.09549277374, 'wo')
        top_fa.add_nuclide("Zr-92", 0.14759527104, 'wo')
        top_fa.add_nuclide("Zr-94", 0.15280552077, 'wo')
        top_fa.add_nuclide("Zr-96", 0.02511169542, 'wo')
        top_fa.add_s_alpha_beta('HH2O', '71t')

        bot_fa = openmc.Material(name='Bottom of fuel assemblies',
                                 material_id=12)
        bot_fa.set_density('g/cm3', 1.762)
        bot_fa.add_nuclide("H-1", 0.0292856, 'wo')
        bot_fa.add_nuclide("O-16", 0.2323919, 'wo')
        bot_fa.add_nuclide("B-10", 9.44159e-5, 'wo')
        bot_fa.add_nuclide("B-11", 4.30120e-4, 'wo')
        bot_fa.add_nuclide("Zr-90", 0.3741373658, 'wo')
        bot_fa.add_nuclide("Zr-91", 0.0824858164, 'wo')
        bot_fa.add_nuclide("Zr-92", 0.1274914944, 'wo')
        bot_fa.add_nuclide("Zr-94", 0.1319920622, 'wo')
        bot_fa.add_nuclide("Zr-96", 0.0216912612, 'wo')
        bot_fa.add_s_alpha_beta('HH2O', '71t')

        # Define the materials file.
        self.materials.default_xs = '71c'
        self.materials += (fuel, clad, cold_water, hot_water, rpv_steel,
                           lower_rad_ref, upper_rad_ref, bot_plate,
                           bot_nozzle, top_nozzle, top_fa, bot_fa)

        # Define surfaces.
        s1 = openmc.ZCylinder(R=0.41, surface_id=1)
        s2 = openmc.ZCylinder(R=0.475, surface_id=2)
        s3 = openmc.ZCylinder(R=0.56, surface_id=3)
        s4 = openmc.ZCylinder(R=0.62, surface_id=4)
        s5 = openmc.ZCylinder(R=187.6, surface_id=5)
        s6 = openmc.ZCylinder(R=209.0, surface_id=6)
        s7 = openmc.ZCylinder(R=229.0, surface_id=7)
        s8 = openmc.ZCylinder(R=249.0, surface_id=8)
        s8.boundary_type = 'vacuum'

        s31 = openmc.ZPlane(z0=-229.0, surface_id=31)
        s31.boundary_type = 'vacuum'
        s32 = openmc.ZPlane(z0=-199.0, surface_id=32)
        s33 = openmc.ZPlane(z0=-193.0, surface_id=33)
        s34 = openmc.ZPlane(z0=-183.0, surface_id=34)
        s35 = openmc.ZPlane(z0=0.0, surface_id=35)
        s36 = openmc.ZPlane(z0=183.0, surface_id=36)
        s37 = openmc.ZPlane(z0=203.0, surface_id=37)
        s38 = openmc.ZPlane(z0=215.0, surface_id=38)
        s39 = openmc.ZPlane(z0=223.0, surface_id=39)
        s39.boundary_type = 'vacuum'

        # Define pin cells.
        fuel_cold = openmc.Universe(name='Fuel pin, cladding, cold water',
                                    universe_id=1)
        c21 = openmc.Cell(cell_id=21)
        c21.region = -s1
        c21.fill = fuel
        c22 = openmc.Cell(cell_id=22)
        c22.region = +s1 & -s2
        c22.fill = clad
        c23 = openmc.Cell(cell_id=23)
        c23.region = +s2
        c23.fill = cold_water
        fuel_cold.add_cells((c21, c22, c23))

        tube_cold = openmc.Universe(name='Instrumentation guide tube, '
                                    'cold water', universe_id=2)
        c24 = openmc.Cell(cell_id=24)
        c24.region = -s3
        c24.fill = cold_water
        c25 = openmc.Cell(cell_id=25)
        c25.region = +s3 & -s4
        c25.fill = clad
        c26 = openmc.Cell(cell_id=26)
        c26.region = +s4
        c26.fill = cold_water
        tube_cold.add_cells((c24, c25, c26))

        fuel_hot = openmc.Universe(name='Fuel pin, cladding, hot water',
                                   universe_id=3)
        c27 = openmc.Cell(cell_id=27)
        c27.region = -s1
        c27.fill = fuel
        c28 = openmc.Cell(cell_id=28)
        c28.region = +s1 & -s2
        c28.fill = clad
        c29 = openmc.Cell(cell_id=29)
        c29.region = +s2
        c29.fill = hot_water
        fuel_hot.add_cells((c27, c28, c29))

        tube_hot = openmc.Universe(name='Instrumentation guide tube, hot water',
                                   universe_id=4)
        c30 = openmc.Cell(cell_id=30)
        c30.region = -s3
        c30.fill = hot_water
        c31 = openmc.Cell(cell_id=31)
        c31.region = +s3 & -s4
        c31.fill = clad
        c32 = openmc.Cell(cell_id=32)
        c32.region = +s4
        c32.fill = hot_water
        tube_hot.add_cells((c30, c31, c32))

        # Define fuel lattices.
        l100 = openmc.RectLattice(name='Fuel assembly (lower half)',
                                  lattice_id=100)
        l100.lower_left = (-10.71, -10.71)
        l100.pitch = (1.26, 1.26)
        l100.universes = [
             [fuel_cold]*17,
             [fuel_cold]*17,
             [fuel_cold]*5 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*5,
             [fuel_cold]*3 + [tube_cold] + [fuel_cold]*9 + [tube_cold]
                + [fuel_cold]*3,
             [fuel_cold]*17,
             [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2,
             [fuel_cold]*17,
             [fuel_cold]*17,
             [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2,
             [fuel_cold]*17,
             [fuel_cold]*17,
             [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*2,
             [fuel_cold]*17,
             [fuel_cold]*3 + [tube_cold] + [fuel_cold]*9 + [tube_cold]
                + [fuel_cold]*3,
             [fuel_cold]*5 + [tube_cold] + [fuel_cold]*2 + [tube_cold]
                + [fuel_cold]*2 + [tube_cold] + [fuel_cold]*5,
             [fuel_cold]*17,
             [fuel_cold]*17 ]

        l101 = openmc.RectLattice(name='Fuel assembly (upper half)',
                                  lattice_id=101)
        l101.lower_left = (-10.71, -10.71)
        l101.pitch = (1.26, 1.26)
        l101.universes = [
             [fuel_hot]*17,
             [fuel_hot]*17,
             [fuel_hot]*5 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*5,
             [fuel_hot]*3 + [tube_hot] + [fuel_hot]*9 + [tube_hot]
                + [fuel_hot]*3,
             [fuel_hot]*17,
             [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2,
             [fuel_hot]*17,
             [fuel_hot]*17,
             [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2,
             [fuel_hot]*17,
             [fuel_hot]*17,
             [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*2,
             [fuel_hot]*17,
             [fuel_hot]*3 + [tube_hot] + [fuel_hot]*9 + [tube_hot]
                + [fuel_hot]*3,
             [fuel_hot]*5 + [tube_hot] + [fuel_hot]*2 + [tube_hot]
                + [fuel_hot]*2 + [tube_hot] + [fuel_hot]*5,
             [fuel_hot]*17,
             [fuel_hot]*17 ]

        # Define assemblies.
        fa_cw = openmc.Universe(name='Water assembly (cold)', universe_id=5)
        c50 = openmc.Cell(cell_id=50)
        c50.region = +s34 & -s35
        c50.fill = cold_water
        fa_cw.add_cells((c50, ))

        fa_hw = openmc.Universe(name='Water assembly (hot)', universe_id=7)
        c70 = openmc.Cell(cell_id=70)
        c70.region = +s35 & -s36
        c70.fill = hot_water
        fa_hw.add_cells((c70, ))

        fa_cold = openmc.Universe(name='Fuel assembly (cold)', universe_id=6)
        c60 = openmc.Cell(cell_id=60)
        c60.region = +s34 & -s35
        c60.fill = l100
        fa_cold.add_cells((c60, ))

        fa_hot = openmc.Universe(name='Fuel assembly (hot)', universe_id=8)
        c80 = openmc.Cell(cell_id=80)
        c80.region = +s35 & -s36
        c80.fill = l101
        fa_hot.add_cells((c80, ))

        # Define core lattices
        l200 = openmc.RectLattice(name='Core lattice (lower half)',
                                  lattice_id=200)
        l200.lower_left = (-224.91, -224.91)
        l200.pitch = (21.42, 21.42)
        l200.universes = [
             [fa_cw]*21,
             [fa_cw]*21,
             [fa_cw]*7 + [fa_cold]*7 + [fa_cw]*7,
             [fa_cw]*5 + [fa_cold]*11 + [fa_cw]*5,
             [fa_cw]*4 + [fa_cold]*13 + [fa_cw]*4,
             [fa_cw]*3 + [fa_cold]*15 + [fa_cw]*3,
             [fa_cw]*3 + [fa_cold]*15 + [fa_cw]*3,
             [fa_cw]*2 + [fa_cold]*17 + [fa_cw]*2,
             [fa_cw]*2 + [fa_cold]*17 + [fa_cw]*2,
             [fa_cw]*2 + [fa_cold]*17 + [fa_cw]*2,
             [fa_cw]*2 + [fa_cold]*17 + [fa_cw]*2,
             [fa_cw]*2 + [fa_cold]*17 + [fa_cw]*2,
             [fa_cw]*2 + [fa_cold]*17 + [fa_cw]*2,
             [fa_cw]*2 + [fa_cold]*17 + [fa_cw]*2,
             [fa_cw]*3 + [fa_cold]*15 + [fa_cw]*3,
             [fa_cw]*3 + [fa_cold]*15 + [fa_cw]*3,
             [fa_cw]*4 + [fa_cold]*13 + [fa_cw]*4,
             [fa_cw]*5 + [fa_cold]*11 + [fa_cw]*5,
             [fa_cw]*7 + [fa_cold]*7 + [fa_cw]*7,
             [fa_cw]*21,
             [fa_cw]*21]

        l201 = openmc.RectLattice(name='Core lattice (lower half)',
                                  lattice_id=201)
        l201.lower_left = (-224.91, -224.91)
        l201.pitch = (21.42, 21.42)
        l201.universes = [
             [fa_hw]*21,
             [fa_hw]*21,
             [fa_hw]*7 + [fa_hot]*7 + [fa_hw]*7,
             [fa_hw]*5 + [fa_hot]*11 + [fa_hw]*5,
             [fa_hw]*4 + [fa_hot]*13 + [fa_hw]*4,
             [fa_hw]*3 + [fa_hot]*15 + [fa_hw]*3,
             [fa_hw]*3 + [fa_hot]*15 + [fa_hw]*3,
             [fa_hw]*2 + [fa_hot]*17 + [fa_hw]*2,
             [fa_hw]*2 + [fa_hot]*17 + [fa_hw]*2,
             [fa_hw]*2 + [fa_hot]*17 + [fa_hw]*2,
             [fa_hw]*2 + [fa_hot]*17 + [fa_hw]*2,
             [fa_hw]*2 + [fa_hot]*17 + [fa_hw]*2,
             [fa_hw]*2 + [fa_hot]*17 + [fa_hw]*2,
             [fa_hw]*2 + [fa_hot]*17 + [fa_hw]*2,
             [fa_hw]*3 + [fa_hot]*15 + [fa_hw]*3,
             [fa_hw]*3 + [fa_hot]*15 + [fa_hw]*3,
             [fa_hw]*4 + [fa_hot]*13 + [fa_hw]*4,
             [fa_hw]*5 + [fa_hot]*11 + [fa_hw]*5,
             [fa_hw]*7 + [fa_hot]*7 + [fa_hw]*7,
             [fa_hw]*21,
             [fa_hw]*21]

        # Define root universe.
        root = openmc.Universe(universe_id=0, name='root universe')
        c1 = openmc.Cell(cell_id=1)
        c1.region = -s6 & +s34 & -s35
        c1.fill = l200

        c2 = openmc.Cell(cell_id=2)
        c2.region = -s6 & +s35 & -s36
        c2.fill = l201

        c3 = openmc.Cell(cell_id=3)
        c3.region = -s7 & +s31 & -s32
        c3.fill = bot_plate

        c4 = openmc.Cell(cell_id=4)
        c4.region = -s5 & +s32 & -s33
        c4.fill = bot_nozzle

        c5 = openmc.Cell(cell_id=5)
        c5.region = -s5 & +s33 & -s34
        c5.fill = bot_fa

        c6 = openmc.Cell(cell_id=6)
        c6.region = -s5 & +s36 & -s37
        c6.fill = top_fa

        c7 = openmc.Cell(cell_id=7)
        c7.region = -s5 & +s37 & -s38
        c7.fill = top_nozzle

        c8 = openmc.Cell(cell_id=8)
        c8.region = -s7 & +s38 & -s39
        c8.fill = upper_rad_ref

        c9 = openmc.Cell(cell_id=9)
        c9.region = +s6 & -s7 & +s32 & -s38
        c9.fill = bot_nozzle

        c10 = openmc.Cell(cell_id=10)
        c10.region = +s7 & -s8 & +s31 & -s39
        c10.fill = rpv_steel

        c11 = openmc.Cell(cell_id=11)
        c11.region = +s5 & -s6 & +s32 & -s34
        c11.fill = lower_rad_ref

        c12 = openmc.Cell(cell_id=12)
        c12.region = +s5 & -s6 & +s36 & -s38
        c12.fill = upper_rad_ref

        root.add_cells((c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12))

        # Assign root universe to geometry
        self.geometry.root_universe = root

    def build_default_settings(self):
        self.settings.batches = 10
        self.settings.inactive = 5
        self.settings.particles = 100
        self.settings.source = Source(space=Box(
            [-160, -160, -183], [160, 160, 183]))

    def build_defualt_plots(self):
        plot = openmc.Plot()
        plot.filename = 'mat'
        plot.origin = (125, 125, 0)
        plot.width = (250, 250)
        plot.pixels = (3000, 3000)
        plot.color = 'mat'

        self.plots.add_plot(plot)


class PinCellInputSet(object):
    def __init__(self):
        self.settings = openmc.Settings()
        self.materials = openmc.Materials()
        self.geometry = openmc.Geometry()
        self.tallies = None
        self.plots = None

    def export(self):
        self.settings.export_to_xml()
        self.materials.export_to_xml()
        self.geometry.export_to_xml()
        if self.tallies is not None:
            self.tallies.export_to_xml()
        if self.plots is not None:
            self.plots.export_to_xml()

    def build_default_materials_and_geometry(self):
        # Define materials.
        fuel = openmc.Material(name='Fuel')
        fuel.set_density('g/cm3', 10.29769)
        fuel.add_nuclide("U-234", 4.4843e-6)
        fuel.add_nuclide("U-235", 5.5815e-4)
        fuel.add_nuclide("U-238", 2.2408e-2)
        fuel.add_nuclide("O-16", 4.5829e-2)

        clad = openmc.Material(name='Cladding')
        clad.set_density('g/cm3', 6.55)
        clad.add_nuclide("Zr-90", 2.1827e-2)
        clad.add_nuclide("Zr-91", 4.7600e-3)
        clad.add_nuclide("Zr-92", 7.2758e-3)
        clad.add_nuclide("Zr-94", 7.3734e-3)
        clad.add_nuclide("Zr-96", 1.1879e-3)

        hot_water = openmc.Material(name='Hot borated water')
        hot_water.set_density('g/cm3', 0.740582)
        hot_water.add_nuclide("H-1", 4.9457e-2)
        hot_water.add_nuclide("O-16", 2.4672e-2)
        hot_water.add_nuclide("B-10", 8.0042e-6)
        hot_water.add_nuclide("B-11", 3.2218e-5)
        hot_water.add_s_alpha_beta('HH2O', '71t')

        # Define the materials file.
        self.materials.default_xs = '71c'
        self.materials += (fuel, clad, hot_water)

        # Instantiate ZCylinder surfaces
        fuel_or = openmc.ZCylinder(x0=0, y0=0, R=0.39218, name='Fuel OR')
        clad_or = openmc.ZCylinder(x0=0, y0=0, R=0.45720, name='Clad OR')
        left = openmc.XPlane(x0=-0.63, name='left')
        right = openmc.XPlane(x0=0.63, name='right')
        bottom = openmc.YPlane(y0=-0.63, name='bottom')
        top = openmc.YPlane(y0=0.63, name='top')

        left.boundary_type = 'reflective'
        right.boundary_type = 'reflective'
        top.boundary_type = 'reflective'
        bottom.boundary_type = 'reflective'

        # Instantiate Cells
        fuel_pin = openmc.Cell(name='cell 1')
        cladding = openmc.Cell(name='cell 3')
        water = openmc.Cell(name='cell 2')

        # Use surface half-spaces to define regions
        fuel_pin.region = -fuel_or
        cladding.region = +fuel_or & -clad_or
        water.region = +clad_or & +left & -right & +bottom & -top

        # Register Materials with Cells
        fuel_pin.fill = fuel
        cladding.fill = clad
        water.fill = hot_water

        # Instantiate Universe
        root = openmc.Universe(universe_id=0, name='root universe')

        # Register Cells with Universe
        root.add_cells([fuel_pin, cladding, water])

        # Instantiate a Geometry, register the root Universe, and export to XML
        self.geometry.root_universe = root

    def build_default_settings(self):
        self.settings.batches = 10
        self.settings.inactive = 5
        self.settings.particles = 100
        self.settings.source = Source(space=Box([-0.63, -0.63, -1],
                                                [0.63, 0.63, 1],
                                                only_fissionable=True))

    def build_defualt_plots(self):
        plot = openmc.Plot()
        plot.filename = 'mat'
        plot.origin = (0.0, 0.0, 0)
        plot.width = (1.26, 1.26)
        plot.pixels = (300, 300)
        plot.color = 'mat'

        self.plots.add_plot(plot)


class MGInputSet(InputSet):
    def build_default_materials_and_geometry(self):
        # Define materials needed for 1D/1G slab problem
        uo2_data = openmc.Macroscopic('uo2_iso', '71c')
        uo2 = openmc.Material(name='UO2', material_id=1)
        uo2.set_density('macro', 1.0)
        uo2.add_macroscopic(uo2_data)

        clad_data = openmc.Macroscopic('clad_ang_mu', '71c')
        clad = openmc.Material(name='Clad', material_id=2)
        clad.set_density('macro', 1.0)
        clad.add_macroscopic(clad_data)

        water_data = openmc.Macroscopic('lwtr_iso_mu', '71c')
        water = openmc.Material(name='LWTR', material_id=3)
        water.set_density('macro', 1.0)
        water.add_macroscopic(water_data)

        # Define the materials file.
        self.materials.default_xs = '71c'
        self.materials += (uo2, clad, water)

        # Define surfaces.

        # Assembly/Problem Boundary
        left = openmc.XPlane(x0=0.0, surface_id=200,
                             boundary_type='reflective')
        right = openmc.XPlane(x0=10.0, surface_id=201,
                              boundary_type='reflective')
        bottom = openmc.YPlane(y0=0.0, surface_id=300,
                               boundary_type='reflective')
        top = openmc.YPlane(y0=10.0, surface_id=301,
                            boundary_type='reflective')

        down = openmc.ZPlane(z0=0.0, surface_id=0,
                             boundary_type='reflective')
        fuel_clad_intfc = openmc.ZPlane(z0=2.0, surface_id=1)
        clad_lwtr_intfc = openmc.ZPlane(z0=2.4, surface_id=2)
        up = openmc.ZPlane(z0=5.0, surface_id=3,
                           boundary_type='reflective')

        # Define cells
        c1 = openmc.Cell(cell_id=1)
        c1.region = +left & -right & +bottom & -top & +down & -fuel_clad_intfc
        c1.fill = uo2
        c2 = openmc.Cell(cell_id=2)
        c2.region = +left & -right & +bottom & -top & +fuel_clad_intfc & -clad_lwtr_intfc
        c2.fill = clad
        c3 = openmc.Cell(cell_id=3)
        c3.region = +left & -right & +bottom & -top & +clad_lwtr_intfc & -up
        c3.fill = water

        # Define root universe.
        root = openmc.Universe(universe_id=0, name='root universe')

        root.add_cells((c1, c2, c3))

        # Assign root universe to geometry
        self.geometry.root_universe = root

    def build_default_settings(self):
        self.settings.batches = 10
        self.settings.inactive = 5
        self.settings.particles = 100
        self.settings.source = Source(space=Box([0.0, 0.0, 0.0],
                                                [10.0, 10.0, 2.0]))
        self.settings.energy_mode = "multi-group"
        self.settings.cross_sections = "../1d_mgxs.xml"

    def build_defualt_plots(self):
        plot = openmc.Plot()
        plot.filename = 'mat'
        plot.origin = (5.0, 5.0, 2.5)
        plot.width = (2.5, 2.5)
        plot.basis = 'xz'
        plot.pixels = (3000, 3000)
        plot.color = 'mat'

        self.plots.add_plot(plot)
