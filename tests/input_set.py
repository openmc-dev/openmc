import numpy as np

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
        fuel.add_nuclide("U234", 4.9476e-6)
        fuel.add_nuclide("U235", 4.8218e-4)
        fuel.add_nuclide("U238", 2.1504e-2)
        fuel.add_nuclide("Xe135", 1.0801e-8)
        fuel.add_nuclide("O16", 4.5737e-2)

        clad = openmc.Material(name='Cladding', material_id=2)
        clad.set_density('g/cm3', 5.77)
        clad.add_nuclide("Zr90", 0.5145)
        clad.add_nuclide("Zr91", 0.1122)
        clad.add_nuclide("Zr92", 0.1715)
        clad.add_nuclide("Zr94", 0.1738)
        clad.add_nuclide("Zr96", 0.0280)

        cold_water = openmc.Material(name='Cold borated water', material_id=3)
        cold_water.set_density('atom/b-cm', 0.07416)
        cold_water.add_nuclide("H1", 2.0)
        cold_water.add_nuclide("O16", 1.0)
        cold_water.add_nuclide("B10", 6.490e-4)
        cold_water.add_nuclide("B11", 2.689e-3)
        cold_water.add_s_alpha_beta('c_H_in_H2O')

        hot_water = openmc.Material(name='Hot borated water', material_id=4)
        hot_water.set_density('atom/b-cm', 0.06614)
        hot_water.add_nuclide("H1", 2.0)
        hot_water.add_nuclide("O16", 1.0)
        hot_water.add_nuclide("B10", 6.490e-4)
        hot_water.add_nuclide("B11", 2.689e-3)
        hot_water.add_s_alpha_beta('c_H_in_H2O')

        rpv_steel = openmc.Material(name='Reactor pressure vessel steel',
                                    material_id=5)
        rpv_steel.set_density('g/cm3', 7.9)
        rpv_steel.add_nuclide("Fe54", 0.05437098, 'wo')
        rpv_steel.add_nuclide("Fe56", 0.88500663, 'wo')
        rpv_steel.add_nuclide("Fe57", 0.0208008, 'wo')
        rpv_steel.add_nuclide("Fe58", 0.00282159, 'wo')
        rpv_steel.add_nuclide("Ni58", 0.0067198, 'wo')
        rpv_steel.add_nuclide("Ni60", 0.0026776, 'wo')
        rpv_steel.add_nuclide("Mn55", 0.01, 'wo')
        rpv_steel.add_nuclide("Cr52", 0.002092475, 'wo')
        rpv_steel.add_nuclide("C0", 0.0025, 'wo')
        rpv_steel.add_nuclide("Cu63", 0.0013696, 'wo')

        lower_rad_ref = openmc.Material(name='Lower radial reflector',
                                        material_id=6)
        lower_rad_ref.set_density('g/cm3', 4.32)
        lower_rad_ref.add_nuclide("H1", 0.0095661, 'wo')
        lower_rad_ref.add_nuclide("O16", 0.0759107, 'wo')
        lower_rad_ref.add_nuclide("B10", 3.08409e-5, 'wo')
        lower_rad_ref.add_nuclide("B11", 1.40499e-4, 'wo')
        lower_rad_ref.add_nuclide("Fe54", 0.035620772088, 'wo')
        lower_rad_ref.add_nuclide("Fe56", 0.579805982228, 'wo')
        lower_rad_ref.add_nuclide("Fe57", 0.01362750048, 'wo')
        lower_rad_ref.add_nuclide("Fe58", 0.001848545204, 'wo')
        lower_rad_ref.add_nuclide("Ni58", 0.055298376566, 'wo')
        lower_rad_ref.add_nuclide("Mn55", 0.0182870, 'wo')
        lower_rad_ref.add_nuclide("Cr52", 0.145407678031, 'wo')
        lower_rad_ref.add_s_alpha_beta('c_H_in_H2O')

        upper_rad_ref = openmc.Material(name='Upper radial reflector /'
                                             'Top plate region', material_id=7)
        upper_rad_ref.set_density('g/cm3', 4.28)
        upper_rad_ref.add_nuclide("H1", 0.0086117, 'wo')
        upper_rad_ref.add_nuclide("O16", 0.0683369, 'wo')
        upper_rad_ref.add_nuclide("B10", 2.77638e-5, 'wo')
        upper_rad_ref.add_nuclide("B11", 1.26481e-4, 'wo')
        upper_rad_ref.add_nuclide("Fe54", 0.035953677186, 'wo')
        upper_rad_ref.add_nuclide("Fe56", 0.585224740891, 'wo')
        upper_rad_ref.add_nuclide("Fe57", 0.01375486056, 'wo')
        upper_rad_ref.add_nuclide("Fe58", 0.001865821363, 'wo')
        upper_rad_ref.add_nuclide("Ni58", 0.055815129186, 'wo')
        upper_rad_ref.add_nuclide("Mn55", 0.0184579, 'wo')
        upper_rad_ref.add_nuclide("Cr52", 0.146766614995, 'wo')
        upper_rad_ref.add_s_alpha_beta('c_H_in_H2O')

        bot_plate = openmc.Material(name='Bottom plate region', material_id=8)
        bot_plate.set_density('g/cm3', 7.184)
        bot_plate.add_nuclide("H1", 0.0011505, 'wo')
        bot_plate.add_nuclide("O16", 0.0091296, 'wo')
        bot_plate.add_nuclide("B10", 3.70915e-6, 'wo')
        bot_plate.add_nuclide("B11", 1.68974e-5, 'wo')
        bot_plate.add_nuclide("Fe54", 0.03855611055, 'wo')
        bot_plate.add_nuclide("Fe56", 0.627585036425, 'wo')
        bot_plate.add_nuclide("Fe57", 0.014750478, 'wo')
        bot_plate.add_nuclide("Fe58", 0.002000875025, 'wo')
        bot_plate.add_nuclide("Ni58", 0.059855207342, 'wo')
        bot_plate.add_nuclide("Mn55", 0.0197940, 'wo')
        bot_plate.add_nuclide("Cr52", 0.157390026871, 'wo')
        bot_plate.add_s_alpha_beta('c_H_in_H2O')

        bot_nozzle = openmc.Material(name='Bottom nozzle region',
                                     material_id=9)
        bot_nozzle.set_density('g/cm3', 2.53)
        bot_nozzle.add_nuclide("H1", 0.0245014, 'wo')
        bot_nozzle.add_nuclide("O16", 0.1944274, 'wo')
        bot_nozzle.add_nuclide("B10", 7.89917e-5, 'wo')
        bot_nozzle.add_nuclide("B11", 3.59854e-4, 'wo')
        bot_nozzle.add_nuclide("Fe54", 0.030411411144, 'wo')
        bot_nozzle.add_nuclide("Fe56", 0.495012237964, 'wo')
        bot_nozzle.add_nuclide("Fe57", 0.01163454624, 'wo')
        bot_nozzle.add_nuclide("Fe58", 0.001578204652, 'wo')
        bot_nozzle.add_nuclide("Ni58", 0.047211231662, 'wo')
        bot_nozzle.add_nuclide("Mn55", 0.0156126, 'wo')
        bot_nozzle.add_nuclide("Cr52", 0.124142524198, 'wo')
        bot_nozzle.add_s_alpha_beta('c_H_in_H2O')

        top_nozzle = openmc.Material(name='Top nozzle region', material_id=10)
        top_nozzle.set_density('g/cm3', 1.746)
        top_nozzle.add_nuclide("H1", 0.0358870, 'wo')
        top_nozzle.add_nuclide("O16", 0.2847761, 'wo')
        top_nozzle.add_nuclide("B10", 1.15699e-4, 'wo')
        top_nozzle.add_nuclide("B11", 5.27075e-4, 'wo')
        top_nozzle.add_nuclide("Fe54", 0.02644016154, 'wo')
        top_nozzle.add_nuclide("Fe56", 0.43037146399, 'wo')
        top_nozzle.add_nuclide("Fe57", 0.0101152584, 'wo')
        top_nozzle.add_nuclide("Fe58", 0.00137211607, 'wo')
        top_nozzle.add_nuclide("Ni58", 0.04104621835, 'wo')
        top_nozzle.add_nuclide("Mn55", 0.0135739, 'wo')
        top_nozzle.add_nuclide("Cr52", 0.107931450781, 'wo')
        top_nozzle.add_s_alpha_beta('c_H_in_H2O')

        top_fa = openmc.Material(name='Top of fuel assemblies', material_id=11)
        top_fa.set_density('g/cm3', 3.044)
        top_fa.add_nuclide("H1", 0.0162913, 'wo')
        top_fa.add_nuclide("O16", 0.1292776, 'wo')
        top_fa.add_nuclide("B10", 5.25228e-5, 'wo')
        top_fa.add_nuclide("B11", 2.39272e-4, 'wo')
        top_fa.add_nuclide("Zr90", 0.43313403903, 'wo')
        top_fa.add_nuclide("Zr91", 0.09549277374, 'wo')
        top_fa.add_nuclide("Zr92", 0.14759527104, 'wo')
        top_fa.add_nuclide("Zr94", 0.15280552077, 'wo')
        top_fa.add_nuclide("Zr96", 0.02511169542, 'wo')
        top_fa.add_s_alpha_beta('c_H_in_H2O')

        bot_fa = openmc.Material(name='Bottom of fuel assemblies',
                                 material_id=12)
        bot_fa.set_density('g/cm3', 1.762)
        bot_fa.add_nuclide("H1", 0.0292856, 'wo')
        bot_fa.add_nuclide("O16", 0.2323919, 'wo')
        bot_fa.add_nuclide("B10", 9.44159e-5, 'wo')
        bot_fa.add_nuclide("B11", 4.30120e-4, 'wo')
        bot_fa.add_nuclide("Zr90", 0.3741373658, 'wo')
        bot_fa.add_nuclide("Zr91", 0.0824858164, 'wo')
        bot_fa.add_nuclide("Zr92", 0.1274914944, 'wo')
        bot_fa.add_nuclide("Zr94", 0.1319920622, 'wo')
        bot_fa.add_nuclide("Zr96", 0.0216912612, 'wo')
        bot_fa.add_s_alpha_beta('c_H_in_H2O')

        # Define the materials file.
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
        fuel.add_nuclide("U234", 4.4843e-6)
        fuel.add_nuclide("U235", 5.5815e-4)
        fuel.add_nuclide("U238", 2.2408e-2)
        fuel.add_nuclide("O16", 4.5829e-2)

        clad = openmc.Material(name='Cladding')
        clad.set_density('g/cm3', 6.55)
        clad.add_nuclide("Zr90", 2.1827e-2)
        clad.add_nuclide("Zr91", 4.7600e-3)
        clad.add_nuclide("Zr92", 7.2758e-3)
        clad.add_nuclide("Zr94", 7.3734e-3)
        clad.add_nuclide("Zr96", 1.1879e-3)

        hot_water = openmc.Material(name='Hot borated water')
        hot_water.set_density('g/cm3', 0.740582)
        hot_water.add_nuclide("H1", 4.9457e-2)
        hot_water.add_nuclide("O16", 2.4672e-2)
        hot_water.add_nuclide("B10", 8.0042e-6)
        hot_water.add_nuclide("B11", 3.2218e-5)
        hot_water.add_s_alpha_beta('c_H_in_H2O')

        # Define the materials file.
        self.materials += (fuel, clad, hot_water)

        # Instantiate ZCylinder surfaces
        fuel_or = openmc.ZCylinder(x0=0, y0=0, R=0.39218, name='Fuel OR')
        clad_or = openmc.ZCylinder(x0=0, y0=0, R=0.45720, name='Clad OR')
        left = openmc.XPlane(x0=-0.63, name='left', boundary_type='reflective')
        right = openmc.XPlane(x0=0.63, name='right', boundary_type='reflective')
        bottom = openmc.YPlane(y0=-0.63, name='bottom',
                               boundary_type='reflective')
        top = openmc.YPlane(y0=0.63, name='top', boundary_type='reflective')

        # Instantiate Cells
        fuel_pin = openmc.Cell(name='cell 1', fill=fuel)
        cladding = openmc.Cell(name='cell 3', fill=clad)
        water = openmc.Cell(name='cell 2', fill=hot_water)

        # Use surface half-spaces to define regions
        fuel_pin.region = -fuel_or
        cladding.region = +fuel_or & -clad_or
        water.region = +clad_or & +left & -right & +bottom & -top

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


class AssemblyInputSet(object):
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
        fuel.add_nuclide("U234", 4.4843e-6)
        fuel.add_nuclide("U235", 5.5815e-4)
        fuel.add_nuclide("U238", 2.2408e-2)
        fuel.add_nuclide("O16", 4.5829e-2)

        clad = openmc.Material(name='Cladding')
        clad.set_density('g/cm3', 6.55)
        clad.add_nuclide("Zr90", 2.1827e-2)
        clad.add_nuclide("Zr91", 4.7600e-3)
        clad.add_nuclide("Zr92", 7.2758e-3)
        clad.add_nuclide("Zr94", 7.3734e-3)
        clad.add_nuclide("Zr96", 1.1879e-3)

        hot_water = openmc.Material(name='Hot borated water')
        hot_water.set_density('g/cm3', 0.740582)
        hot_water.add_nuclide("H1", 4.9457e-2)
        hot_water.add_nuclide("O16", 2.4672e-2)
        hot_water.add_nuclide("B10", 8.0042e-6)
        hot_water.add_nuclide("B11", 3.2218e-5)
        hot_water.add_s_alpha_beta('c_H_in_H2O')

        # Define the materials file.
        self.materials += (fuel, clad, hot_water)

        # Instantiate ZCylinder surfaces
        fuel_or = openmc.ZCylinder(x0=0, y0=0, R=0.39218, name='Fuel OR')
        clad_or = openmc.ZCylinder(x0=0, y0=0, R=0.45720, name='Clad OR')

        # Create boundary planes to surround the geometry
        min_x = openmc.XPlane(x0=-10.71, boundary_type='reflective')
        max_x = openmc.XPlane(x0=+10.71, boundary_type='reflective')
        min_y = openmc.YPlane(y0=-10.71, boundary_type='reflective')
        max_y = openmc.YPlane(y0=+10.71, boundary_type='reflective')

        # Create a Universe to encapsulate a fuel pin
        fuel_pin_universe = openmc.Universe(name='Fuel Pin')

        # Create fuel Cell
        fuel_cell = openmc.Cell(name='fuel')
        fuel_cell.fill = fuel
        fuel_cell.region = -fuel_or
        fuel_pin_universe.add_cell(fuel_cell)

        # Create a clad Cell
        clad_cell = openmc.Cell(name='clad')
        clad_cell.fill = clad
        clad_cell.region = +fuel_or & -clad_or
        fuel_pin_universe.add_cell(clad_cell)

        # Create a moderator Cell
        hot_water_cell = openmc.Cell(name='hot water')
        hot_water_cell.fill = hot_water
        hot_water_cell.region = +clad_or
        fuel_pin_universe.add_cell(hot_water_cell)

        # Create a Universe to encapsulate a control rod guide tube
        guide_tube_universe = openmc.Universe(name='Guide Tube')

        # Create guide tube inner Cell
        gt_inner_cell = openmc.Cell(name='guide tube inner water')
        gt_inner_cell.fill = hot_water
        gt_inner_cell.region = -fuel_or
        guide_tube_universe.add_cell(gt_inner_cell)

        # Create a clad Cell
        gt_clad_cell = openmc.Cell(name='guide tube clad')
        gt_clad_cell.fill = clad
        gt_clad_cell.region = +fuel_or & -clad_or
        guide_tube_universe.add_cell(gt_clad_cell)

        # Create a guide tube outer Cell
        gt_outer_cell = openmc.Cell(name='guide tube outer water')
        gt_outer_cell.fill = hot_water
        gt_outer_cell.region = +clad_or
        guide_tube_universe.add_cell(gt_outer_cell)

        # Create fuel assembly Lattice
        assembly = openmc.RectLattice(name='Fuel Assembly')
        assembly.pitch = (1.26, 1.26)
        assembly.lower_left = [-1.26 * 17. / 2.0] * 2

        # Create array indices for guide tube locations in lattice
        template_x = np.array([5, 8, 11, 3, 13, 2, 5, 8, 11, 14, 2, 5, 8,
                               11, 14, 2, 5, 8, 11, 14, 3, 13, 5, 8, 11])
        template_y = np.array([2, 2, 2, 3, 3, 5, 5, 5, 5, 5, 8, 8, 8, 8,
                               8, 11, 11, 11, 11, 11, 13, 13, 14, 14, 14])

        # Initialize an empty 17x17 array of the lattice universes
        universes = np.empty((17, 17), dtype=openmc.Universe)

        # Fill the array with the fuel pin and guide tube universes
        universes[:,:] = fuel_pin_universe
        universes[template_x, template_y] = guide_tube_universe

        # Store the array of universes in the lattice
        assembly.universes = universes

        # Create root Cell
        root_cell = openmc.Cell(name='root cell')
        root_cell.fill = assembly

        # Add boundary planes
        root_cell.region = +min_x & -max_x & +min_y & -max_y

        # Create root Universe
        root_universe = openmc.Universe(universe_id=0, name='root universe')
        root_universe.add_cell(root_cell)

        # Instantiate a Geometry, register the root Universe, and export to XML
        self.geometry.root_universe = root_universe

    def build_default_settings(self):
        self.settings.batches = 10
        self.settings.inactive = 5
        self.settings.particles = 100
        self.settings.source = Source(space=Box([-10.71, -10.71, -1],
                                                [10.71, 10.71, 1],
                                                only_fissionable=True))

    def build_defualt_plots(self):
        plot = openmc.Plot()
        plot.filename = 'mat'
        plot.origin = (0.0, 0.0, 0)
        plot.width = (21.42, 21.42)
        plot.pixels = (300, 300)
        plot.color = 'mat'

        self.plots.add_plot(plot)


class MGInputSet(InputSet):
    def build_default_materials_and_geometry(self, reps=None, as_macro=True):
        # Define materials needed for 1D/1G slab problem
        mat_names = ['uo2', 'clad', 'lwtr']
        mgxs_reps = ['ang', 'ang_mu', 'iso', 'iso_mu']

        if reps is None:
            reps = mgxs_reps

        xs = []
        mats = []
        i = 0
        for mat in mat_names:
            for rep in reps:
                i += 1
                if as_macro:
                    xs.append(openmc.Macroscopic(mat + '_' + rep))
                    mats.append(openmc.Material(name=str(i)))
                    mats[-1].set_density('macro', 1.)
                    mats[-1].add_macroscopic(xs[-1])
                else:
                    xs.append(openmc.Nuclide(mat + '_' + rep))
                    mats.append(openmc.Material(name=str(i)))
                    mats[-1].set_density('atom/b-cm', 1.)
                    mats[-1].add_nuclide(xs[-1].name, 1.0, 'ao')

        # Define the materials file
        self.xs_data = xs
        self.materials += mats
        self.materials.cross_sections = "../1d_mgxs.h5"

        # Define surfaces.
        # Assembly/Problem Boundary
        left = openmc.XPlane(x0=0.0, boundary_type='reflective')
        right = openmc.XPlane(x0=10.0, boundary_type='reflective')
        bottom = openmc.YPlane(y0=0.0, boundary_type='reflective')
        top = openmc.YPlane(y0=10.0, boundary_type='reflective')
        # for each material add a plane
        planes = [openmc.ZPlane(z0=0.0, boundary_type='reflective')]
        dz = round(5. / float(len(mats)), 4)
        for i in range(len(mats) - 1):
            planes.append(openmc.ZPlane(z0=dz * float(i + 1)))
        planes.append(openmc.ZPlane(z0=5.0, boundary_type='reflective'))

        # Define cells for each material
        cells = []
        xy = +left & -right & +bottom & -top
        for i, mat in enumerate(mats):
            cells.append(openmc.Cell())
            cells[-1].region = xy & +planes[i] & -planes[i + 1]
            cells[-1].fill = mat

        # Define root universe.
        root = openmc.Universe(universe_id=0, name='root universe')
        root.add_cells(cells)

        # Assign root universe to geometry
        self.geometry.root_universe = root

    def build_default_settings(self):
        self.settings.batches = 10
        self.settings.inactive = 5
        self.settings.particles = 100
        self.settings.source = Source(space=Box([0.0, 0.0, 0.0],
                                                [10.0, 10.0, 5.]))
        self.settings.energy_mode = "multi-group"

    def build_defualt_plots(self):
        plot = openmc.Plot()
        plot.filename = 'mat'
        plot.origin = (5.0, 5.0, 2.5)
        plot.width = (2.5, 2.5)
        plot.basis = 'xz'
        plot.pixels = (3000, 3000)
        plot.color = 'mat'

        self.plots.add_plot(plot)
