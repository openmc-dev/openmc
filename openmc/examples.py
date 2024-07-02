from numbers import Integral

import numpy as np

import openmc
import openmc.model


def pwr_pin_cell():
    """Create a PWR pin-cell model.

    This model is a single fuel pin with 2.4 w/o enriched UO2 corresponding to a
    beginning-of-cycle condition and borated water. The specifications are from
    the `BEAVRS <http://crpg.mit.edu/research/beavrs>`_ benchmark. Note that the
    number of particles/batches is initially set very low for testing purposes.

    Returns
    -------
    model : openmc.model.Model
        A PWR pin-cell model

    """
    model = openmc.model.Model()

    # Define materials.
    fuel = openmc.Material(name='UO2 (2.4%)')
    fuel.set_density('g/cm3', 10.29769)
    fuel.add_nuclide('U234', 4.4843e-6)
    fuel.add_nuclide('U235', 5.5815e-4)
    fuel.add_nuclide('U238', 2.2408e-2)
    fuel.add_nuclide('O16', 4.5829e-2)

    clad = openmc.Material(name='Zircaloy')
    clad.set_density('g/cm3', 6.55)
    clad.add_nuclide('Zr90', 2.1827e-2)
    clad.add_nuclide('Zr91', 4.7600e-3)
    clad.add_nuclide('Zr92', 7.2758e-3)
    clad.add_nuclide('Zr94', 7.3734e-3)
    clad.add_nuclide('Zr96', 1.1879e-3)

    hot_water = openmc.Material(name='Hot borated water')
    hot_water.set_density('g/cm3', 0.740582)
    hot_water.add_nuclide('H1', 4.9457e-2)
    hot_water.add_nuclide('O16', 2.4672e-2)
    hot_water.add_nuclide('B10', 8.0042e-6)
    hot_water.add_nuclide('B11', 3.2218e-5)
    hot_water.add_s_alpha_beta('c_H_in_H2O')

    # Define the materials file.
    model.materials = (fuel, clad, hot_water)

    # Instantiate ZCylinder surfaces
    pitch = 1.26
    fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.39218, name='Fuel OR')
    clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.45720, name='Clad OR')
    left = openmc.XPlane(x0=-pitch/2, name='left', boundary_type='reflective')
    right = openmc.XPlane(x0=pitch/2, name='right', boundary_type='reflective')
    bottom = openmc.YPlane(y0=-pitch/2, name='bottom',
                           boundary_type='reflective')
    top = openmc.YPlane(y0=pitch/2, name='top', boundary_type='reflective')

    # Instantiate Cells
    fuel_pin = openmc.Cell(name='Fuel', fill=fuel)
    cladding = openmc.Cell(name='Cladding', fill=clad)
    water = openmc.Cell(name='Water', fill=hot_water)

    # Use surface half-spaces to define regions
    fuel_pin.region = -fuel_or
    cladding.region = +fuel_or & -clad_or
    water.region = +clad_or & +left & -right & +bottom & -top

    # Create root universe
    model.geometry.root_universe = openmc.Universe(0, name='root universe')
    model.geometry.root_universe.add_cells([fuel_pin, cladding, water])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 100
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box([-pitch/2, -pitch/2, -1],
                               [pitch/2, pitch/2, 1]),
        constraints={'fissionable': True}
    )

    plot = openmc.Plot.from_geometry(model.geometry)
    plot.pixels = (300, 300)
    plot.color_by = 'material'
    model.plots.append(plot)

    return model


def pwr_core():
    """Create a PWR full-core model.

    This model is the OECD/NEA Monte Carlo Performance benchmark which is a
    grossly simplified pressurized water reactor (PWR) with 241 fuel
    assemblies. Note that the number of particles/batches is initially set very
    low for testing purposes.

    Returns
    -------
    model : openmc.model.Model
        Full-core PWR model

    """
    model = openmc.model.Model()

    # Define materials.
    fuel = openmc.Material(1, name='UOX fuel')
    fuel.set_density('g/cm3', 10.062)
    fuel.add_nuclide('U234', 4.9476e-6)
    fuel.add_nuclide('U235', 4.8218e-4)
    fuel.add_nuclide('U238', 2.1504e-2)
    fuel.add_nuclide('Xe135', 1.0801e-8)
    fuel.add_nuclide('O16', 4.5737e-2)

    clad = openmc.Material(2, name='Zircaloy')
    clad.set_density('g/cm3', 5.77)
    clad.add_nuclide('Zr90', 0.5145)
    clad.add_nuclide('Zr91', 0.1122)
    clad.add_nuclide('Zr92', 0.1715)
    clad.add_nuclide('Zr94', 0.1738)
    clad.add_nuclide('Zr96', 0.0280)

    cold_water = openmc.Material(3, name='Cold borated water')
    cold_water.set_density('atom/b-cm', 0.07416)
    cold_water.add_nuclide('H1', 2.0)
    cold_water.add_nuclide('O16', 1.0)
    cold_water.add_nuclide('B10', 6.490e-4)
    cold_water.add_nuclide('B11', 2.689e-3)
    cold_water.add_s_alpha_beta('c_H_in_H2O')

    hot_water = openmc.Material(4, name='Hot borated water')
    hot_water.set_density('atom/b-cm', 0.06614)
    hot_water.add_nuclide('H1', 2.0)
    hot_water.add_nuclide('O16', 1.0)
    hot_water.add_nuclide('B10', 6.490e-4)
    hot_water.add_nuclide('B11', 2.689e-3)
    hot_water.add_s_alpha_beta('c_H_in_H2O')

    rpv_steel = openmc.Material(5, name='Reactor pressure vessel steel')
    rpv_steel.set_density('g/cm3', 7.9)
    rpv_steel.add_nuclide('Fe54', 0.05437098, 'wo')
    rpv_steel.add_nuclide('Fe56', 0.88500663, 'wo')
    rpv_steel.add_nuclide('Fe57', 0.0208008, 'wo')
    rpv_steel.add_nuclide('Fe58', 0.00282159, 'wo')
    rpv_steel.add_nuclide('Ni58', 0.0067198, 'wo')
    rpv_steel.add_nuclide('Ni60', 0.0026776, 'wo')
    rpv_steel.add_nuclide('Mn55', 0.01, 'wo')
    rpv_steel.add_nuclide('Cr52', 0.002092475, 'wo')
    rpv_steel.add_nuclide('C0', 0.0025, 'wo')
    rpv_steel.add_nuclide('Cu63', 0.0013696, 'wo')

    lower_rad_ref = openmc.Material(6, name='Lower radial reflector')
    lower_rad_ref.set_density('g/cm3', 4.32)
    lower_rad_ref.add_nuclide('H1', 0.0095661, 'wo')
    lower_rad_ref.add_nuclide('O16', 0.0759107, 'wo')
    lower_rad_ref.add_nuclide('B10', 3.08409e-5, 'wo')
    lower_rad_ref.add_nuclide('B11', 1.40499e-4, 'wo')
    lower_rad_ref.add_nuclide('Fe54', 0.035620772088, 'wo')
    lower_rad_ref.add_nuclide('Fe56', 0.579805982228, 'wo')
    lower_rad_ref.add_nuclide('Fe57', 0.01362750048, 'wo')
    lower_rad_ref.add_nuclide('Fe58', 0.001848545204, 'wo')
    lower_rad_ref.add_nuclide('Ni58', 0.055298376566, 'wo')
    lower_rad_ref.add_nuclide('Mn55', 0.0182870, 'wo')
    lower_rad_ref.add_nuclide('Cr52', 0.145407678031, 'wo')
    lower_rad_ref.add_s_alpha_beta('c_H_in_H2O')

    upper_rad_ref = openmc.Material(
        7, name='Upper radial reflector / Top plate region')
    upper_rad_ref.set_density('g/cm3', 4.28)
    upper_rad_ref.add_nuclide('H1', 0.0086117, 'wo')
    upper_rad_ref.add_nuclide('O16', 0.0683369, 'wo')
    upper_rad_ref.add_nuclide('B10', 2.77638e-5, 'wo')
    upper_rad_ref.add_nuclide('B11', 1.26481e-4, 'wo')
    upper_rad_ref.add_nuclide('Fe54', 0.035953677186, 'wo')
    upper_rad_ref.add_nuclide('Fe56', 0.585224740891, 'wo')
    upper_rad_ref.add_nuclide('Fe57', 0.01375486056, 'wo')
    upper_rad_ref.add_nuclide('Fe58', 0.001865821363, 'wo')
    upper_rad_ref.add_nuclide('Ni58', 0.055815129186, 'wo')
    upper_rad_ref.add_nuclide('Mn55', 0.0184579, 'wo')
    upper_rad_ref.add_nuclide('Cr52', 0.146766614995, 'wo')
    upper_rad_ref.add_s_alpha_beta('c_H_in_H2O')

    bot_plate = openmc.Material(8, name='Bottom plate region')
    bot_plate.set_density('g/cm3', 7.184)
    bot_plate.add_nuclide('H1', 0.0011505, 'wo')
    bot_plate.add_nuclide('O16', 0.0091296, 'wo')
    bot_plate.add_nuclide('B10', 3.70915e-6, 'wo')
    bot_plate.add_nuclide('B11', 1.68974e-5, 'wo')
    bot_plate.add_nuclide('Fe54', 0.03855611055, 'wo')
    bot_plate.add_nuclide('Fe56', 0.627585036425, 'wo')
    bot_plate.add_nuclide('Fe57', 0.014750478, 'wo')
    bot_plate.add_nuclide('Fe58', 0.002000875025, 'wo')
    bot_plate.add_nuclide('Ni58', 0.059855207342, 'wo')
    bot_plate.add_nuclide('Mn55', 0.0197940, 'wo')
    bot_plate.add_nuclide('Cr52', 0.157390026871, 'wo')
    bot_plate.add_s_alpha_beta('c_H_in_H2O')

    bot_nozzle = openmc.Material(9, name='Bottom nozzle region')
    bot_nozzle.set_density('g/cm3', 2.53)
    bot_nozzle.add_nuclide('H1', 0.0245014, 'wo')
    bot_nozzle.add_nuclide('O16', 0.1944274, 'wo')
    bot_nozzle.add_nuclide('B10', 7.89917e-5, 'wo')
    bot_nozzle.add_nuclide('B11', 3.59854e-4, 'wo')
    bot_nozzle.add_nuclide('Fe54', 0.030411411144, 'wo')
    bot_nozzle.add_nuclide('Fe56', 0.495012237964, 'wo')
    bot_nozzle.add_nuclide('Fe57', 0.01163454624, 'wo')
    bot_nozzle.add_nuclide('Fe58', 0.001578204652, 'wo')
    bot_nozzle.add_nuclide('Ni58', 0.047211231662, 'wo')
    bot_nozzle.add_nuclide('Mn55', 0.0156126, 'wo')
    bot_nozzle.add_nuclide('Cr52', 0.124142524198, 'wo')
    bot_nozzle.add_s_alpha_beta('c_H_in_H2O')

    top_nozzle = openmc.Material(10, name='Top nozzle region')
    top_nozzle.set_density('g/cm3', 1.746)
    top_nozzle.add_nuclide('H1', 0.0358870, 'wo')
    top_nozzle.add_nuclide('O16', 0.2847761, 'wo')
    top_nozzle.add_nuclide('B10', 1.15699e-4, 'wo')
    top_nozzle.add_nuclide('B11', 5.27075e-4, 'wo')
    top_nozzle.add_nuclide('Fe54', 0.02644016154, 'wo')
    top_nozzle.add_nuclide('Fe56', 0.43037146399, 'wo')
    top_nozzle.add_nuclide('Fe57', 0.0101152584, 'wo')
    top_nozzle.add_nuclide('Fe58', 0.00137211607, 'wo')
    top_nozzle.add_nuclide('Ni58', 0.04104621835, 'wo')
    top_nozzle.add_nuclide('Mn55', 0.0135739, 'wo')
    top_nozzle.add_nuclide('Cr52', 0.107931450781, 'wo')
    top_nozzle.add_s_alpha_beta('c_H_in_H2O')

    top_fa = openmc.Material(11, name='Top of fuel assemblies')
    top_fa.set_density('g/cm3', 3.044)
    top_fa.add_nuclide('H1', 0.0162913, 'wo')
    top_fa.add_nuclide('O16', 0.1292776, 'wo')
    top_fa.add_nuclide('B10', 5.25228e-5, 'wo')
    top_fa.add_nuclide('B11', 2.39272e-4, 'wo')
    top_fa.add_nuclide('Zr90', 0.43313403903, 'wo')
    top_fa.add_nuclide('Zr91', 0.09549277374, 'wo')
    top_fa.add_nuclide('Zr92', 0.14759527104, 'wo')
    top_fa.add_nuclide('Zr94', 0.15280552077, 'wo')
    top_fa.add_nuclide('Zr96', 0.02511169542, 'wo')
    top_fa.add_s_alpha_beta('c_H_in_H2O')

    bot_fa = openmc.Material(12, name='Bottom of fuel assemblies')
    bot_fa.set_density('g/cm3', 1.762)
    bot_fa.add_nuclide('H1', 0.0292856, 'wo')
    bot_fa.add_nuclide('O16', 0.2323919, 'wo')
    bot_fa.add_nuclide('B10', 9.44159e-5, 'wo')
    bot_fa.add_nuclide('B11', 4.30120e-4, 'wo')
    bot_fa.add_nuclide('Zr90', 0.3741373658, 'wo')
    bot_fa.add_nuclide('Zr91', 0.0824858164, 'wo')
    bot_fa.add_nuclide('Zr92', 0.1274914944, 'wo')
    bot_fa.add_nuclide('Zr94', 0.1319920622, 'wo')
    bot_fa.add_nuclide('Zr96', 0.0216912612, 'wo')
    bot_fa.add_s_alpha_beta('c_H_in_H2O')

    # Define the materials file.
    model.materials = (fuel, clad, cold_water, hot_water, rpv_steel,
                       lower_rad_ref, upper_rad_ref, bot_plate,
                       bot_nozzle, top_nozzle, top_fa, bot_fa)

    # Define surfaces.
    s1 = openmc.ZCylinder(r=0.41, surface_id=1)
    s2 = openmc.ZCylinder(r=0.475, surface_id=2)
    s3 = openmc.ZCylinder(r=0.56, surface_id=3)
    s4 = openmc.ZCylinder(r=0.62, surface_id=4)
    s5 = openmc.ZCylinder(r=187.6, surface_id=5)
    s6 = openmc.ZCylinder(r=209.0, surface_id=6)
    s7 = openmc.ZCylinder(r=229.0, surface_id=7)
    s8 = openmc.ZCylinder(r=249.0, surface_id=8, boundary_type='vacuum')

    s31 = openmc.ZPlane(z0=-229.0, surface_id=31, boundary_type='vacuum')
    s32 = openmc.ZPlane(z0=-199.0, surface_id=32)
    s33 = openmc.ZPlane(z0=-193.0, surface_id=33)
    s34 = openmc.ZPlane(z0=-183.0, surface_id=34)
    s35 = openmc.ZPlane(z0=0.0, surface_id=35)
    s36 = openmc.ZPlane(z0=183.0, surface_id=36)
    s37 = openmc.ZPlane(z0=203.0, surface_id=37)
    s38 = openmc.ZPlane(z0=215.0, surface_id=38)
    s39 = openmc.ZPlane(z0=223.0, surface_id=39, boundary_type='vacuum')

    # Define pin cells.
    fuel_cold = openmc.Universe(name='Fuel pin, cladding, cold water',
                                universe_id=1)
    c21 = openmc.Cell(cell_id=21, fill=fuel, region=-s1)
    c22 = openmc.Cell(cell_id=22, fill=clad, region=+s1 & -s2)
    c23 = openmc.Cell(cell_id=23, fill=cold_water, region=+s2)
    fuel_cold.add_cells((c21, c22, c23))

    tube_cold = openmc.Universe(name='Instrumentation guide tube, '
                                'cold water', universe_id=2)
    c24 = openmc.Cell(cell_id=24, fill=cold_water, region=-s3)
    c25 = openmc.Cell(cell_id=25, fill=clad, region=+s3 & -s4)
    c26 = openmc.Cell(cell_id=26, fill=cold_water, region=+s4)
    tube_cold.add_cells((c24, c25, c26))

    fuel_hot = openmc.Universe(name='Fuel pin, cladding, hot water',
                               universe_id=3)
    c27 = openmc.Cell(cell_id=27, fill=fuel, region=-s1)
    c28 = openmc.Cell(cell_id=28, fill=clad, region=+s1 & -s2)
    c29 = openmc.Cell(cell_id=29, fill=hot_water, region=+s2)
    fuel_hot.add_cells((c27, c28, c29))

    tube_hot = openmc.Universe(name='Instrumentation guide tube, hot water',
                               universe_id=4)
    c30 = openmc.Cell(cell_id=30, fill=hot_water, region=-s3)
    c31 = openmc.Cell(cell_id=31, fill=clad, region=+s3 & -s4)
    c32 = openmc.Cell(cell_id=32, fill=hot_water, region=+s4)
    tube_hot.add_cells((c30, c31, c32))

    # Set positions occupied by guide tubes
    tube_x = np.array([5, 8, 11, 3, 13, 2, 5, 8, 11, 14, 2, 5, 8, 11, 14,
                       2, 5, 8, 11, 14, 3, 13, 5, 8, 11])
    tube_y = np.array([2, 2, 2, 3, 3, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8,
                       11, 11, 11, 11, 11, 13, 13, 14, 14, 14])

    # Define fuel lattices.
    l100 = openmc.RectLattice(
        name='Fuel assembly (lower half)', lattice_id=100)
    l100.lower_left = (-10.71, -10.71)
    l100.pitch = (1.26, 1.26)
    l100.universes = np.tile(fuel_cold, (17, 17))
    l100.universes[tube_x, tube_y] = tube_cold

    l101 = openmc.RectLattice(
        name='Fuel assembly (upper half)', lattice_id=101)
    l101.lower_left = (-10.71, -10.71)
    l101.pitch = (1.26, 1.26)
    l101.universes = np.tile(fuel_hot, (17, 17))
    l101.universes[tube_x, tube_y] = tube_hot

    # Define assemblies.
    fa_cw = openmc.Universe(name='Water assembly (cold)', universe_id=5)
    c50 = openmc.Cell(cell_id=50, fill=cold_water, region=+s34 & -s35)
    fa_cw.add_cell(c50)

    fa_hw = openmc.Universe(name='Water assembly (hot)', universe_id=7)
    c70 = openmc.Cell(cell_id=70, fill=hot_water, region=+s35 & -s36)
    fa_hw.add_cell(c70)

    fa_cold = openmc.Universe(name='Fuel assembly (cold)', universe_id=6)
    c60 = openmc.Cell(cell_id=60, fill=l100, region=+s34 & -s35)
    fa_cold.add_cell(c60)

    fa_hot = openmc.Universe(name='Fuel assembly (hot)', universe_id=8)
    c80 = openmc.Cell(cell_id=80, fill=l101, region=+s35 & -s36)
    fa_hot.add_cell(c80)

    # Define core lattices
    l200 = openmc.RectLattice(name='Core lattice (lower half)', lattice_id=200)
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

    l201 = openmc.RectLattice(name='Core lattice (lower half)', lattice_id=201)
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
    c1 = openmc.Cell(cell_id=1, fill=l200, region=-s6 & +s34 & -s35)
    c2 = openmc.Cell(cell_id=2, fill=l201, region=-s6 & +s35 & -s36)
    c3 = openmc.Cell(cell_id=3, fill=bot_plate, region=-s7 & +s31 & -s32)
    c4 = openmc.Cell(cell_id=4, fill=bot_nozzle, region=-s5 & +s32 & -s33)
    c5 = openmc.Cell(cell_id=5, fill=bot_fa, region=-s5 & +s33 & -s34)
    c6 = openmc.Cell(cell_id=6, fill=top_fa, region=-s5 & +s36 & -s37)
    c7 = openmc.Cell(cell_id=7, fill=top_nozzle, region=-s5 & +s37 & -s38)
    c8 = openmc.Cell(cell_id=8, fill=upper_rad_ref, region=-s7 & +s38 & -s39)
    c9 = openmc.Cell(cell_id=9, fill=bot_nozzle,
                     region=+s6 & -s7 & +s32 & -s38)
    c10 = openmc.Cell(cell_id=10, fill=rpv_steel,
                      region=+s7 & -s8 & +s31 & -s39)
    c11 = openmc.Cell(cell_id=11, fill=lower_rad_ref,
                      region=+s5 & -s6 & +s32 & -s34)
    c12 = openmc.Cell(cell_id=12, fill=upper_rad_ref,
                      region=+s5 & -s6 & +s36 & -s38)
    root.add_cells((c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12))

    # Assign root universe to geometry
    model.geometry.root_universe = root

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 100
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Box(
        [-160, -160, -183], [160, 160, 183]))

    plot = openmc.Plot()
    plot.origin = (125, 125, 0)
    plot.width = (250, 250)
    plot.pixels = (3000, 3000)
    plot.color_by = 'material'
    model.plots.append(plot)

    return model


def pwr_assembly():
    """Create a PWR assembly model.

    This model is a reflected 17x17 fuel assembly from the the `BEAVRS
    <http://crpg.mit.edu/research/beavrs>`_ benchmark. The fuel is 2.4 w/o
    enriched UO2 corresponding to a beginning-of-cycle condition. Note that the
    number of particles/batches is initially set very low for testing purposes.

    Returns
    -------
    model : openmc.model.Model
        A PWR assembly model

    """

    model = openmc.model.Model()

    # Define materials.
    fuel = openmc.Material(name='Fuel')
    fuel.set_density('g/cm3', 10.29769)
    fuel.add_nuclide('U234', 4.4843e-6)
    fuel.add_nuclide('U235', 5.5815e-4)
    fuel.add_nuclide('U238', 2.2408e-2)
    fuel.add_nuclide('O16', 4.5829e-2)

    clad = openmc.Material(name='Cladding')
    clad.set_density('g/cm3', 6.55)
    clad.add_nuclide('Zr90', 2.1827e-2)
    clad.add_nuclide('Zr91', 4.7600e-3)
    clad.add_nuclide('Zr92', 7.2758e-3)
    clad.add_nuclide('Zr94', 7.3734e-3)
    clad.add_nuclide('Zr96', 1.1879e-3)

    hot_water = openmc.Material(name='Hot borated water')
    hot_water.set_density('g/cm3', 0.740582)
    hot_water.add_nuclide('H1', 4.9457e-2)
    hot_water.add_nuclide('O16', 2.4672e-2)
    hot_water.add_nuclide('B10', 8.0042e-6)
    hot_water.add_nuclide('B11', 3.2218e-5)
    hot_water.add_s_alpha_beta('c_H_in_H2O')

    # Define the materials file.
    model.materials = (fuel, clad, hot_water)

    # Instantiate ZCylinder surfaces
    fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.39218, name='Fuel OR')
    clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.45720, name='Clad OR')

    # Create boundary planes to surround the geometry
    pitch = 21.42
    min_x = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
    max_x = openmc.XPlane(x0=+pitch/2, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
    max_y = openmc.YPlane(y0=+pitch/2, boundary_type='reflective')

    # Create a fuel pin universe
    fuel_pin_universe = openmc.Universe(name='Fuel Pin')
    fuel_cell = openmc.Cell(name='fuel', fill=fuel, region=-fuel_or)
    clad_cell = openmc.Cell(name='clad', fill=clad, region=+fuel_or & -clad_or)
    hot_water_cell = openmc.Cell(
        name='hot water', fill=hot_water, region=+clad_or)
    fuel_pin_universe.add_cells([fuel_cell, clad_cell, hot_water_cell])

    # Create a control rod guide tube universe
    guide_tube_universe = openmc.Universe(name='Guide Tube')
    gt_inner_cell = openmc.Cell(name='guide tube inner water', fill=hot_water,
                                region=-fuel_or)
    gt_clad_cell = openmc.Cell(name='guide tube clad', fill=clad,
                               region=+fuel_or & -clad_or)
    gt_outer_cell = openmc.Cell(name='guide tube outer water', fill=hot_water,
                                region=+clad_or)
    guide_tube_universe.add_cells([gt_inner_cell, gt_clad_cell, gt_outer_cell])

    # Create fuel assembly Lattice
    assembly = openmc.RectLattice(name='Fuel Assembly')
    assembly.pitch = (pitch/17, pitch/17)
    assembly.lower_left = (-pitch/2, -pitch/2)

    # Create array indices for guide tube locations in lattice
    template_x = np.array([5, 8, 11, 3, 13, 2, 5, 8, 11, 14, 2, 5, 8,
                           11, 14, 2, 5, 8, 11, 14, 3, 13, 5, 8, 11])
    template_y = np.array([2, 2, 2, 3, 3, 5, 5, 5, 5, 5, 8, 8, 8, 8,
                           8, 11, 11, 11, 11, 11, 13, 13, 14, 14, 14])

    # Create 17x17 array of universes
    assembly.universes = np.tile(fuel_pin_universe, (17, 17))
    assembly.universes[template_x, template_y] = guide_tube_universe

    # Create root Cell
    root_cell = openmc.Cell(name='root cell', fill=assembly)
    root_cell.region = +min_x & -max_x & +min_y & -max_y

    # Create root Universe
    model.geometry.root_universe = openmc.Universe(name='root universe')
    model.geometry.root_universe.add_cell(root_cell)

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 100
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box([-pitch/2, -pitch/2, -1],
                               [pitch/2, pitch/2, 1]),
        constraints={'fissionable': True}
    )

    plot = openmc.Plot()
    plot.origin = (0.0, 0.0, 0)
    plot.width = (21.42, 21.42)
    plot.pixels = (300, 300)
    plot.color_by = 'material'
    model.plots.append(plot)

    return model


def slab_mg(num_regions=1, mat_names=None, mgxslib_name='2g.h5'):
    """Create a 1D slab model.

    Parameters
    ----------
    num_regions : int, optional
        Number of regions in the problem, each with a unique MGXS dataset.
        Defaults to 1.

    mat_names : Iterable of str, optional
        List of the material names to use; defaults to ['mat_1', 'mat_2',...].

    mgxslib_name : str, optional
        MGXS Library file to use; defaults to '2g.h5'.

    Returns
    -------
    model : openmc.model.Model
        One-group, 1D slab model

    """

    openmc.check_type('num_regions', num_regions, Integral)
    openmc.check_greater_than('num_regions', num_regions, 0)
    if mat_names is not None:
        openmc.check_length('mat_names', mat_names, num_regions)
        openmc.check_iterable_type('mat_names', mat_names, str)
    else:
        mat_names = []
        for i in range(num_regions):
            mat_names.append('mat_' + str(i + 1))

    # # Make Materials
    materials_file = openmc.Materials()
    macros = []
    mats = []
    for i in range(len(mat_names)):
        macros.append(openmc.Macroscopic('mat_' + str(i + 1)))
        mats.append(openmc.Material(name=mat_names[i]))
        mats[-1].set_density('macro', 1.0)
        mats[-1].add_macroscopic(macros[-1])

    materials_file += mats

    materials_file.cross_sections = mgxslib_name

    # # Make Geometry
    rad_outer = 929.45
    # Set a cell boundary to exist for every material above (exclude the 0)
    rads = np.linspace(0., rad_outer, len(mats) + 1, endpoint=True)[1:]

    # Instantiate Universe
    root = openmc.Universe(universe_id=0, name='root universe')
    cells = []

    surfs = []
    surfs.append(openmc.XPlane(x0=0., boundary_type='reflective'))
    for r, rad in enumerate(rads):
        if r == len(rads) - 1:
            surfs.append(openmc.XPlane(x0=rad, boundary_type='vacuum'))
        else:
            surfs.append(openmc.XPlane(x0=rad))

    # Instantiate Cells
    cells = []
    for c in range(len(surfs) - 1):
        cells.append(openmc.Cell())
        cells[-1].region = (+surfs[c] & -surfs[c + 1])
        cells[-1].fill = mats[c]

    # Register Cells with Universe
    root.add_cells(cells)

    # Instantiate a Geometry, register the root Universe, and export to XML
    geometry_file = openmc.Geometry(root)

    # # Make Settings
    # Instantiate a Settings object, set all runtime parameters
    settings_file = openmc.Settings()
    settings_file.energy_mode = 'multi-group'
    settings_file.tabular_legendre = {'enable': False}
    settings_file.batches = 10
    settings_file.inactive = 5
    settings_file.particles = 1000

    # Build source distribution
    INF = 1000.
    bounds = [0., -INF, -INF, rads[0], INF, INF]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    settings_file.source = openmc.IndependentSource(space=uniform_dist)

    settings_file.output = {'summary': False}

    model = openmc.model.Model()
    model.geometry = geometry_file
    model.materials = materials_file
    model.settings = settings_file
    model.xs_data = macros

    return model


def random_ray_lattice():
    """Create a 2x2 PWR pincell asymmetrical lattic eexample.

    This model is a 2x2 reflective lattice of fuel pins with one of the lattice
    locations having just moderator instead of a fuel pin. It uses 7 group
    cross section data.

    Returns
    -------
    model : openmc.model.Model
        A PWR 2x2 lattice model

    """
    model = openmc.model.Model()

    ###########################################################################
    # Create MGXS data for the problem

    # Instantiate the energy group data
    group_edges = [1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges)

    # Instantiate the 7-group (C5G7) cross section data
    uo2_xsdata = openmc.XSdata('UO2', groups)
    uo2_xsdata.order = 0
    uo2_xsdata.set_total(
        [0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678,
         0.5644058])
    uo2_xsdata.set_absorption([8.0248e-03, 3.7174e-03, 2.6769e-02, 9.6236e-02,
                               3.0020e-02, 1.1126e-01, 2.8278e-01])
    scatter_matrix = np.array(
        [[[0.1275370, 0.0423780, 0.0000094, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.3244560, 0.0016314, 0.0000000,
              0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.4509400, 0.0026792,
              0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.4525650,
              0.0055664, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.0001253,
              0.2714010, 0.0102550, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000,
              0.0012968, 0.2658020, 0.0168090],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0085458, 0.2730800]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    uo2_xsdata.set_scatter_matrix(scatter_matrix)
    uo2_xsdata.set_fission([7.21206e-03, 8.19301e-04, 6.45320e-03,
                            1.85648e-02, 1.78084e-02, 8.30348e-02,
                            2.16004e-01])
    uo2_xsdata.set_nu_fission([2.005998e-02, 2.027303e-03, 1.570599e-02,
                               4.518301e-02, 4.334208e-02, 2.020901e-01,
                               5.257105e-01])
    uo2_xsdata.set_chi([5.8791e-01, 4.1176e-01, 3.3906e-04, 1.1761e-07, 0.0000e+00,
                        0.0000e+00, 0.0000e+00])

    h2o_xsdata = openmc.XSdata('LWTR', groups)
    h2o_xsdata.order = 0
    h2o_xsdata.set_total([0.15920605, 0.412969593, 0.59030986, 0.58435,
                          0.718, 1.2544497, 2.650379])
    h2o_xsdata.set_absorption([6.0105e-04, 1.5793e-05, 3.3716e-04,
                               1.9406e-03, 5.7416e-03, 1.5001e-02,
                               3.7239e-02])
    scatter_matrix = np.array(
        [[[0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0000000, 0.0000000],
          [0.0000000, 0.2823340, 0.1299400, 0.0006234,
              0.0000480, 0.0000074, 0.0000010],
          [0.0000000, 0.0000000, 0.3452560, 0.2245700,
              0.0169990, 0.0026443, 0.0005034],
          [0.0000000, 0.0000000, 0.0000000, 0.0910284,
              0.4155100, 0.0637320, 0.0121390],
          [0.0000000, 0.0000000, 0.0000000, 0.0000714,
              0.1391380, 0.5118200, 0.0612290],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000,
              0.0022157, 0.6999130, 0.5373200],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1324400, 2.4807000]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    h2o_xsdata.set_scatter_matrix(scatter_matrix)

    mg_cross_sections = openmc.MGXSLibrary(groups)
    mg_cross_sections.add_xsdatas([uo2_xsdata, h2o_xsdata])
    mg_cross_sections.export_to_hdf5('mgxs.h5')

    ###########################################################################
    # Create materials for the problem

    # Instantiate some Materials and register the appropriate macroscopic data
    uo2 = openmc.Material(name='UO2 fuel')
    uo2.set_density('macro', 1.0)
    uo2.add_macroscopic('UO2')

    water = openmc.Material(name='Water')
    water.set_density('macro', 1.0)
    water.add_macroscopic('LWTR')

    # Instantiate a Materials collection and export to XML
    materials = openmc.Materials([uo2, water])
    materials.cross_sections = "mgxs.h5"

    ###########################################################################
    # Define problem geometry

    ########################################
    # Define an unbounded pincell universe

    pitch = 1.26

    # Create a surface for the fuel outer radius
    fuel_or = openmc.ZCylinder(r=0.54, name='Fuel OR')
    inner_ring_a = openmc.ZCylinder(r=0.33, name='inner ring a')
    inner_ring_b = openmc.ZCylinder(r=0.45, name='inner ring b')
    outer_ring_a = openmc.ZCylinder(r=0.60, name='outer ring a')
    outer_ring_b = openmc.ZCylinder(r=0.69, name='outer ring b')

    # Instantiate Cells
    fuel_a = openmc.Cell(fill=uo2, region=-inner_ring_a, name='fuel inner a')
    fuel_b = openmc.Cell(fill=uo2, region=+inner_ring_a & -
                         inner_ring_b, name='fuel inner b')
    fuel_c = openmc.Cell(fill=uo2, region=+inner_ring_b & -
                         fuel_or, name='fuel inner c')
    moderator_a = openmc.Cell(
        fill=water, region=+fuel_or & -outer_ring_a, name='moderator inner a')
    moderator_b = openmc.Cell(
        fill=water, region=+outer_ring_a & -outer_ring_b, name='moderator outer b')
    moderator_c = openmc.Cell(
        fill=water, region=+outer_ring_b, name='moderator outer c')

    # Create pincell universe
    pincell_base = openmc.Universe()

    # Register Cells with Universe
    pincell_base.add_cells(
        [fuel_a, fuel_b, fuel_c, moderator_a, moderator_b, moderator_c])

    # Create planes for azimuthal sectors
    azimuthal_planes = []
    for i in range(8):
        angle = 2 * i * openmc.pi / 8
        normal_vector = (-openmc.sin(angle), openmc.cos(angle), 0)
        azimuthal_planes.append(openmc.Plane(
            a=normal_vector[0], b=normal_vector[1], c=normal_vector[2], d=0))

    # Create a cell for each azimuthal sector
    azimuthal_cells = []
    for i in range(8):
        azimuthal_cell = openmc.Cell(name=f'azimuthal_cell_{i}')
        azimuthal_cell.fill = pincell_base
        azimuthal_cell.region = + \
            azimuthal_planes[i] & -azimuthal_planes[(i+1) % 8]
        azimuthal_cells.append(azimuthal_cell)

    # Create a geometry with the azimuthal universes
    pincell = openmc.Universe(cells=azimuthal_cells)

    ########################################
    # Define a moderator lattice universe

    moderator_infinite = openmc.Cell(fill=water, name='moderator infinite')
    mu = openmc.Universe()
    mu.add_cells([moderator_infinite])

    lattice = openmc.RectLattice()
    lattice.lower_left = [-pitch/2.0, -pitch/2.0]
    lattice.pitch = [pitch/10.0, pitch/10.0]
    lattice.universes = np.full((10, 10), mu)

    mod_lattice_cell = openmc.Cell(fill=lattice)

    mod_lattice_uni = openmc.Universe()

    mod_lattice_uni.add_cells([mod_lattice_cell])

    ########################################
    # Define 2x2 outer lattice
    lattice2x2 = openmc.RectLattice()
    lattice2x2.lower_left = (-pitch, -pitch)
    lattice2x2.pitch = (pitch, pitch)
    lattice2x2.universes = [
        [pincell, pincell],
        [pincell, mod_lattice_uni]
    ]

    ########################################
    # Define cell containing lattice and other stuff
    box = openmc.model.RectangularPrism(
        pitch*2, pitch*2, boundary_type='reflective')

    assembly = openmc.Cell(fill=lattice2x2, region=-box, name='assembly')

    # Create a geometry with the top-level cell
    geometry = openmc.Geometry([assembly])

    ###########################################################################
    # Define problem settings

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.batches = 10
    settings.inactive = 5
    settings.particles = 100

    # Create an initial uniform spatial source distribution over fissionable zones
    lower_left = (-pitch, -pitch, -1)
    upper_right = (pitch, pitch, 1)
    uniform_dist = openmc.stats.Box(lower_left, upper_right)
    rr_source = openmc.IndependentSource(space=uniform_dist)

    settings.random_ray['distance_active'] = 100.0
    settings.random_ray['distance_inactive'] = 20.0
    settings.random_ray['ray_source'] = rr_source
    settings.random_ray['volume_normalized_flux_tallies'] = True

    ###########################################################################
    # Define tallies

    # Create a mesh that will be used for tallying
    mesh = openmc.RegularMesh()
    mesh.dimension = (2, 2)
    mesh.lower_left = (-pitch, -pitch)
    mesh.upper_right = (pitch, pitch)

    # Create a mesh filter that can be used in a tally
    mesh_filter = openmc.MeshFilter(mesh)

    # Create an energy group filter as well
    group_edges = [1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6]
    energy_filter = openmc.EnergyFilter(group_edges)

    # Now use the mesh filter in a tally and indicate what scores are desired
    tally = openmc.Tally(name="Mesh tally")
    tally.filters = [mesh_filter, energy_filter]
    tally.scores = ['flux', 'fission', 'nu-fission']
    tally.estimator = 'analog'

    # Instantiate a Tallies collection and export to XML
    tallies = openmc.Tallies([tally])

    ###########################################################################
    #                   Exporting to OpenMC model
    ###########################################################################

    model.geometry = geometry
    model.materials = materials
    model.settings = settings
    model.tallies = tallies
    return model


def random_ray_three_region_cube():
    """Create a three region cube model.

    This is a simple monoenergetic problem of a cube with three concentric cubic
    regions. The innermost region is near void (with Sigma_t around 10^-5) and
    contains an external isotropic source term, the middle region is void (with
    Sigma_t around 10^-4), and the outer region of the cube is an absorber
    (with Sigma_t around 1).

    Returns
    -------
    model : openmc.model.Model
        A three region cube model

    """

    model = openmc.model.Model()

    ###########################################################################
    # Helper function creates a 3 region cube with different fills in each region
    def fill_cube(N, n_1, n_2, fill_1, fill_2, fill_3):
        cube = [[[0 for _ in range(N)] for _ in range(N)] for _ in range(N)]
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    if i < n_1 and j >= (N-n_1) and k < n_1:
                        cube[i][j][k] = fill_1
                    elif i < n_2 and j >= (N-n_2) and k < n_2:
                        cube[i][j][k] = fill_2
                    else:
                        cube[i][j][k] = fill_3
        return cube

    ###########################################################################
    # Create multigroup data

    # Instantiate the energy group data
    ebins = [1e-5, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges=ebins)

    void_sigma_a = 4.0e-6
    void_sigma_s = 3.0e-4
    void_mat_data = openmc.XSdata('void', groups)
    void_mat_data.order = 0
    void_mat_data.set_total([void_sigma_a + void_sigma_s])
    void_mat_data.set_absorption([void_sigma_a])
    void_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[void_sigma_s]]]), 0, 3))

    absorber_sigma_a = 0.75
    absorber_sigma_s = 0.25
    absorber_mat_data = openmc.XSdata('absorber', groups)
    absorber_mat_data.order = 0
    absorber_mat_data.set_total([absorber_sigma_a + absorber_sigma_s])
    absorber_mat_data.set_absorption([absorber_sigma_a])
    absorber_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[absorber_sigma_s]]]), 0, 3))

    multiplier = 0.1
    source_sigma_a = void_sigma_a * multiplier
    source_sigma_s = void_sigma_s * multiplier
    source_mat_data = openmc.XSdata('source', groups)
    source_mat_data.order = 0
    source_mat_data.set_total([source_sigma_a + source_sigma_s])
    source_mat_data.set_absorption([source_sigma_a])
    source_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[source_sigma_s]]]), 0, 3))

    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas(
        [source_mat_data, void_mat_data, absorber_mat_data])
    mg_cross_sections_file.export_to_hdf5()

    ###########################################################################
    # Create materials for the problem

    # Instantiate some Macroscopic Data
    source_data = openmc.Macroscopic('source')
    void_data = openmc.Macroscopic('void')
    absorber_data = openmc.Macroscopic('absorber')

    # Instantiate some Materials and register the appropriate Macroscopic objects
    source_mat = openmc.Material(name='source')
    source_mat.set_density('macro', 1.0)
    source_mat.add_macroscopic(source_data)

    void_mat = openmc.Material(name='void')
    void_mat.set_density('macro', 1.0)
    void_mat.add_macroscopic(void_data)

    absorber_mat = openmc.Material(name='absorber')
    absorber_mat.set_density('macro', 1.0)
    absorber_mat.add_macroscopic(absorber_data)

    # Instantiate a Materials collection and export to XML
    materials_file = openmc.Materials([source_mat, void_mat, absorber_mat])
    materials_file.cross_sections = "mgxs.h5"

    ###########################################################################
    # Define problem geometry

    source_cell = openmc.Cell(fill=source_mat, name='infinite source region')
    void_cell = openmc.Cell(fill=void_mat, name='infinite void region')
    absorber_cell = openmc.Cell(
        fill=absorber_mat, name='infinite absorber region')

    source_universe = openmc.Universe(name='source universe')
    source_universe.add_cells([source_cell])

    void_universe = openmc.Universe()
    void_universe.add_cells([void_cell])

    absorber_universe = openmc.Universe()
    absorber_universe.add_cells([absorber_cell])

    absorber_width = 30.0
    n_base = 6

    # This variable can be increased above 1 to refine the FSR mesh resolution further
    refinement_level = 2

    n = n_base * refinement_level
    pitch = absorber_width / n

    pattern = fill_cube(n, 1*refinement_level, 5*refinement_level,
                        source_universe, void_universe, absorber_universe)

    lattice = openmc.RectLattice()
    lattice.lower_left = [0.0, 0.0, 0.0]
    lattice.pitch = [pitch, pitch, pitch]
    lattice.universes = pattern

    lattice_cell = openmc.Cell(fill=lattice)

    lattice_uni = openmc.Universe()
    lattice_uni.add_cells([lattice_cell])

    x_low = openmc.XPlane(x0=0.0, boundary_type='reflective')
    x_high = openmc.XPlane(x0=absorber_width, boundary_type='vacuum')

    y_low = openmc.YPlane(y0=0.0, boundary_type='reflective')
    y_high = openmc.YPlane(y0=absorber_width, boundary_type='vacuum')

    z_low = openmc.ZPlane(z0=0.0, boundary_type='reflective')
    z_high = openmc.ZPlane(z0=absorber_width, boundary_type='vacuum')

    full_domain = openmc.Cell(fill=lattice_uni, region=+x_low & -
                              x_high & +y_low & -y_high & +z_low & -z_high, name='full domain')

    root = openmc.Universe(name='root universe')
    root.add_cell(full_domain)

    # Create a geometry with the two cells and export to XML
    geometry = openmc.Geometry(root)

    ###########################################################################
    # Define problem settings

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.inactive = 5
    settings.batches = 10
    settings.particles = 90
    settings.run_mode = 'fixed source'

    # Create an initial uniform spatial source for ray integration
    lower_left_ray = [0.0, 0.0, 0.0]
    upper_right_ray = [absorber_width, absorber_width, absorber_width]
    uniform_dist_ray = openmc.stats.Box(
        lower_left_ray, upper_right_ray, only_fissionable=False)
    rr_source = openmc.IndependentSource(space=uniform_dist_ray)

    settings.random_ray['distance_active'] = 500.0
    settings.random_ray['distance_inactive'] = 100.0
    settings.random_ray['ray_source'] = rr_source
    settings.random_ray['volume_normalized_flux_tallies'] = True

    # Create the neutron source in the bottom right of the moderator
    # Good - fast group appears largest (besides most thermal)
    strengths = [1.0]
    midpoints = [100.0]
    energy_distribution = openmc.stats.Discrete(x=midpoints, p=strengths)

    source = openmc.IndependentSource(energy=energy_distribution, constraints={
                                      'domains': [source_universe]}, strength=3.14)

    settings.source = [source]

    ###########################################################################
    # Define tallies

    estimator = 'tracklength'

    absorber_filter = openmc.MaterialFilter(absorber_mat)
    absorber_tally = openmc.Tally(name="Absorber Tally")
    absorber_tally.filters = [absorber_filter]
    absorber_tally.scores = ['flux']
    absorber_tally.estimator = estimator

    void_filter = openmc.MaterialFilter(void_mat)
    void_tally = openmc.Tally(name="Void Tally")
    void_tally.filters = [void_filter]
    void_tally.scores = ['flux']
    void_tally.estimator = estimator

    source_filter = openmc.MaterialFilter(source_mat)
    source_tally = openmc.Tally(name="Source Tally")
    source_tally.filters = [source_filter]
    source_tally.scores = ['flux']
    source_tally.estimator = estimator

    # Instantiate a Tallies collection and export to XML
    tallies = openmc.Tallies([source_tally, void_tally, absorber_tally])

    ###########################################################################
    # Assmble Model

    model.geometry = geometry
    model.materials = materials_file
    model.settings = settings
    model.tallies = tallies

    return model
