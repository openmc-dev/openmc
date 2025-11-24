"""
HEU-SOL-THERM-038-1:  WINCO SLAB TANKS; 3.692" Gap
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Tank #1 Solution
mat11 = openmc.Material(material_id=11)
mat11.set_density("sum")
mat11.add_nuclide("U234", 8.747700e-06)
mat11.add_nuclide("U235", 9.633800e-04)
mat11.add_nuclide("U236", 2.874700e-06)
mat11.add_nuclide("U238", 5.902700e-05)
mat11.add_nuclide("H1", 5.784200e-02)
mat11.add_element("N", 2.260800e-03)
mat11.add_nuclide("O16", 3.767500e-02)
mat11.add_s_alpha_beta("c_H_in_H2O")

# Tank #2 Solution
mat12 = openmc.Material(material_id=12)
mat12.set_density("sum")
mat12.add_nuclide("U234", 8.767300e-06)
mat12.add_nuclide("U235", 9.640400e-04)
mat12.add_nuclide("U236", 2.878700e-06)
mat12.add_nuclide("U238", 5.906900e-05)
mat12.add_nuclide("H1", 5.787100e-02)
mat12.add_element("N", 2.262200e-03)
mat12.add_nuclide("O16", 3.769500e-02)
mat12.add_s_alpha_beta("c_H_in_H2O")

# Hoop #1 SST
mat21 = openmc.Material(material_id=21)
mat21.set_density("sum")
mat21.add_element("P", 7.889400e+00)
mat21.add_element("Fe", 7.019400e+01)
mat21.add_element("Cr", 1.830000e+01)
mat21.add_element("Ni", 8.210000e+00)
mat21.add_element("Mn", 1.910000e+00)
mat21.add_element("Si", 6.700000e-01)
mat21.add_element("Cu", 3.100000e-01)
mat21.add_element("Mo", 2.600000e-01)
mat21.add_element("N", 6.000000e-02)
mat21.add_element("P", 3.100000e-02)
mat21.add_element("C", 5.000000e-02)
mat21.add_element("S", 5.000000e-03)

# Hoop #2 SST
mat22 = openmc.Material(material_id=22)
mat22.set_density("sum")
mat22.add_element("P", 7.893200e+00)
mat22.add_element("Fe", 7.019400e+01)
mat22.add_element("Cr", 1.830000e+01)
mat22.add_element("Ni", 8.210000e+00)
mat22.add_element("Mn", 1.910000e+00)
mat22.add_element("Si", 6.700000e-01)
mat22.add_element("Cu", 3.100000e-01)
mat22.add_element("Mo", 2.600000e-01)
mat22.add_element("N", 6.000000e-02)
mat22.add_element("P", 3.100000e-02)
mat22.add_element("C", 5.000000e-02)
mat22.add_element("S", 5.000000e-03)

# Tank Plates SST
mat23 = openmc.Material(material_id=23)
mat23.set_density("sum")
mat23.add_element("P", 7.883000e+00)
mat23.add_element("Fe", 7.019400e+01)
mat23.add_element("Cr", 1.830000e+01)
mat23.add_element("Ni", 8.210000e+00)
mat23.add_element("Mn", 1.910000e+00)
mat23.add_element("Si", 6.700000e-01)
mat23.add_element("Cu", 3.100000e-01)
mat23.add_element("Mo", 2.600000e-01)
mat23.add_element("N", 6.000000e-02)
mat23.add_element("P", 3.100000e-02)
mat23.add_element("C", 5.000000e-02)
mat23.add_element("S", 5.000000e-03)

# Center Support SST
mat24 = openmc.Material(material_id=24)
mat24.set_density("sum")
mat24.add_element("P", 7.829700e+00)
mat24.add_element("Fe", 7.019400e+01)
mat24.add_element("Cr", 1.830000e+01)
mat24.add_element("Ni", 8.210000e+00)
mat24.add_element("Mn", 1.910000e+00)
mat24.add_element("Si", 6.700000e-01)
mat24.add_element("Cu", 3.100000e-01)
mat24.add_element("Mo", 2.600000e-01)
mat24.add_element("N", 6.000000e-02)
mat24.add_element("P", 3.100000e-02)
mat24.add_element("C", 5.000000e-02)
mat24.add_element("S", 5.000000e-03)

# Al-6061 Supports
mat30 = openmc.Material(material_id=30)
mat30.set_density("sum")
mat30.add_element("P", 2.690000e+00)
mat30.add_element("Al", 9.732500e+01)
mat30.add_element("Mg", 1.000000e+00)
mat30.add_element("Si", 6.000000e-01)
mat30.add_element("Fe", 3.500000e-01)
mat30.add_element("Cu", 2.500000e-01)
mat30.add_element("Cr", 2.000000e-01)
mat30.add_element("Zn", 1.250000e-01)
mat30.add_element("Mn", 7.500000e-02)
mat30.add_element("Ti", 7.500000e-02)

materials = openmc.Materials([mat11, mat12, mat21, mat22, mat23, mat24, mat30])

# ==============================================================================
# Geometry
# ==============================================================================

# Support Post
surf1 = openmc.ZCylinder(surface_id=1, r=1.26873)
# Hoop #1/Inner
surf2 = openmc.ZCylinder(surface_id=2, r=35.81654)
# Tank #1/Outer
surf3 = openmc.ZCylinder(surface_id=3, r=37.90696)
# Plate C/Upper
# surf4: Unsupported surface type "analytic" with params ['1.', 'z', '-0.967994', 'constant']
# Plate B/Lower
# surf5: Unsupported surface type "analytic" with params ['1.', 'z', '-9.868664', 'constant']
# Part A
surf11 = openmc.model.RectangularParallelepiped(-39.217600000000004, -38.5826, -38.5826, 38.5826, 0.6700000000000004, 8.290000000000001)
# Part B
surf12 = openmc.model.RectangularParallelepiped(38.5826, 39.217600000000004, -38.5826, 38.5826, 0.6700000000000004, 8.290000000000001)
# Part C
surf13 = openmc.model.RectangularParallelepiped(-58.42, 58.42, -39.217600000000004, -38.5826, 0.6700000000000004, 8.290000000000001)
# Part D
surf14 = openmc.model.RectangularParallelepiped(-58.42, 58.42, 38.5826, 39.217600000000004, 0.6700000000000004, 8.290000000000001)
# Part E
surf15 = openmc.model.RectangularParallelepiped(48.260000000000005, 58.42, -58.42, 58.42, -1.87, 0.67)
# Part F
surf16 = openmc.model.RectangularParallelepiped(-58.42, -48.260000000000005, -58.42, 58.42, -1.87, 0.67)
# Hoop #2/Inner
surf22 = openmc.ZCylinder(surface_id=22, r=35.78733)
# Tank #2/Outer
surf23 = openmc.ZCylinder(surface_id=23, x0=-10.84199, y0=0.0, r=37.95903)
# Plate B/Upper
# surf24: Unsupported surface type "analytic" with params ['1.', 'z', '+9.885426', 'constant', 'tr', '0.', '0.', '-9.37768']
# Plate A/Lower
# surf25: Unsupported surface type "analytic" with params ['1.', 'z', '+0.973836', 'constant', 'tr', '0.', '0.', '-9.37768']
surf31 = openmc.ZCylinder(surface_id=31, x0=tr, y0=-26.4922, r=1.905)
surf32 = openmc.ZCylinder(surface_id=32, x0=-46.40199, y0=-10.84199, r=2.54)
surf33 = openmc.ZCylinder(surface_id=33, x0=tr, y0=+26.4922, r=1.905)
surf34 = openmc.ZCylinder(surface_id=34, x0=-46.40199, y0=-10.84199, r=2.54)
surf35 = openmc.ZCylinder(surface_id=35, x0=tr, y0=-26.4922, r=1.905)
surf36 = openmc.ZCylinder(surface_id=36, x0=-46.40199, y0=-10.84199, r=2.54)
surf37 = openmc.ZCylinder(surface_id=37, x0=tr, y0=+26.4922, r=1.905)
surf38 = openmc.ZCylinder(surface_id=38, x0=-46.40199, y0=-10.84199, r=2.54)
# Support Plate G with
surf40 = openmc.model.RectangularParallelepiped(-39.37, 39.37, -39.37, 39.37, -57.04967, -55.77967)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Post
cell1 = openmc.Cell(cell_id=1, fill=mat24)
cell1.region = -surf1 & +surf4 & -surf5

# Soln--1
cell2 = openmc.Cell(cell_id=2, fill=mat11)
cell2.region = +surf1 & -surf2 & +surf4 & -surf5

# Hoop--1
cell3 = openmc.Cell(cell_id=3, fill=mat21)
cell3.region = +surf2 & -surf3 & +surf4 & -surf5

# Plate-B
cell4 = openmc.Cell(cell_id=4, fill=mat23)
cell4.region = -surf3 & +surf5

# Plate-C
cell5 = openmc.Cell(cell_id=5, fill=mat23)
cell5.region = -surf3 & -surf4

# Part--A
cell6 = openmc.Cell(cell_id=6, fill=mat30)
cell6.region = -surf11

# Part--B
cell7 = openmc.Cell(cell_id=7, fill=mat30)
cell7.region = -surf12

# Part--C
cell8 = openmc.Cell(cell_id=8, fill=mat30)
cell8.region = -surf13

# Part--D
cell9 = openmc.Cell(cell_id=9, fill=mat30)
cell9.region = -surf14

# Part--E
cell10 = openmc.Cell(cell_id=10, fill=mat30)
cell10.region = -surf15

# Part--F
cell11 = openmc.Cell(cell_id=11, fill=mat30)
cell11.region = -surf16

# Post
cell12 = openmc.Cell(cell_id=12, fill=mat24)
cell12.region = -surf1 & +surf24 & -surf25

# Soln--2
cell13 = openmc.Cell(cell_id=13, fill=mat12)
cell13.region = +surf1 & -surf22 & +surf24 & -surf25

# Hoop--2
cell14 = openmc.Cell(cell_id=14, fill=mat22)
cell14.region = +surf22 & -surf23 & +surf24 & -surf25

# Plate-A
cell15 = openmc.Cell(cell_id=15, fill=mat23)
cell15.region = -surf23 & +surf25

# Plate-B
cell16 = openmc.Cell(cell_id=16, fill=mat23)
cell16.region = -surf23 & -surf24

# Air1
cell17 = openmc.Cell(cell_id=17)
cell17.region = +surf23 & -surf31 & -surf32 & +surf40

# Leg1
cell18 = openmc.Cell(cell_id=18, fill=mat30)
cell18.region = +surf23 & +surf31 & -surf32 & +surf40

# Air2
cell19 = openmc.Cell(cell_id=19)
cell19.region = +surf23 & -surf33 & -surf34 & +surf40

# Leg2
cell20 = openmc.Cell(cell_id=20, fill=mat30)
cell20.region = +surf23 & +surf33 & -surf34 & +surf40

# Air3
cell21 = openmc.Cell(cell_id=21)
cell21.region = +surf23 & -surf35 & -surf36 & +surf40

# Leg3
cell22 = openmc.Cell(cell_id=22, fill=mat30)
cell22.region = +surf23 & +surf35 & -surf36 & +surf40

# Air4
cell23 = openmc.Cell(cell_id=23)
cell23.region = +surf23 & -surf37 & -surf38 & +surf40

# Leg4
cell24 = openmc.Cell(cell_id=24, fill=mat30)
cell24.region = +surf23 & +surf37 & -surf38 & +surf40

# Plate-G
cell25 = openmc.Cell(cell_id=25, fill=mat30)
cell25.region = -surf40

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24, cell25])
geometry = openmc.Geometry(root_universe)

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 15000
settings.batches = 4400
settings.inactive = 100
settings.run_mode = "eigenvalue"

source = openmc.IndependentSource()
source.space = openmc.stats.Box((-4.0, -4.0, -14.0), (4.0, 4.0, 4.0))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
