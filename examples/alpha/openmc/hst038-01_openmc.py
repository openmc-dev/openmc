"""
HEU-SOL-THERM-038-1:  WINCO SLAB TANKS; 3.692" Gap
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Tank #1 Solution
mat11 = openmc.Material(material_id=11, name="Tank #1 Solution")
mat11.set_density("atom/b-cm", 1.034029e-03)
mat11.add_nuclide("U234", 8.7477e-6)
mat11.add_nuclide("U235", 9.6338e-4)
mat11.add_nuclide("U236", 2.8747e-6)
mat11.add_nuclide("U238", 5.9027e-5)

# Tank #2 Solution
mat12 = openmc.Material(material_id=12, name="Tank #2 Solution")
mat12.set_density("atom/b-cm", 1.034755e-03)
mat12.add_nuclide("U234", 8.7673e-6)
mat12.add_nuclide("U235", 9.6404e-4)
mat12.add_nuclide("U236", 2.8787e-6)
mat12.add_nuclide("U238", 5.9069e-5)

# Hoop #1 SST
mat21 = openmc.Material(material_id=21, name="Hoop #1 SST")
mat21.set_density("atom/b-cm", 1.074834e+02)
mat21.add_element("P", 7.8894)
mat21.add_element("Fe", 70.194)
mat21.add_element("Cr", 18.30)
mat21.add_element("Ni", 8.21)
mat21.add_element("Mn", 1.91)
mat21.add_element("Si", 0.67)
mat21.add_element("Cu", 0.31)

# Hoop #2 SST
mat22 = openmc.Material(material_id=22, name="Hoop #2 SST")
mat22.set_density("atom/b-cm", 1.074872e+02)
mat22.add_element("P", 7.8932)
mat22.add_element("Fe", 70.194)
mat22.add_element("Cr", 18.30)
mat22.add_element("Ni", 8.21)
mat22.add_element("Mn", 1.91)
mat22.add_element("Si", 0.67)
mat22.add_element("Cu", 0.31)

# Tank Plates SST
mat23 = openmc.Material(material_id=23, name="Tank Plates SST")
mat23.set_density("atom/b-cm", 1.074770e+02)
mat23.add_element("P", 7.8830)
mat23.add_element("Fe", 70.194)
mat23.add_element("Cr", 18.30)
mat23.add_element("Ni", 8.21)
mat23.add_element("Mn", 1.91)
mat23.add_element("Si", 0.67)
mat23.add_element("Cu", 0.31)

# Center Support SST
mat24 = openmc.Material(material_id=24, name="Center Support SST")
mat24.set_density("atom/b-cm", 1.074237e+02)
mat24.add_element("P", 7.8297)
mat24.add_element("Fe", 70.194)
mat24.add_element("Cr", 18.30)
mat24.add_element("Ni", 8.21)
mat24.add_element("Mn", 1.91)
mat24.add_element("Si", 0.67)
mat24.add_element("Cu", 0.31)

# Al-6061 Supports
mat30 = openmc.Material(material_id=30, name="Al-6061 Supports")
mat30.set_density("atom/b-cm", 1.024150e+02)
mat30.add_element("P", 2.69)
mat30.add_element("Al", 97.325)
mat30.add_element("Mg", 1.0)
mat30.add_element("Si", 0.6)
mat30.add_element("Fe", 0.35)
mat30.add_element("Cu", 0.25)
mat30.add_element("Cr", 0.2)

materials = openmc.Materials([mat11, mat12, mat21, mat22, mat23, mat24, mat30])
materials.export_to_xml()

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
# COG surface type "analytic" with parameters: 1. z -0.967994 constant
# This surface type requires manual translation to OpenMC
surf4 = openmc.Sphere(surface_id=4, r=1.0)  # PLACEHOLDER - REPLACE THIS

# Plate B/Lower
# COG surface type "analytic" with parameters: 1. z -9.868664 constant
# This surface type requires manual translation to OpenMC
surf5 = openmc.Sphere(surface_id=5, r=1.0)  # PLACEHOLDER - REPLACE THIS

# Part A
surf11 = openmc.model.RectangularParallelepiped(
    -39.217600000000004, -38.5826, -38.5826, 38.5826, 0.6700000000000004, 8.290000000000001, surface_id=11)

# Part B
surf12 = openmc.model.RectangularParallelepiped(
    38.5826, 39.217600000000004, -38.5826, 38.5826, 0.6700000000000004, 8.290000000000001, surface_id=12)

# Part C
surf13 = openmc.model.RectangularParallelepiped(
    -58.42, 58.42, -39.217600000000004, -38.5826, 0.6700000000000004, 8.290000000000001, surface_id=13)

# Part D
surf14 = openmc.model.RectangularParallelepiped(
    -58.42, 58.42, 38.5826, 39.217600000000004, 0.6700000000000004, 8.290000000000001, surface_id=14)

# Part E
surf15 = openmc.model.RectangularParallelepiped(
    48.260000000000005, 58.42, -58.42, 58.42, -1.87, 0.67, surface_id=15)

# Part F
surf16 = openmc.model.RectangularParallelepiped(
    -58.42, -48.260000000000005, -58.42, 58.42, -1.87, 0.67, surface_id=16)

# Hoop #2/Inner
surf22 = openmc.ZCylinder(surface_id=22, r=35.78733)

# Tank #2/Outer
surf23 = openmc.ZCylinder(surface_id=23, r=37.95903)

# Plate B/Upper
# COG surface type "analytic" with parameters: 1. z +9.885426 constant tr 0. 0. -9.37768
# This surface type requires manual translation to OpenMC
surf24 = openmc.Sphere(surface_id=24, r=1.0)  # PLACEHOLDER - REPLACE THIS

# Plate A/Lower
# COG surface type "analytic" with parameters: 1. z +0.973836 constant tr 0. 0. -9.37768
# This surface type requires manual translation to OpenMC
surf25 = openmc.Sphere(surface_id=25, r=1.0)  # PLACEHOLDER - REPLACE THIS

surf31 = openmc.ZCylinder(surface_id=31, r=1.905)

surf32 = openmc.ZCylinder(surface_id=32, r=2.540)

surf33 = openmc.ZCylinder(surface_id=33, r=1.905)

surf34 = openmc.ZCylinder(surface_id=34, r=2.540)

surf35 = openmc.ZCylinder(surface_id=35, r=1.905)

surf36 = openmc.ZCylinder(surface_id=36, r=2.540)

surf37 = openmc.ZCylinder(surface_id=37, r=1.905)

surf38 = openmc.ZCylinder(surface_id=38, r=2.540)

# Support Plate G with
surf40 = openmc.model.RectangularParallelepiped(
    -39.37, 39.37, -39.37, 39.37, -57.04967, -55.77967, surface_id=40)


# Cell: Air1
cell0 = openmc.Cell(cell_id=0, fill=mat0, name="Air1")
cell0.region = +surf23 & -surf31 & -surf32 & +surf40

# Cell: Air2
cell1 = openmc.Cell(cell_id=1, fill=mat0, name="Air2")
cell1.region = +surf23 & -surf33 & -surf34 & +surf40

# Cell: Air3
cell2 = openmc.Cell(cell_id=2, fill=mat0, name="Air3")
cell2.region = +surf23 & -surf35 & -surf36 & +surf40

# Cell: Air4
cell3 = openmc.Cell(cell_id=3, fill=mat0, name="Air4")
cell3.region = +surf23 & -surf37 & -surf38 & +surf40

# Cell: Soln--1
cell4 = openmc.Cell(cell_id=4, fill=mat11, name="Soln--1")
cell4.region = +surf1 & -surf2 & +surf4 & -surf5

# Cell: Soln--2
cell5 = openmc.Cell(cell_id=5, fill=mat12, name="Soln--2")
cell5.region = +surf1 & -surf22 & +surf24 & -surf25

# Cell: Hoop--1
cell6 = openmc.Cell(cell_id=6, fill=mat21, name="Hoop--1")
cell6.region = +surf2 & -surf3 & +surf4 & -surf5

# Cell: Hoop--2
cell7 = openmc.Cell(cell_id=7, fill=mat22, name="Hoop--2")
cell7.region = +surf22 & -surf23 & +surf24 & -surf25

# Cell: Plate-B
cell8 = openmc.Cell(cell_id=8, fill=mat23, name="Plate-B")
cell8.region = -surf3 & +surf5

# Cell: Plate-C
cell9 = openmc.Cell(cell_id=9, fill=mat23, name="Plate-C")
cell9.region = -surf3 & -surf4

# Cell: Plate-A
cell10 = openmc.Cell(cell_id=10, fill=mat23, name="Plate-A")
cell10.region = -surf23 & +surf25

# Cell: Plate-B
cell11 = openmc.Cell(cell_id=11, fill=mat23, name="Plate-B")
cell11.region = -surf23 & -surf24

# Cell: Post
cell12 = openmc.Cell(cell_id=12, fill=mat24, name="Post")
cell12.region = -surf1 & +surf4 & -surf5

# Cell: Post
cell13 = openmc.Cell(cell_id=13, fill=mat24, name="Post")
cell13.region = -surf1 & +surf24 & -surf25

# Cell: Part--A
cell14 = openmc.Cell(cell_id=14, fill=mat30, name="Part--A")
cell14.region = -surf11

# Cell: Part--B
cell15 = openmc.Cell(cell_id=15, fill=mat30, name="Part--B")
cell15.region = -surf12

# Cell: Part--C
cell16 = openmc.Cell(cell_id=16, fill=mat30, name="Part--C")
cell16.region = -surf13

# Cell: Part--D
cell17 = openmc.Cell(cell_id=17, fill=mat30, name="Part--D")
cell17.region = -surf14

# Cell: Part--E
cell18 = openmc.Cell(cell_id=18, fill=mat30, name="Part--E")
cell18.region = -surf15

# Cell: Part--F
cell19 = openmc.Cell(cell_id=19, fill=mat30, name="Part--F")
cell19.region = -surf16

# Cell: Leg1
cell20 = openmc.Cell(cell_id=20, fill=mat30, name="Leg1")
cell20.region = +surf23 & +surf31 & -surf32 & +surf40

# Cell: Leg2
cell21 = openmc.Cell(cell_id=21, fill=mat30, name="Leg2")
cell21.region = +surf23 & +surf33 & -surf34 & +surf40

# Cell: Leg3
cell22 = openmc.Cell(cell_id=22, fill=mat30, name="Leg3")
cell22.region = +surf23 & +surf35 & -surf36 & +surf40

# Cell: Leg4
cell23 = openmc.Cell(cell_id=23, fill=mat30, name="Leg4")
cell23.region = +surf23 & +surf37 & -surf38 & +surf40

# Cell: Plate-G
cell24 = openmc.Cell(cell_id=24, fill=mat30, name="Plate-G")
cell24.region = -surf40

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary
# TODO: Adjust dimensions to encompass your entire geometry
boundary_box = openmc.model.RectangularParallelepiped(
    -200, 200, -200, 200, -200, 200,  # xmin, xmax, ymin, ymax, zmin, zmax
    boundary_type="vacuum")

# Create outer void cell (everything outside geometry but inside boundary)
# Particles are killed at the vacuum boundary
outer_region = -boundary_box
outer_region = outer_region & ~cell0.region
outer_region = outer_region & ~cell1.region
outer_region = outer_region & ~cell2.region
outer_region = outer_region & ~cell3.region
outer_region = outer_region & ~cell4.region
outer_region = outer_region & ~cell5.region
outer_region = outer_region & ~cell6.region
outer_region = outer_region & ~cell7.region
outer_region = outer_region & ~cell8.region
outer_region = outer_region & ~cell9.region
outer_region = outer_region & ~cell10.region
outer_region = outer_region & ~cell11.region
outer_region = outer_region & ~cell12.region
outer_region = outer_region & ~cell13.region
outer_region = outer_region & ~cell14.region
outer_region = outer_region & ~cell15.region
outer_region = outer_region & ~cell16.region
outer_region = outer_region & ~cell17.region
outer_region = outer_region & ~cell18.region
outer_region = outer_region & ~cell19.region
outer_region = outer_region & ~cell20.region
outer_region = outer_region & ~cell21.region
outer_region = outer_region & ~cell22.region
outer_region = outer_region & ~cell23.region
outer_region = outer_region & ~cell24.region
outer_cell = openmc.Cell(cell_id=25, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24, outer_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 15000
settings.batches = 4400
settings.inactive = 100
settings.run_mode = "eigenvalue"

# Source definition
source = openmc.IndependentSource()
source.space = openmc.stats.Point((3.0, 0.0, 3.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)
settings.source = source

# Enable delayed neutron kinetics and alpha eigenvalue calculations
settings.calculate_prompt_k = True
settings.calculate_alpha = True

settings.export_to_xml()

# ==============================================================================
# Tallies
# ==============================================================================

tallies = openmc.Tallies()
tallies.export_to_xml()

# ==============================================================================
# Run OpenMC
# ==============================================================================

openmc.run()
