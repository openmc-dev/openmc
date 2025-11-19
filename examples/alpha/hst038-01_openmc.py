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
# Box (6 planes): xmin=-39.217600000000004, xmax=-38.5826, ymin=-38.5826, ymax=38.5826, zmin=0.6700000000000004, zmax=8.290000000000001
surf11_xmin = openmc.XPlane(surface_id=10000, x0=-39.217600000000004)
surf11_xmax = openmc.XPlane(surface_id=10001, x0=-38.5826)
surf11_ymin = openmc.YPlane(surface_id=10002, y0=-38.5826)
surf11_ymax = openmc.YPlane(surface_id=10003, y0=38.5826)
surf11_zmin = openmc.ZPlane(surface_id=10004, z0=0.6700000000000004)
surf11_zmax = openmc.ZPlane(surface_id=10005, z0=8.290000000000001)

# Part B
# Box (6 planes): xmin=38.5826, xmax=39.217600000000004, ymin=-38.5826, ymax=38.5826, zmin=0.6700000000000004, zmax=8.290000000000001
surf12_xmin = openmc.XPlane(surface_id=10006, x0=38.5826)
surf12_xmax = openmc.XPlane(surface_id=10007, x0=39.217600000000004)
surf12_ymin = openmc.YPlane(surface_id=10008, y0=-38.5826)
surf12_ymax = openmc.YPlane(surface_id=10009, y0=38.5826)
surf12_zmin = openmc.ZPlane(surface_id=10010, z0=0.6700000000000004)
surf12_zmax = openmc.ZPlane(surface_id=10011, z0=8.290000000000001)

# Part C
# Box (6 planes): xmin=-58.42, xmax=58.42, ymin=-39.217600000000004, ymax=-38.5826, zmin=0.6700000000000004, zmax=8.290000000000001
surf13_xmin = openmc.XPlane(surface_id=10012, x0=-58.42)
surf13_xmax = openmc.XPlane(surface_id=10013, x0=58.42)
surf13_ymin = openmc.YPlane(surface_id=10014, y0=-39.217600000000004)
surf13_ymax = openmc.YPlane(surface_id=10015, y0=-38.5826)
surf13_zmin = openmc.ZPlane(surface_id=10016, z0=0.6700000000000004)
surf13_zmax = openmc.ZPlane(surface_id=10017, z0=8.290000000000001)

# Part D
# Box (6 planes): xmin=-58.42, xmax=58.42, ymin=38.5826, ymax=39.217600000000004, zmin=0.6700000000000004, zmax=8.290000000000001
surf14_xmin = openmc.XPlane(surface_id=10018, x0=-58.42)
surf14_xmax = openmc.XPlane(surface_id=10019, x0=58.42)
surf14_ymin = openmc.YPlane(surface_id=10020, y0=38.5826)
surf14_ymax = openmc.YPlane(surface_id=10021, y0=39.217600000000004)
surf14_zmin = openmc.ZPlane(surface_id=10022, z0=0.6700000000000004)
surf14_zmax = openmc.ZPlane(surface_id=10023, z0=8.290000000000001)

# Part E
# Box (6 planes): xmin=48.260000000000005, xmax=58.42, ymin=-58.42, ymax=58.42, zmin=-1.87, zmax=0.67
surf15_xmin = openmc.XPlane(surface_id=10024, x0=48.260000000000005)
surf15_xmax = openmc.XPlane(surface_id=10025, x0=58.42)
surf15_ymin = openmc.YPlane(surface_id=10026, y0=-58.42)
surf15_ymax = openmc.YPlane(surface_id=10027, y0=58.42)
surf15_zmin = openmc.ZPlane(surface_id=10028, z0=-1.87)
surf15_zmax = openmc.ZPlane(surface_id=10029, z0=0.67)

# Part F
# Box (6 planes): xmin=-58.42, xmax=-48.260000000000005, ymin=-58.42, ymax=58.42, zmin=-1.87, zmax=0.67
surf16_xmin = openmc.XPlane(surface_id=10030, x0=-58.42)
surf16_xmax = openmc.XPlane(surface_id=10031, x0=-48.260000000000005)
surf16_ymin = openmc.YPlane(surface_id=10032, y0=-58.42)
surf16_ymax = openmc.YPlane(surface_id=10033, y0=58.42)
surf16_zmin = openmc.ZPlane(surface_id=10034, z0=-1.87)
surf16_zmax = openmc.ZPlane(surface_id=10035, z0=0.67)

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
# Box (6 planes): xmin=-39.37, xmax=39.37, ymin=-39.37, ymax=39.37, zmin=-57.04967, zmax=-55.77967
surf40_xmin = openmc.XPlane(surface_id=10036, x0=-39.37)
surf40_xmax = openmc.XPlane(surface_id=10037, x0=39.37)
surf40_ymin = openmc.YPlane(surface_id=10038, y0=-39.37)
surf40_ymax = openmc.YPlane(surface_id=10039, y0=39.37)
surf40_zmin = openmc.ZPlane(surface_id=10040, z0=-57.04967)
surf40_zmax = openmc.ZPlane(surface_id=10041, z0=-55.77967)


# Cell: Air1
cell0 = openmc.Cell(cell_id=0, fill=None, name="Air1")  # mat0 undefined, using void
cell0.region = +surf23 & -surf31 & -surf32 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

# Cell: Air2
cell1 = openmc.Cell(cell_id=1, fill=None, name="Air2")  # mat0 undefined, using void
cell1.region = +surf23 & -surf33 & -surf34 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

# Cell: Air3
cell2 = openmc.Cell(cell_id=2, fill=None, name="Air3")  # mat0 undefined, using void
cell2.region = +surf23 & -surf35 & -surf36 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

# Cell: Air4
cell3 = openmc.Cell(cell_id=3, fill=None, name="Air4")  # mat0 undefined, using void
cell3.region = +surf23 & -surf37 & -surf38 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

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
cell14.region = (+surf11_xmin & -surf11_xmax & +surf11_ymin & -surf11_ymax & +surf11_zmin & -surf11_zmax)

# Cell: Part--B
cell15 = openmc.Cell(cell_id=15, fill=mat30, name="Part--B")
cell15.region = (+surf12_xmin & -surf12_xmax & +surf12_ymin & -surf12_ymax & +surf12_zmin & -surf12_zmax)

# Cell: Part--C
cell16 = openmc.Cell(cell_id=16, fill=mat30, name="Part--C")
cell16.region = (+surf13_xmin & -surf13_xmax & +surf13_ymin & -surf13_ymax & +surf13_zmin & -surf13_zmax)

# Cell: Part--D
cell17 = openmc.Cell(cell_id=17, fill=mat30, name="Part--D")
cell17.region = (+surf14_xmin & -surf14_xmax & +surf14_ymin & -surf14_ymax & +surf14_zmin & -surf14_zmax)

# Cell: Part--E
cell18 = openmc.Cell(cell_id=18, fill=mat30, name="Part--E")
cell18.region = (+surf15_xmin & -surf15_xmax & +surf15_ymin & -surf15_ymax & +surf15_zmin & -surf15_zmax)

# Cell: Part--F
cell19 = openmc.Cell(cell_id=19, fill=mat30, name="Part--F")
cell19.region = (+surf16_xmin & -surf16_xmax & +surf16_ymin & -surf16_ymax & +surf16_zmin & -surf16_zmax)

# Cell: Leg1
cell20 = openmc.Cell(cell_id=20, fill=mat30, name="Leg1")
cell20.region = +surf23 & +surf31 & -surf32 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

# Cell: Leg2
cell21 = openmc.Cell(cell_id=21, fill=mat30, name="Leg2")
cell21.region = +surf23 & +surf33 & -surf34 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

# Cell: Leg3
cell22 = openmc.Cell(cell_id=22, fill=mat30, name="Leg3")
cell22.region = +surf23 & +surf35 & -surf36 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

# Cell: Leg4
cell23 = openmc.Cell(cell_id=23, fill=mat30, name="Leg4")
cell23.region = +surf23 & +surf37 & -surf38 & (-surf40_xmin | +surf40_xmax | -surf40_ymin | +surf40_ymax | -surf40_zmin | +surf40_zmax)

# Cell: Plate-G
cell24 = openmc.Cell(cell_id=24, fill=mat30, name="Plate-G")
cell24.region = (+surf40_xmin & -surf40_xmax & +surf40_ymin & -surf40_ymax & +surf40_zmin & -surf40_zmax)

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary (6 planes)
# TODO: Adjust dimensions to encompass your entire geometry
boundary_xmin = openmc.XPlane(surface_id=10042, x0=-200, boundary_type="vacuum")
boundary_xmax = openmc.XPlane(surface_id=10043, x0=200, boundary_type="vacuum")
boundary_ymin = openmc.YPlane(surface_id=10044, y0=-200, boundary_type="vacuum")
boundary_ymax = openmc.YPlane(surface_id=10045, y0=200, boundary_type="vacuum")
boundary_zmin = openmc.ZPlane(surface_id=10046, z0=-200, boundary_type="vacuum")
boundary_zmax = openmc.ZPlane(surface_id=10047, z0=200, boundary_type="vacuum")

# Create outer void cell (everything outside geometry but inside boundary)
# Particles are killed at the vacuum boundary
outer_region = +boundary_xmin & -boundary_xmax & +boundary_ymin & -boundary_ymax & +boundary_zmin & -boundary_zmax
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
