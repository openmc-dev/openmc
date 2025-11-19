"""
HEU-MET-INTER-006-1:  125.6 kg U(93.2) @ C/X=51.2 (1st Zeus Exp't)
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Oralloy
mat1 = openmc.Material(material_id=1, name="Oralloy")
mat1.set_density("atom/b-cm", 4.814600e-02)
mat1.add_nuclide("U234", 4.9483e-4)
mat1.add_nuclide("U235", 4.4918e-2)
mat1.add_nuclide("U236", 1.5917e-4)
mat1.add_nuclide("U238", 2.5740e-3)

# Graphite
mat2 = openmc.Material(material_id=2, name="Graphite")
mat2.set_density("atom/b-cm", 1.702900e+00)
mat2.add_element("C", 1.7029)
mat2.add_s_alpha_beta("c_Graphite")

# Copper
mat3 = openmc.Material(material_id=3, name="Copper")
mat3.set_density("atom/b-cm", 8.735100e+00)
mat3.add_element("Cu", 8.7351)

# Al-6061
mat4 = openmc.Material(material_id=4, name="Al-6061")
mat4.set_density("atom/b-cm", 1.020657e+02)
mat4.add_element("P", 2.6657)
mat4.add_element("Al", 97.175)
mat4.add_element("Mg", 1.00)
mat4.add_element("Si", 0.60)
mat4.add_element("Fe", 0.35)
mat4.add_element("Cu", 0.275)

materials = openmc.Materials([mat1, mat2, mat3, mat4])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Oy/outer/1
surf1 = openmc.ZCylinder(surface_id=1, r=26.67)

# Oy/outer/2
surf2 = openmc.ZCylinder(surface_id=2, r=26.67)

# Oy/outer/3
surf3 = openmc.ZCylinder(surface_id=3, r=26.67)

# Oy/outer/4
surf4 = openmc.ZCylinder(surface_id=4, r=26.67)

# Oy/outer/5
surf5 = openmc.ZCylinder(surface_id=5, r=26.67)

# Oy/outer/6
surf6 = openmc.ZCylinder(surface_id=6, r=26.67)

# Oy/outer/7
surf7 = openmc.ZCylinder(surface_id=7, r=26.67)

# Oy/outer/8
surf8 = openmc.ZCylinder(surface_id=8, r=26.67)

# Oy/outer/9
surf9 = openmc.ZCylinder(surface_id=9, r=26.67)

# Oy/outer/10
surf10 = openmc.ZCylinder(surface_id=10, r=26.67)

# Oy/inner/1-4
surf11 = openmc.ZCylinder(surface_id=11, r=3.175)

# Graphite/outer
surf12 = openmc.ZCylinder(surface_id=12, r=26.67)

# Copper/ram/outer
surf13 = openmc.ZCylinder(surface_id=13, r=26.67)

# Copper/inner
surf14 = openmc.ZCylinder(surface_id=14, r=26.797)

# Copper/void
# Box (6 planes): xmin=-27.94, xmax=27.94, ymin=-27.94, ymax=27.94, zmin=123.0376, zmax=125.0376
surf15_xmin = openmc.XPlane(surface_id=10000, x0=-27.94)
surf15_xmax = openmc.XPlane(surface_id=10001, x0=27.94)
surf15_ymin = openmc.YPlane(surface_id=10002, y0=-27.94)
surf15_ymax = openmc.YPlane(surface_id=10003, y0=27.94)
surf15_zmin = openmc.ZPlane(surface_id=10004, z0=123.0376)
surf15_zmax = openmc.ZPlane(surface_id=10005, z0=125.0376)

# Copper/outer
# Box (6 planes): xmin=-44.1452, xmax=44.1452, ymin=-44.1452, ymax=44.1452, zmin=0.0, zmax=123.9012
surf16_xmin = openmc.XPlane(surface_id=10006, x0=-44.1452)
surf16_xmax = openmc.XPlane(surface_id=10007, x0=44.1452)
surf16_ymin = openmc.YPlane(surface_id=10008, y0=-44.1452)
surf16_ymax = openmc.YPlane(surface_id=10009, y0=44.1452)
surf16_zmin = openmc.ZPlane(surface_id=10010, z0=0.0)
surf16_zmax = openmc.ZPlane(surface_id=10011, z0=123.9012)

# Platen/Inner
surf17 = openmc.ZCylinder(surface_id=17, r=4.7625)

# Platen/Outer
surf18 = openmc.ZCylinder(surface_id=18, r=26.67)

# Alignment Tube/Inner
surf19 = openmc.ZCylinder(surface_id=19, r=2.54)

# Alignment Tube/Outer
surf20 = openmc.ZCylinder(surface_id=20, r=3.1496)


# Cell: HEU
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="HEU")
cell0.region = -surf1 & +surf11 & -surf12

# Cell: HEU
cell1 = openmc.Cell(cell_id=1, fill=mat1, name="HEU")
cell1.region = -surf2 & +surf11 & -surf12

# Cell: HEU
cell2 = openmc.Cell(cell_id=2, fill=mat1, name="HEU")
cell2.region = -surf3 & +surf11 & -surf12

# Cell: HEU
cell3 = openmc.Cell(cell_id=3, fill=mat1, name="HEU")
cell3.region = -surf4 & +surf11 & -surf12

# Cell: HEU
cell4 = openmc.Cell(cell_id=4, fill=mat1, name="HEU")
cell4.region = -surf5 & +surf11 & -surf12

# Cell: HEU
cell5 = openmc.Cell(cell_id=5, fill=mat1, name="HEU")
cell5.region = -surf6 & +surf11 & -surf12

# Cell: HEU
cell6 = openmc.Cell(cell_id=6, fill=mat1, name="HEU")
cell6.region = -surf7 & +surf11 & -surf12

# Cell: HEU
cell7 = openmc.Cell(cell_id=7, fill=mat1, name="HEU")
cell7.region = -surf8 & +surf11 & -surf12

# Cell: HEU
cell8 = openmc.Cell(cell_id=8, fill=mat1, name="HEU")
cell8.region = -surf9 & +surf11 & -surf12

# Cell: HEU
cell9 = openmc.Cell(cell_id=9, fill=mat1, name="HEU")
cell9.region = -surf10 & +surf11 & -surf12

# Cell: C
cell10 = openmc.Cell(cell_id=10, fill=mat2, name="C")
cell10.region = +surf11 & -surf12 & +surf1 & +surf2 & +surf3 & +surf4 & +surf5 & +surf6 & +surf7 & +surf8 & +surf9 & +surf10

# Cell: Cu
cell11 = openmc.Cell(cell_id=11, fill=mat3, name="Cu")
cell11.region = +surf11 & -surf13

# Cell: Cu
cell12 = openmc.Cell(cell_id=12, fill=mat3, name="Cu")
cell12.region = +surf14 & (-surf15_xmin | +surf15_xmax | -surf15_ymin | +surf15_ymax | -surf15_zmin | +surf15_zmax) & (+surf16_xmin & -surf16_xmax & +surf16_ymin & -surf16_ymax & +surf16_zmin & -surf16_zmax)

# Cell: Al
cell13 = openmc.Cell(cell_id=13, fill=mat4, name="Al")
cell13.region = +surf17 & -surf18

# Cell: Al
cell14 = openmc.Cell(cell_id=14, fill=mat4, name="Al")
cell14.region = +surf19 & -surf20

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary (6 planes)
# TODO: Adjust dimensions to encompass your entire geometry
boundary_xmin = openmc.XPlane(surface_id=10012, x0=-200, boundary_type="vacuum")
boundary_xmax = openmc.XPlane(surface_id=10013, x0=200, boundary_type="vacuum")
boundary_ymin = openmc.YPlane(surface_id=10014, y0=-200, boundary_type="vacuum")
boundary_ymax = openmc.YPlane(surface_id=10015, y0=200, boundary_type="vacuum")
boundary_zmin = openmc.ZPlane(surface_id=10016, z0=-200, boundary_type="vacuum")
boundary_zmax = openmc.ZPlane(surface_id=10017, z0=200, boundary_type="vacuum")

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
outer_cell = openmc.Cell(cell_id=15, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, outer_cell])
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
source.space = openmc.stats.Point((0.0, 0.0, 70.2))
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
