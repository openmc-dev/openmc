"""
HEU-MET-INTER-006-1:  125.6 kg U(93.2) @ C/X=51.2 (1st Zeus Exp't)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Oralloy
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 4.948300e-04)
mat1.add_nuclide("U235", 4.491800e-02)
mat1.add_nuclide("U236", 1.591700e-04)
mat1.add_nuclide("U238", 2.574000e-03)

# Graphite
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("C", 1.702900e+00)

# Copper
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_element("Cu", 8.735100e+00)

# Al-6061
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_element("P", 2.665700e+00)
mat4.add_element("Al", 9.717500e+01)
mat4.add_element("Mg", 1.000000e+00)
mat4.add_element("Si", 6.000000e-01)
mat4.add_element("Fe", 3.500000e-01)
mat4.add_element("Cu", 2.750000e-01)
mat4.add_element("Cr", 2.500000e-01)
mat4.add_element("Zn", 1.250000e-01)
mat4.add_element("Mg", 7.500000e-02)
mat4.add_element("Ti", 7.500000e-02)
mat4.add_nuclide("O16", 7.500000e-02)

materials = openmc.Materials([mat1, mat2, mat3, mat4])

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
surf15 = openmc.model.RectangularParallelepiped(-27.94, 27.94, -27.94, 27.94, 123.0376, 125.0376)
# Copper/outer
surf16 = openmc.model.RectangularParallelepiped(-44.1452, 44.1452, -44.1452, 44.1452, 0.0, 123.9012)
# Platen/Inner
surf17 = openmc.ZCylinder(surface_id=17, r=4.7625)
# Platen/Outer
surf18 = openmc.ZCylinder(surface_id=18, r=26.67)
# Alignment Tube/Inner
surf19 = openmc.ZCylinder(surface_id=19, r=2.54)
# Alignment Tube/Outer
surf20 = openmc.ZCylinder(surface_id=20, r=3.1496)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# HEU
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1 & +surf11 & -surf12

# HEU
cell2 = openmc.Cell(cell_id=2, fill=mat1)
cell2.region = -surf2 & +surf11 & -surf12

# HEU
cell3 = openmc.Cell(cell_id=3, fill=mat1)
cell3.region = -surf3 & +surf11 & -surf12

# HEU
cell4 = openmc.Cell(cell_id=4, fill=mat1)
cell4.region = -surf4 & +surf11 & -surf12

# HEU
cell5 = openmc.Cell(cell_id=5, fill=mat1)
cell5.region = -surf5 & +surf11 & -surf12

# HEU
cell6 = openmc.Cell(cell_id=6, fill=mat1)
cell6.region = -surf6 & +surf11 & -surf12

# HEU
cell7 = openmc.Cell(cell_id=7, fill=mat1)
cell7.region = -surf7 & +surf11 & -surf12

# HEU
cell8 = openmc.Cell(cell_id=8, fill=mat1)
cell8.region = -surf8 & +surf11 & -surf12

# HEU
cell9 = openmc.Cell(cell_id=9, fill=mat1)
cell9.region = -surf9 & +surf11 & -surf12

# HEU
cell10 = openmc.Cell(cell_id=10, fill=mat1)
cell10.region = -surf10 & +surf11 & -surf12

# C
cell11 = openmc.Cell(cell_id=11, fill=mat2)
cell11.region = +surf11 & -surf12 & +surf1 & +surf2 & +surf3 & +surf4 & +surf5 & +surf6 & +surf7 & +surf8 & +surf9 & +surf10

# Cu
cell12 = openmc.Cell(cell_id=12, fill=mat3)
cell12.region = +surf11 & -surf13

# Cu
cell13 = openmc.Cell(cell_id=13, fill=mat3)
cell13.region = +surf14 & +surf15 & -surf16

# Al
cell14 = openmc.Cell(cell_id=14, fill=mat4)
cell14.region = +surf17 & -surf18

# Al
cell15 = openmc.Cell(cell_id=15, fill=mat4)
cell15.region = +surf19 & -surf20

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15])
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
source.space = openmc.stats.Point((0.0, 0.0, 70.2))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
