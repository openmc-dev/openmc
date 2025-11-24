"""
HMF073-1: ZEUS - Unmoderated
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# HEU outer rigns
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 4.870700e-04)
mat1.add_nuclide("U235", 4.457400e-02)
mat1.add_nuclide("U236", 2.067500e-04)
mat1.add_nuclide("U238", 2.542400e-03)

# HEU inner discs
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("U234", 5.037700e-04)
mat2.add_nuclide("U235", 4.538400e-02)
mat2.add_nuclide("U236", 1.133700e-04)
mat2.add_nuclide("U238", 2.621100e-03)

# Al-6061
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_element("Si", 3.429500e-04)
mat3.add_element("Fe", 1.006100e-04)
mat3.add_element("Cu", 6.947100e-05)
mat3.add_element("Mn", 2.191500e-05)
mat3.add_element("Mg", 6.604900e-04)
mat3.add_element("Cr", 7.718500e-05)
mat3.add_element("Zn", 3.068700e-05)
mat3.add_element("Ti", 2.514600e-05)
mat3.add_element("Al", 5.781600e-02)

# SS304
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_element("C", 2.063700e-04)
mat4.add_element("N", 1.702900e-04)
mat4.add_element("Si", 1.015800e-03)
mat4.add_element("P", 4.227800e-05)
mat4.add_element("S", 5.833200e-06)
mat4.add_element("Cr", 1.644200e-02)
mat4.add_element("Mn", 1.455700e-03)
mat4.add_element("Fe", 5.955400e-02)
mat4.add_element("Ni", 6.454600e-03)
mat4.add_element("Cu", 1.545500e-05)
mat4.add_element("Mo", 1.364900e-05)

# Copper top    reflector
mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")
mat5.add_element("Cu", 8.339400e-02)

# Copper bottom reflector
mat6 = openmc.Material(material_id=6)
mat6.set_density("sum")
mat6.add_element("Cu", 8.331500e-02)

# Copper corner reflector
mat7 = openmc.Material(material_id=7)
mat7.set_density("sum")
mat7.add_element("Cu", 8.295300e-02)

# Copper side   reflector
mat8 = openmc.Material(material_id=8)
mat8.set_density("sum")
mat8.add_element("Cu", 8.278400e-02)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8])

# ==============================================================================
# Geometry
# ==============================================================================

# Cu top    reflector
surf1 = openmc.model.RectangularParallelepiped(-27.94, 27.94, -27.94, 27.94, 59.24296, 73.67016)
# inner contour
surf2 = openmc.ZCylinder(surface_id=2, r=3.175)
# Cu bottom reflector, outer
surf3 = openmc.ZCylinder(surface_id=3, r=26.67)
# Cu corner reflector, inner
surf4 = openmc.ZCylinder(surface_id=4, r=26.797)
# Cu corner reflector, outer
surf5 = openmc.model.RectangularParallelepiped(-27.94, 27.94, -27.94, 27.94, 0.0, 59.24296)
# Cu side   reflector, inner
surf6 = openmc.model.RectangularParallelepiped(-27.94, 27.94, -27.94, 27.94, 0.0, 103.251)
# Cu side   reflector, outer
surf7 = openmc.model.RectangularParallelepiped(-44.1452, 44.1452, -44.1452, 44.1452, 0.0, 103.251)
# SS304 diaphragm
surf8 = openmc.model.RectangularParallelepiped(-27.94, 27.94, -27.94, 27.94, 57.7088, 57.97296)
# inner contour
surf10 = openmc.ZCylinder(surface_id=10, r=7.62635)
# HEU inner disc
surf11 = openmc.ZCylinder(surface_id=11, r=19.05)
# HEU inner disc
surf12 = openmc.ZCylinder(surface_id=12, r=19.05)
# HEU inner disc
surf13 = openmc.ZCylinder(surface_id=13, r=19.05)
# HEU outer disc
surf21 = openmc.ZCylinder(surface_id=21, r=26.67)
# HEU outer disc
surf22 = openmc.ZCylinder(surface_id=22, r=26.67)
# HEU outer disc
surf23 = openmc.ZCylinder(surface_id=23, r=26.67)
# Al-6061 alignment tube, inner
surf31 = openmc.ZCylinder(surface_id=31, r=2.54)
# Al-6061 alignment tube, inner
surf32 = openmc.ZCylinder(surface_id=32, r=3.1496)
# Al-6061 platen, inner
surf33 = openmc.ZCylinder(surface_id=33, r=4.7625)
# Al-6061 platen, outer
surf34 = openmc.ZCylinder(surface_id=34, r=26.67)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Cu
cell1 = openmc.Cell(cell_id=1, fill=mat5)
cell1.region = -surf1

# Cu
cell2 = openmc.Cell(cell_id=2, fill=mat6)
cell2.region = +surf2 & -surf3

# Cu
cell3 = openmc.Cell(cell_id=3, fill=mat7)
cell3.region = +surf1 & +surf4 & -surf5 & +surf8

# Cu
cell4 = openmc.Cell(cell_id=4, fill=mat8)
cell4.region = +surf1 & +surf5 & +surf6 & -surf7 & +surf8

# SS304
cell5 = openmc.Cell(cell_id=5, fill=mat4)
cell5.region = -surf8

# HEU
cell6 = openmc.Cell(cell_id=6, fill=mat1)
cell6.region = +surf8 & +surf11 & -surf21

# HEU
cell7 = openmc.Cell(cell_id=7, fill=mat1)
cell7.region = +surf8 & +surf10 & +surf12 & -surf22 & +surf23

# HEU
cell8 = openmc.Cell(cell_id=8, fill=mat1)
cell8.region = +surf2 & +surf3 & +surf12 & +surf13 & -surf23

# HEU
cell9 = openmc.Cell(cell_id=9, fill=mat2)
cell9.region = +surf8 & -surf11

# HEU
cell10 = openmc.Cell(cell_id=10, fill=mat2)
cell10.region = +surf8 & +surf10 & -surf12

# HEU
cell11 = openmc.Cell(cell_id=11, fill=mat2)
cell11.region = +surf2 & +surf3 & +surf12 & -surf13

# Al6061
cell12 = openmc.Cell(cell_id=12, fill=mat3)
cell12.region = +surf8 & +surf31 & -surf32

# Al6061
cell13 = openmc.Cell(cell_id=13, fill=mat3)
cell13.region = +surf3 & +surf33 & -surf34

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13])
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
source.space = openmc.stats.Box((-9.0, -9.0, 55.9), (9.0, 9.0, 59.5))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
