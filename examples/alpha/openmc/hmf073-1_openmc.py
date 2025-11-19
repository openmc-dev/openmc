"""
HMF073-1: ZEUS - Unmoderated
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# HEU outer rigns
mat1 = openmc.Material(material_id=1, name="HEU outer rigns")
mat1.set_density("atom/b-cm", 4.781022e-02)
mat1.add_nuclide("U234", 4.8707e-4)
mat1.add_nuclide("U235", 4.4574e-2)
mat1.add_nuclide("U236", 2.0675e-4)
mat1.add_nuclide("U238", 2.5424e-3)

# HEU inner discs
mat2 = openmc.Material(material_id=2, name="HEU inner discs")
mat2.set_density("atom/b-cm", 4.862224e-02)
mat2.add_nuclide("U234", 5.0377e-4)
mat2.add_nuclide("U235", 4.5384e-2)
mat2.add_nuclide("U236", 1.1337e-4)
mat2.add_nuclide("U238", 2.6211e-3)

# Al-6061
mat3 = openmc.Material(material_id=3, name="Al-6061")
mat3.set_density("atom/b-cm", 1.272621e-03)
mat3.add_element("Si", 3.4295e-4)
mat3.add_element("Fe", 1.0061e-4)
mat3.add_element("Cu", 6.9471e-5)
mat3.add_element("Mn", 2.1915e-5)
mat3.add_element("Mg", 6.6049e-4)
mat3.add_element("Cr", 7.7185e-5)

# SS304
mat4 = openmc.Material(material_id=4, name="SS304")
mat4.set_density("atom/b-cm", 1.788257e-02)
mat4.add_element("C", 2.0637e-4)
mat4.add_element("N", 1.7029e-4)
mat4.add_element("Si", 1.0158e-3)
mat4.add_element("P", 4.2278e-5)
mat4.add_element("S", 5.8332e-6)
mat4.add_element("Cr", 1.6442e-2)
mat4.add_s_alpha_beta("c_Graphite")

# Copper top    reflector
mat5 = openmc.Material(material_id=5, name="Copper top    reflector")
mat5.set_density("atom/b-cm", 8.339400e-02)
mat5.add_element("Cu", 8.3394e-2)

# Copper bottom reflector
mat6 = openmc.Material(material_id=6, name="Copper bottom reflector")
mat6.set_density("atom/b-cm", 8.331500e-02)
mat6.add_element("Cu", 8.3315e-2)

# Copper corner reflector
mat7 = openmc.Material(material_id=7, name="Copper corner reflector")
mat7.set_density("atom/b-cm", 8.295300e-02)
mat7.add_element("Cu", 8.2953e-2)

# Copper side   reflector
mat8 = openmc.Material(material_id=8, name="Copper side   reflector")
mat8.set_density("atom/b-cm", 8.278400e-02)
mat8.add_element("Cu", 8.2784e-2)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Cu top    reflector
surf1 = openmc.model.RectangularParallelepiped(
    -27.94, 27.94, -27.94, 27.94, 59.24296, 73.67016, surface_id=1)

# inner contour
surf2 = openmc.ZCylinder(surface_id=2, r=3.175)

# Cu bottom reflector, outer
surf3 = openmc.ZCylinder(surface_id=3, r=26.67)

# Cu corner reflector, inner
surf4 = openmc.ZCylinder(surface_id=4, r=26.797)

# Cu corner reflector, outer
surf5 = openmc.model.RectangularParallelepiped(
    -27.94, 27.94, -27.94, 27.94, 0.0, 59.24296, surface_id=5)

# Cu side   reflector, inner
surf6 = openmc.model.RectangularParallelepiped(
    -27.94, 27.94, -27.94, 27.94, 0.0, 103.251, surface_id=6)

# Cu side   reflector, outer
surf7 = openmc.model.RectangularParallelepiped(
    -44.1452, 44.1452, -44.1452, 44.1452, 0.0, 103.251, surface_id=7)

# SS304 diaphragm
surf8 = openmc.model.RectangularParallelepiped(
    -27.94, 27.94, -27.94, 27.94, 57.7088, 57.97296, surface_id=8)

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


# Cell: HEU
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="HEU")
cell0.region = +surf8 & +surf11 & -surf21

# Cell: HEU
cell1 = openmc.Cell(cell_id=1, fill=mat1, name="HEU")
cell1.region = +surf8 & +surf10 & +surf12 & -surf22 & +surf23

# Cell: HEU
cell2 = openmc.Cell(cell_id=2, fill=mat1, name="HEU")
cell2.region = +surf2 & +surf3 & +surf12 & +surf13 & -surf23

# Cell: HEU
cell3 = openmc.Cell(cell_id=3, fill=mat2, name="HEU")
cell3.region = +surf8 & -surf11

# Cell: HEU
cell4 = openmc.Cell(cell_id=4, fill=mat2, name="HEU")
cell4.region = +surf8 & +surf10 & -surf12

# Cell: HEU
cell5 = openmc.Cell(cell_id=5, fill=mat2, name="HEU")
cell5.region = +surf2 & +surf3 & +surf12 & -surf13

# Cell: Al6061
cell6 = openmc.Cell(cell_id=6, fill=mat3, name="Al6061")
cell6.region = +surf8 & +surf31 & -surf32

# Cell: Al6061
cell7 = openmc.Cell(cell_id=7, fill=mat3, name="Al6061")
cell7.region = +surf3 & +surf33 & -surf34

# Cell: SS304
cell8 = openmc.Cell(cell_id=8, fill=mat4, name="SS304")
cell8.region = -surf8

# Cell: Cu
cell9 = openmc.Cell(cell_id=9, fill=mat5, name="Cu")
cell9.region = -surf1

# Cell: Cu
cell10 = openmc.Cell(cell_id=10, fill=mat6, name="Cu")
cell10.region = +surf2 & -surf3

# Cell: Cu
cell11 = openmc.Cell(cell_id=11, fill=mat7, name="Cu")
cell11.region = +surf1 & +surf4 & -surf5 & +surf8

# Cell: Cu
cell12 = openmc.Cell(cell_id=12, fill=mat8, name="Cu")
cell12.region = +surf1 & +surf5 & +surf6 & -surf7 & +surf8

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
outer_cell = openmc.Cell(cell_id=13, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, outer_cell])
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
source.space = openmc.stats.Point((0.0, 0.0, 58.5))
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
