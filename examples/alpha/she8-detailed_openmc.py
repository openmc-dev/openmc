"""
SHE-8
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# fuel disk
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U235", 8.446100e-05)
mat1.add_nuclide("U238", 3.367100e-04)
mat1.add_nuclide("B10", 2.949400e-09)
mat1.add_nuclide("B11", 1.282300e-08)
mat1.add_nuclide("C12", 9.477700e-02)
mat1.add_nuclide("O16", 8.947900e-04)
mat1.add_element("Si", 2.024800e-07)
mat1.add_element("Al", 6.322800e-08)
mat1.add_element("Fe", 6.109400e-08)
mat1.add_element("Cr", 2.187300e-08)
mat1.add_element("Mn", 2.070100e-08)
mat1.add_element("Ni", 1.937100e-08)
mat1.add_element("Mo", 2.370900e-08)
mat1.add_element("W", 3.093200e-08)
mat1.add_element("Cd", 2.023700e-10)
mat1.add_element("H", 7.965500e-05)

# graphite sleeve
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("B10", 1.675200e-09)
mat2.add_nuclide("B11", 7.516500e-09)
mat2.add_nuclide("C12", 8.071900e-02)
mat2.add_nuclide("O16", 3.288800e-05)
mat2.add_element("N", 2.043900e-06)
mat2.add_element("H", 6.467700e-05)

# absorber (boron)
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("B10", 8.097100e-05)
mat3.add_nuclide("B11", 3.520300e-04)
mat3.add_nuclide("C12", 8.743700e-02)
mat3.add_nuclide("O16", 3.499700e-05)
mat3.add_element("H", 6.998900e-05)

# aluminum
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_element("Al", 1.000000e+00)

materials = openmc.Materials([mat1, mat2, mat3, mat4])

# ==============================================================================
# Geometry
# ==============================================================================

# fuel disk
surf1 = openmc.ZCylinder(surface_id=1, r=2.225)
# gap
surf2 = openmc.ZCylinder(surface_id=2, r=2.26)
# fuel sleeve
surf3 = openmc.ZCylinder(surface_id=3, r=2.735)
# gap
surf4 = openmc.ZCylinder(surface_id=4, r=2.76)
# matrix tube
surf5 = openmc.ZCylinder(surface_id=5, r=3.25)
# replaceable reflector
surf6 = openmc.ZCylinder(surface_id=6, r=2.75)
# gap
surf7 = openmc.ZCylinder(surface_id=7, r=2.76)
# matrix tube
surf8 = openmc.ZCylinder(surface_id=8, r=3.25)
# permanent reflector
surf9 = openmc.ZCylinder(surface_id=9, r=3.25)
surf10 = openmc.ZPlane(surface_id=10, z0=-120)
surf11 = openmc.ZPlane(surface_id=11, z0=-118)
surf12 = openmc.ZPlane(surface_id=12, z0=-117.5)
surf13 = openmc.ZPlane(surface_id=13, z0=-116.5)
surf14 = openmc.ZPlane(surface_id=14, z0=-0.5)
surf15 = openmc.ZPlane(surface_id=15, z0=0)
surf16 = openmc.ZPlane(surface_id=16, z0=0.5)
surf17 = openmc.ZPlane(surface_id=17, z0=116.5)
surf18 = openmc.ZPlane(surface_id=18, z0=117.5)
surf19 = openmc.ZPlane(surface_id=19, z0=118)
surf20 = openmc.ZPlane(surface_id=20, z0=120)
# void
surf21 = openmc.ZCylinder(surface_id=21, r=1.5)
# absorber
surf22 = openmc.ZCylinder(surface_id=22, r=2.5)
# gap
surf23 = openmc.ZCylinder(surface_id=23, r=2.55)
# aluminum tube (approximated due to drawing typo)
surf24 = openmc.ZCylinder(surface_id=24, r=2.65)
# gap
surf25 = openmc.ZCylinder(surface_id=25, r=2.76)
# matrix tube
surf26 = openmc.ZCylinder(surface_id=26, r=3.25)
# absorber
surf27 = openmc.ZCylinder(surface_id=27, r=2.1)
# gap
surf28 = openmc.ZCylinder(surface_id=28, r=2.24)
# aluminum tube
surf29 = openmc.ZCylinder(surface_id=29, r=2.5)
# gap
surf30 = openmc.ZCylinder(surface_id=30, r=2.76)
# matrix tube
surf31 = openmc.ZCylinder(surface_id=31, r=3.25)
surf99 = openmc.model.RectangularParallelepiped(-208.0, 208.0, -121.48, 121.48, -120.0, 120.0)

# ------------------------------------------------------------------------------
# Universes
# ------------------------------------------------------------------------------

u1_cell0 = openmc.Cell(fill=mat1)
u1_cell0.region = -surf1 & +surf13 & -surf14 & -surf1 & +surf16 & -surf17
u1_cell1 = openmc.Cell(fill=mat2)
u1_cell1.region = -surf3 & +surf10 & -surf11 & -surf3 & +surf19 & -surf20
u1_cell2 = openmc.Cell(fill=mat2)
u1_cell2.region = -surf2 & +surf11 & -surf12 & -surf1 & +surf18 & -surf19
u1_cell3 = openmc.Cell(fill=mat2)
u1_cell3.region = -surf1 & +surf14 & -surf16
u1_cell4 = openmc.Cell(fill=mat2)
u1_cell4.region = -surf1 & +surf12 & -surf13 & -surf1 & +surf17 & -surf18
u1_cell5 = openmc.Cell(fill=mat2)
u1_cell5.region = +surf2 & -surf3 & +surf11 & -surf19
u1_cell6 = openmc.Cell(fill=mat2)
u1_cell6.region = +surf4 & -surf5 & +surf11 & -surf19
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2, u1_cell3, u1_cell4, u1_cell5, u1_cell6])

u2_cell0 = openmc.Cell(fill=mat2)
u2_cell0.region = -surf6 & +surf10 & -surf20
u2_cell1 = openmc.Cell(fill=mat2)
u2_cell1.region = +surf4 & -surf5 & +surf11 & -surf19
universe2 = openmc.Universe(universe_id=2, cells=[u2_cell0, u2_cell1])

u3_cell0 = openmc.Cell(fill=mat2)
u3_cell0.region = -surf8 & +surf10 & -surf20
universe3 = openmc.Universe(universe_id=3, cells=[u3_cell0])

u4_cell0 = openmc.Cell(fill=mat3)
u4_cell0.region = +surf21 & -surf22 & +surf10 & -surf20
u4_cell1 = openmc.Cell(fill=mat4)
u4_cell1.region = +surf23 & -surf24 & +surf10 & -surf20
u4_cell2 = openmc.Cell(fill=mat2)
u4_cell2.region = +surf25 & -surf26 & +surf11 & -surf19
universe4 = openmc.Universe(universe_id=4, cells=[u4_cell0, u4_cell1, u4_cell2])

u5_cell0 = openmc.Cell(fill=mat3)
u5_cell0.region = -surf27 & +surf10 & -surf20
u5_cell1 = openmc.Cell(fill=mat4)
u5_cell1.region = +surf28 & -surf29 & +surf10 & -surf20
u5_cell2 = openmc.Cell(fill=mat2)
u5_cell2.region = +surf30 & -surf31 & +surf11 & -surf19
universe5 = openmc.Universe(universe_id=5, cells=[u5_cell0, u5_cell1, u5_cell2])

universe6 = openmc.Universe(universe_id=6, cells=[])

u9_cell0 = openmc.Cell()
u9_cell0.region = -surf8 & +surf10 & -surf20
universe9 = openmc.Universe(universe_id=9, cells=[u9_cell0])

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# array
cell1 = openmc.Cell(cell_id=1, fill=universe6)
cell1.region = -surf99

root_universe = openmc.Universe(cells=[cell1])
geometry = openmc.Geometry(root_universe)

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 20000
settings.batches = 1000
settings.inactive = 50
settings.run_mode = "eigenvalue"

source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
