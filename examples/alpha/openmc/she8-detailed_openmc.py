"""
SHE-8
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

# fuel disk
mat1 = openmc.Material(material_id=1, name="fuel disk")
mat1.set_density("atom/b-cm", 1.0)  # TODO: Verify density

# graphite sleeve
mat2 = openmc.Material(material_id=2, name="graphite sleeve")
mat2.set_density("atom/b-cm", 1.0)  # TODO: Verify density

# absorber (boron)
mat3 = openmc.Material(material_id=3, name="absorber (boron)")
mat3.set_density("atom/b-cm", 1.0)  # TODO: Verify density

# aluminum
mat4 = openmc.Material(material_id=4, name="aluminum")
mat4.set_density("atom/b-cm", 2.700000e+00)
mat4.add_nuclide("F", 2.7)

materials = openmc.Materials([mat1, mat2, mat3, mat4])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# fuel disk
surf1 = openmc.ZCylinder(surface_id=1, r=2.225)

# gap
surf2 = openmc.ZCylinder(surface_id=2, r=2.260)

# fuel sleeve
surf3 = openmc.ZCylinder(surface_id=3, r=2.735)

# gap
surf4 = openmc.ZCylinder(surface_id=4, r=2.760)

# matrix tube
surf5 = openmc.ZCylinder(surface_id=5, r=3.250)

# replaceable reflector
surf6 = openmc.ZCylinder(surface_id=6, r=2.750)

# gap
surf7 = openmc.ZCylinder(surface_id=7, r=2.760)

# matrix tube
surf8 = openmc.ZCylinder(surface_id=8, r=3.250)

# permanent reflector
surf9 = openmc.ZCylinder(surface_id=9, r=3.250)

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

surf99 = openmc.model.RectangularParallelepiped(
    -208.0, 208.0, -121.48, 121.48, -120.0, 120.0, surface_id=99)


# ==============================================================================
# Universes (from COG define unit blocks)
# ==============================================================================

# Unit 1: fuel element
u1_cell0 = openmc.Cell(fill=mat1, name="fuel")
u1_cell0.region = -surf1 & +surf13 & -surf14 & -surf1 & +surf16 & -surf17
u1_cell1 = openmc.Cell(fill=mat2, name="graphite")
u1_cell1.region = -surf3 & +surf10 & -surf11 & -surf3 & +surf19 & -surf20
u1_cell2 = openmc.Cell(fill=mat2, name="graphite")
u1_cell2.region = -surf2 & +surf11 & -surf12 & -surf1 & +surf18 & -surf19
u1_cell3 = openmc.Cell(fill=mat2, name="graphite")
u1_cell3.region = -surf1 & +surf14 & -surf16
u1_cell4 = openmc.Cell(fill=mat2, name="graphite")
u1_cell4.region = -surf1 & +surf12 & -surf13 & -surf1 & +surf17 & -surf18
u1_cell5 = openmc.Cell(fill=mat2, name="graphite")
u1_cell5.region = +surf2 & -surf3 & +surf11 & -surf19
u1_cell6 = openmc.Cell(fill=mat2, name="graphite")
u1_cell6.region = +surf4 & -surf5 & +surf11 & -surf19
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2, u1_cell3, u1_cell4, u1_cell5, u1_cell6])

# Unit 2: replacement reflector element
u2_cell0 = openmc.Cell(fill=mat2, name="graphite")
u2_cell0.region = -surf6 & +surf10 & -surf20
u2_cell1 = openmc.Cell(fill=mat2, name="graphite")
u2_cell1.region = +surf4 & -surf5 & +surf11 & -surf19
universe2 = openmc.Universe(universe_id=2, cells=[u2_cell0, u2_cell1])

# Unit 3: permanent reflector element
u3_cell0 = openmc.Cell(fill=mat2, name="graphite")
u3_cell0.region = -surf8 & +surf10 & -surf20
universe3 = openmc.Universe(universe_id=3, cells=[u3_cell0])

# Unit 4: control rod - black - 50/30
u4_cell0 = openmc.Cell(fill=mat3, name="boron")
u4_cell0.region = +surf21 & -surf22 & +surf10 & -surf20
u4_cell1 = openmc.Cell(fill=mat4, name="aluminum")
u4_cell1.region = +surf23 & -surf24 & +surf10 & -surf20
u4_cell2 = openmc.Cell(fill=mat2, name="graphite")
u4_cell2.region = +surf25 & -surf26 & +surf11 & -surf19
universe4 = openmc.Universe(universe_id=4, cells=[u4_cell0, u4_cell1, u4_cell2])

# Unit 5: control rod - gray - 42-1
u5_cell0 = openmc.Cell(fill=mat3, name="boron")
u5_cell0.region = -surf27 & +surf10 & -surf20
u5_cell1 = openmc.Cell(fill=mat4, name="aluminum")
u5_cell1.region = +surf28 & -surf29 & +surf10 & -surf20
u5_cell2 = openmc.Cell(fill=mat2, name="graphite")
u5_cell2.region = +surf30 & -surf31 & +surf11 & -surf19
universe5 = openmc.Universe(universe_id=5, cells=[u5_cell0, u5_cell1, u5_cell2])

# TODO: Incomplete lattice data for unit 6

# Unit 9: lattice void
u9_cell0 = openmc.Cell(fill=mat0, name="void")
u9_cell0.region = -surf8 & +surf10 & -surf20
universe9 = openmc.Universe(universe_id=9, cells=[u9_cell0])

# Cell using unit 6: array
cell0 = openmc.Cell(cell_id=0, fill=universe6, name="array")
cell0.region = -surf99

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 20000
settings.batches = 1000
settings.inactive = 50
settings.run_mode = "eigenvalue"

# Source definition
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)
settings.source = source

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
