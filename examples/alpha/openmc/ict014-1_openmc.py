"""
ICT014-1: RA-6
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

# Al-U3Si2
mat1 = openmc.Material(material_id=1, name="Al-U3Si2")
mat1.set_density("atom/b-cm", 5.091724e-02)
mat1.add_nuclide("U234", 1.7680e-5)
mat1.add_nuclide("U235", 2.4170e-3)
mat1.add_nuclide("U236", 1.3635e-5)
mat1.add_nuclide("U238", 9.6538e-3)
mat1.add_nuclide("Si", 8.3622e-3)
mat1.add_nuclide("Al", 3.0451e-2)
mat1.add_nuclide("B", 1.9296e-6)

# Al6061
mat2 = openmc.Material(material_id=2, name="Al6061")
mat2.set_density("atom/b-cm", 6.004096e-02)
mat2.add_nuclide("Al", 5.8811e-2)
mat2.add_nuclide("Cu", 6.6527e-5)
mat2.add_nuclide("Cr", 3.4398e-5)
mat2.add_nuclide("Mg", 6.6229e-4)
mat2.add_nuclide("Si", 3.8789e-4)
mat2.add_nuclide("Zn", 2.4866e-7)
mat2.add_nuclide("Fe", 7.8610e-5)

# H2O
mat3 = openmc.Material(material_id=3, name="H2O")
mat3.set_density("atom/b-cm", 3.333800e-02)
mat3.add_nuclide("O16", 3.3338e-2)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Cd
mat4 = openmc.Material(material_id=4, name="Cd")
mat4.set_density("atom/b-cm", 4.634000e-02)
mat4.add_nuclide("Cd", 4.6340e-2)

# Al
mat5 = openmc.Material(material_id=5, name="Al")
mat5.set_density("atom/b-cm", 6.026200e-02)
mat5.add_nuclide("Al", 6.0262e-2)

# Ag-In-Cd
mat6 = openmc.Material(material_id=6, name="Ag-In-Cd")
mat6.set_density("atom/b-cm", 5.617460e-02)
mat6.add_nuclide("Ag", 4.5365e-2)
mat6.add_nuclide("In", 7.9765e-3)
mat6.add_nuclide("Cd", 2.8331e-3)

# Stainless steel
mat7 = openmc.Material(material_id=7, name="Stainless steel")
mat7.set_density("atom/b-cm", 8.728029e-02)
mat7.add_nuclide("Fe", 5.9899e-2)
mat7.add_nuclide("Ni", 8.1984e-3)
mat7.add_nuclide("Cr", 1.7582e-2)
mat7.add_nuclide("Si", 6.4246e-4)
mat7.add_nuclide("Mn", 8.7583e-4)
mat7.add_nuclide("C0", 6.0091e-5)
mat7.add_nuclide("S", 2.2505e-5)

# Al2O3
mat8 = openmc.Material(material_id=8, name="Al2O3")
mat8.set_density("atom/b-cm", 1.169440e-01)
mat8.add_nuclide("Al", 4.6778e-2)
mat8.add_nuclide("O16", 7.0166e-2)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Al6061 side plates: inner
surf11 = openmc.model.RectangularParallelepiped(
    -3.3, 3.3, -4, 4, -38.35, 38.35, surface_id=11)

# Al6061 side plates: outer
surf12 = openmc.model.RectangularParallelepiped(
    -3.8, 3.8, -4, 4, -38.35, 38.35, surface_id=12)

# Al-U3Si2 fuel
surf21 = openmc.model.RectangularParallelepiped(
    -3, 3, -0.02533, 0.02533, -31.15, 31.15, surface_id=21)

# Al6061 clad: inner fuel plates
surf22 = openmc.model.RectangularParallelepiped(
    -3.525, 3.525, -0.0745, 0.0745, -33.55, 33.55, surface_id=22)

# Water slot:  inner fuel plates
surf23 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=23)

# Al6061 clad: outer fuel plates
surf24 = openmc.model.RectangularParallelepiped(
    -3.525, 3.525, -0.0745, 0.0745, -38.35, 35.15, surface_id=24)

# Water slot:  outer fuel plates
surf25 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -38.35, 35.15, surface_id=25)

# Cd wire
surf26 = openmc.ZCylinder(surface_id=26, r=0.02425)

# Al6061 crossbar
surf29 = openmc.XCylinder(surface_id=29, r=0.625)

# Al6061 nozzle: inner
surf30 = openmc.ZCylinder(surface_id=30, r=2.4895)

# Al6061 nozzle: outer
surf31 = openmc.ZCylinder(surface_id=31, r=3.0895)

# Al6061 side plates: inner
surf41 = openmc.model.RectangularParallelepiped(
    -3.3, 3.3, -4, 4, -38.35, 39.65, surface_id=41)

# Al6061 side plates: outer
surf42 = openmc.model.RectangularParallelepiped(
    -3.8, 3.8, -4, 4, -38.35, 39.65, surface_id=42)

# Al6061 internal guide plate
surf43 = openmc.model.RectangularParallelepiped(
    -3.525, 3.525, -0.065, 0.065, -33.55, 33.55, surface_id=43)

# Al6061 external guide plate
surf44 = openmc.model.RectangularParallelepiped(
    -3.525, 3.525, -0.065, 0.065, -38.35, 39.65, surface_id=44)

# In-Ag-Cd
surf50 = openmc.model.RectangularParallelepiped(
    -3.09, 3.09, -0.11, 0.11, 31.55, 94.95, surface_id=50)

# edge
surf51 = openmc.XPlane(surface_id=51, x0=-2.98)

# edge
surf52 = openmc.XPlane(surface_id=52, x0=2.98)

# edge
surf53 = openmc.ZCylinder(surface_id=53, r=0.11)

# edge
surf54 = openmc.ZCylinder(surface_id=54, r=0.11)

# indentation
# COG surface type "s" with parameters: 0.1625 tr 0 0.1725 63.25
# This surface type requires manual translation to OpenMC
surf55 = openmc.Sphere(surface_id=55, r=1.0)  # PLACEHOLDER - REPLACE THIS

# SS304L cladding: inner
surf56 = openmc.model.RectangularParallelepiped(
    -3.175, 3.175, -0.145, 0.145, 31.55, 94.95, surface_id=56)

# SS304L cladding: outer
surf57 = openmc.model.RectangularParallelepiped(
    -3.245, 3.245, -0.215, 0.215, 30.45, 117.15, surface_id=57)

# BNCT
surf60 = openmc.model.RectangularParallelepiped(
    -38.55, 38.55, -199.0, -32.40, -40.75, 41.60, surface_id=60)

# Al/Cd
surf61 = openmc.YPlane(surface_id=61, y0=-49.50)

# Cd/Al
surf62 = openmc.YPlane(surface_id=62, y0=-49.55)

# Al/Cd
surf63 = openmc.YPlane(surface_id=63, y0=-59.55)

# Cd/Al2O3
surf64 = openmc.YPlane(surface_id=64, y0=-59.70)

# Pure aluminum grid plate
surf90 = openmc.model.RectangularParallelepiped(
    -30.8, 30.8, -40.5, 40.5, -60.75, -40.75, surface_id=90)

# Primary water hole
surf91 = openmc.ZCylinder(surface_id=91, r=3.0895)

# 1st (outer) plate and slot
surf101 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -38.35, 35.15, surface_id=101)

# 2nd (inner) plate and slot
surf102 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=102)

# 3rd (inner) plate and slot
surf103 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=103)

# 4th (inner) plate and slot
surf104 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=104)

# 5th (inner) plate and slot
surf105 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=105)

# 6th (inner) plate and slot
surf106 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=106)

# 7th (inner) plate and slot
surf107 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=107)

# 8th (inner) plate and slot
surf108 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=108)

# 9th (inner) plate and slot
surf109 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=109)

# 11th (inner) plate and slot
surf111 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=111)

# 12th (inner) plate and slot
surf112 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=112)

# 13th (inner) plate and slot
surf113 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=113)

# 14th (inner) plate and slot
surf114 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=114)

# 15th (inner) plate and slot
surf115 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=115)

# 16th (inner) plate and slot
surf116 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=116)

# 17th (inner) plate and slot
surf117 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=117)

# 18th (inner) plate and slot
surf118 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 35.15, surface_id=118)

# 19th (outer) plate and slot
surf119 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -38.35, 35.15, surface_id=119)

# Slot for Cd wire: outer fuel plate
surf121 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -38.35, 31.65, surface_id=121)

# Slot for Cd wire: inner fuel plate
surf123 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=123)

# Slot for Cd wire: inner fuel plate
surf125 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=125)

# Slot for Cd wire: inner fuel plate
surf127 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=127)

# Slot for Cd wire: inner fuel plate
surf129 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=129)

# Slot for Cd wire: inner fuel plate
surf131 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=131)

# Slot for Cd wire: inner fuel plate
surf133 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=133)

# Slot for Cd wire: inner fuel plate
surf135 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=135)

# Slot for Cd wire: inner fuel plate
surf137 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 31.65, surface_id=137)

# Slot for Cd wire: outer fuel plate
surf139 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -38.35, 31.65, surface_id=139)

surf201 = openmc.ZCylinder(surface_id=201, r=1.1125)

surf206 = openmc.ZCylinder(surface_id=206, r=1.1125)

surf211 = openmc.ZCylinder(surface_id=211, r=1.1125)

surf216 = openmc.ZCylinder(surface_id=216, r=1.1125)

surf221 = openmc.ZCylinder(surface_id=221, r=1.1125)

surf226 = openmc.ZCylinder(surface_id=226, r=1.1125)

surf231 = openmc.ZCylinder(surface_id=231, r=1.1125)

surf236 = openmc.ZCylinder(surface_id=236, r=1.1125)

surf241 = openmc.ZCylinder(surface_id=241, r=1.1125)

surf246 = openmc.ZCylinder(surface_id=246, r=1.1125)

surf251 = openmc.ZCylinder(surface_id=251, r=1.1125)

surf256 = openmc.ZCylinder(surface_id=256, r=1.1125)

# 1st (inner) plate and slot
surf301 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=301)

# 2nd (inner) plate and slot
surf302 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=302)

# 3rd (inner) plate and slot
surf303 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=303)

# 4th (inner) plate and slot
surf304 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=304)

# 5th (inner) plate and slot
surf305 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=305)

# 6th (inner) plate and slot
surf306 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=306)

# 7th (inner) plate and slot
surf307 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=307)

# 8th (inner) plate and slot
surf308 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=308)

# 9th (inner) plate and slot
surf309 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=309)

# 10th (inner) plate and slot
surf310 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=310)

# 11th (inner) plate and slot
surf311 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=311)

# 12th (inner) plate and slot
surf312 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=312)

# 13th (inner) plate and slot
surf313 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=313)

# 14th (inner) plate and slot
surf314 = openmc.model.RectangularParallelepiped(
    -3.55, 3.55, -0.08, 0.08, -33.85, 39.65, surface_id=314)

# Slot for Cd wire: outer fuel plate
surf321 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -38.35, 39.65, surface_id=321)

# Slot for Cd wire: inner fuel plate
surf323 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 39.65, surface_id=323)

# Slot for Cd wire: inner fuel plate
surf325 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 39.65, surface_id=325)

# Slot for Cd wire: inner fuel plate
surf327 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 39.65, surface_id=327)

# Slot for Cd wire: inner fuel plate
surf329 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 39.65, surface_id=329)

# Slot for Cd wire: inner fuel plate
surf331 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 39.65, surface_id=331)

# Slot for Cd wire: inner fuel plate
surf333 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 39.65, surface_id=333)

# Slot for Cd wire: inner fuel plate
surf335 = openmc.model.RectangularParallelepiped(
    -0.03, 0.03, -0.025, 0.025, -33.85, 39.65, surface_id=335)

# Al6061 external guide plate slot
surf341 = openmc.model.RectangularParallelepiped(
    -3.555, 3.555, -0.065, 0.065, -38.35, 39.65, surface_id=341)

# Al6061 internal guide plate slot
surf342 = openmc.model.RectangularParallelepiped(
    -3.555, 3.555, -0.065, 0.065, -33.85, 39.65, surface_id=342)

# Al6061 internal guide plate slot
surf343 = openmc.model.RectangularParallelepiped(
    -3.555, 3.555, -0.065, 0.065, -33.85, 39.65, surface_id=343)

# Al6061 external guide plate slot
surf344 = openmc.model.RectangularParallelepiped(
    -3.555, 3.555, -0.065, 0.065, -38.35, 39.65, surface_id=344)

# Control element, fully extracted: upper
surf401 = openmc.model.RectangularParallelepiped(
    -3.245, 3.245, -0.215, 0.215, 30.45, 117.15, surface_id=401)

# Control element, fully extracted: lower
surf402 = openmc.model.RectangularParallelepiped(
    -3.245, 3.245, -0.215, 0.215, 30.45, 117.15, surface_id=402)

# Control element 4: upper
surf403 = openmc.model.RectangularParallelepiped(
    -3.245, 3.245, -0.215, 0.215, 30.45, 117.15, surface_id=403)

# Control element 4: lower
surf404 = openmc.model.RectangularParallelepiped(
    -3.245, 3.245, -0.215, 0.215, 30.45, 117.15, surface_id=404)

surf901 = openmc.model.RectangularParallelepiped(
    -15.4, 23.1, -24.3, 16.2, -180, 800, surface_id=901)

surf999 = openmc.ZCylinder(surface_id=999, r=120)


# ==============================================================================
# Universes (from COG define unit blocks)
# ==============================================================================

# Unit 1: Inner fuel plate and slot
u1_cell0 = openmc.Cell(fill=mat1, name="Al-U3Si2")
u1_cell0.region = -surf21
u1_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u1_cell1.region = +surf21 & -surf22
u1_cell2 = openmc.Cell(fill=mat3, name="H2O")
u1_cell2.region = +surf22 & -surf23
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

# Unit 2: Outer fuel plate and slot
u2_cell0 = openmc.Cell(fill=mat1, name="Al-U3Si2")
u2_cell0.region = -surf21
u2_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u2_cell1.region = +surf21 & -surf24
u2_cell2 = openmc.Cell(fill=mat3, name="H2O")
u2_cell2.region = +surf24 & -surf25
universe2 = openmc.Universe(universe_id=2, cells=[u2_cell0, u2_cell1, u2_cell2])

# Unit 3: Cd wire and slot
u3_cell0 = openmc.Cell(fill=mat4, name="Cd")
u3_cell0.region = -surf26
universe3 = openmc.Universe(universe_id=3, cells=[u3_cell0])

# Unit 4: Normal fuel element
u4_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u4_cell0.region = +surf41 & -surf42 & +surf23 & +surf101 & +surf102 & +surf103 & +surf104 & +surf105 & +surf106 & +surf107 & +surf108 & +surf109
u4_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u4_cell1.region = +surf41 & +surf42 & +surf30 & -surf31
universe4 = openmc.Universe(universe_id=4, cells=[u4_cell0, u4_cell1])

# Unit 5: Control fuel element
u5_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u5_cell0.region = +surf41 & -surf42 & +surf301 & +surf302 & +surf303 & +surf304 & +surf305 & +surf306 & +surf307 & +surf308 & +surf309 & +surf310
u5_cell1 = openmc.Cell(fill=mat2, name="Al6061")
u5_cell1.region = +surf41 & +surf42 & +surf30 & -surf31
universe5 = openmc.Universe(universe_id=5, cells=[u5_cell0, u5_cell1])

# Unit 6: Grid plate unit cell with no nozzle
u6_cell0 = openmc.Cell(fill=mat5, name="Al")
u6_cell0.region = -surf90 & +surf91
u6_cell1 = openmc.Cell(fill=mat3, name="H2O")
u6_cell1.region = -surf90 & -surf91
universe6 = openmc.Universe(universe_id=6, cells=[u6_cell0, u6_cell1])

# Unit 7: Grid plate unit cell with a nozzle
u7_cell0 = openmc.Cell(fill=mat5, name="Al")
u7_cell0.region = -surf90 & +surf91
universe7 = openmc.Universe(universe_id=7, cells=[u7_cell0])

# Unit 8: Grid plate
# Lattice for unit 8: 8x10 array
# Pitch: (7.700000, 8.100000) cm
# Lower left: (-30.8, -40.5)

# Creating 8x10 lattice with mixed universes
lattice8 = openmc.RectLattice(lattice_id=8)
lattice8.lower_left = [-30.8, -40.5]
lattice8.pitch = [7.700000, 8.100000]
# TODO: Set up universe array for mixed lattice
lattice8.universes = [[universe6]*8]*10
universe8 = lattice8  # Lattice can be used as universe

# Unit 9: Grid plate with primary and secondary holes
universe9 = openmc.Universe(universe_id=9, cells=[])

# Unit 10: Al6061 internal guide plate in water
u10_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u10_cell0.region = -surf43
universe10 = openmc.Universe(universe_id=10, cells=[u10_cell0])

# Unit 11: Al6061 External guide plate in water
u11_cell0 = openmc.Cell(fill=mat2, name="Al6061")
u11_cell0.region = -surf44
universe11 = openmc.Universe(universe_id=11, cells=[u11_cell0])

# Unit 12: Control element, fully withdrawn
u12_cell0 = openmc.Cell(fill=mat6, name="AgInCd")
u12_cell0.region = -surf50 & +surf51 & -surf52 & +surf55 & -surf56
u12_cell1 = openmc.Cell(fill=mat6, name="AgInCd")
u12_cell1.region = -surf50 & -surf51 & -surf53 & -surf56
u12_cell2 = openmc.Cell(fill=mat6, name="AgInCd")
u12_cell2.region = -surf50 & +surf52 & -surf54 & -surf56
u12_cell3 = openmc.Cell(fill=mat7, name="SS304L")
u12_cell3.region = +surf56 & -surf57
u12_cell4 = openmc.Cell(fill=mat3, name="H2O")
u12_cell4.region = +surf57 & -surf999
universe12 = openmc.Universe(universe_id=12, cells=[u12_cell0, u12_cell1, u12_cell2, u12_cell3, u12_cell4])

# Unit 13: Control fuel elements 1, 2 and 3: fully extracted
universe13 = openmc.Universe(universe_id=13, cells=[])

# Unit 14: Control fuel element 4: partially extracted
universe14 = openmc.Universe(universe_id=14, cells=[])

# Unit 15: Water
u15_cell0 = openmc.Cell(fill=mat3, name="H2O")
u15_cell0.region = -surf999
universe15 = openmc.Universe(universe_id=15, cells=[u15_cell0])

# Unit 16: Core: array of NFEs and CFEs
# Lattice for unit 16: 5x5 array
# Pitch: (7.700000, 8.100000) cm
# Lower left: (-15.4, -24.3)

# Creating 5x5 lattice with mixed universes
lattice16 = openmc.RectLattice(lattice_id=16)
lattice16.lower_left = [-15.4, -24.3]
lattice16.pitch = [7.700000, 8.100000]
# TODO: Set up universe array for mixed lattice
lattice16.universes = [[universe15]*5]*5
universe16 = lattice16  # Lattice can be used as universe

# Cell: H2O
cell0 = openmc.Cell(cell_id=0, fill=mat3, name="H2O")
cell0.region = -surf999 & +surf901 & +surf90 & +surf60

# Cell: H2O
cell1 = openmc.Cell(cell_id=1, fill=mat3, name="H2O")
cell1.region = +surf22 & -surf23

# Cell: H2O
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="H2O")
cell2.region = +surf24 & -surf25

# Cell: H2O
cell3 = openmc.Cell(cell_id=3, fill=mat3, name="H2O")
cell3.region = -surf90 & -surf91

# Cell: H2O
cell4 = openmc.Cell(cell_id=4, fill=mat3, name="H2O")
cell4.region = +surf57 & -surf999

# Cell: H2O
cell5 = openmc.Cell(cell_id=5, fill=mat3, name="H2O")
cell5.region = -surf999

# Cell: Cd
cell6 = openmc.Cell(cell_id=6, fill=mat4, name="Cd")
cell6.region = -surf999 & -surf60 & -surf61 & +surf62

# Cell: Cd
cell7 = openmc.Cell(cell_id=7, fill=mat4, name="Cd")
cell7.region = -surf999 & -surf60 & -surf63 & +surf64

# Cell: Al
cell8 = openmc.Cell(cell_id=8, fill=mat5, name="Al")
cell8.region = -surf999 & -surf60 & +surf61

# Cell: Al
cell9 = openmc.Cell(cell_id=9, fill=mat5, name="Al")
cell9.region = -surf999 & -surf60 & -surf62 & +surf63

# Cell: Al
cell10 = openmc.Cell(cell_id=10, fill=mat5, name="Al")
cell10.region = -surf90 & +surf91

# Cell: Al2O3
cell11 = openmc.Cell(cell_id=11, fill=mat8, name="Al2O3")
cell11.region = -surf999 & -surf60 & -surf64

# Cell using unit 9: GP
cell12 = openmc.Cell(cell_id=12, fill=universe9, name="GP")
cell12.region = -surf999 & -surf90 & +surf60

# Cell using unit 16: Core
cell13 = openmc.Cell(cell_id=13, fill=universe16, name="Core")
cell13.region = -surf999 & -surf901 & +surf90

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 1000
settings.batches = 5050
settings.inactive = 51
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
