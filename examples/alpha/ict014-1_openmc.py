"""
ICT014-1: RA-6
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Al-U3Si2
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 1.768000e-05)
mat1.add_nuclide("U235", 2.417000e-03)
mat1.add_nuclide("U236", 1.363500e-05)
mat1.add_nuclide("U238", 9.653800e-03)
mat1.add_element("Si", 8.362200e-03)
mat1.add_element("Al", 3.045100e-02)
mat1.add_element("B", 1.929600e-06)

# Al6061
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("Al", 5.881100e-02)
mat2.add_element("Cu", 6.652700e-05)
mat2.add_element("Cr", 3.439800e-05)
mat2.add_element("Mg", 6.622900e-04)
mat2.add_element("Si", 3.878900e-04)
mat2.add_element("Zn", 2.486600e-07)
mat2.add_element("Fe", 7.861000e-05)
mat2.add_element("Mn", 2.663700e-05)
mat2.add_element("Ti", 3.395900e-06)
mat2.add_element("B", 1.504000e-06)
mat2.add_element("Cd", 7.232300e-08)
mat2.add_element("Co", 1.655400e-06)
mat2.add_element("Li", 2.342600e-07)

# H2O
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("H1", 6.667500e-02)
mat3.add_nuclide("O16", 3.333800e-02)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Cd
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_element("Cd", 4.634000e-02)

# Al
mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")
mat5.add_element("Al", 6.026200e-02)

# Ag-In-Cd
mat6 = openmc.Material(material_id=6)
mat6.set_density("sum")
mat6.add_element("Ag", 4.536500e-02)
mat6.add_element("In", 7.976500e-03)
mat6.add_element("Cd", 2.833100e-03)

# Stainless steel
mat7 = openmc.Material(material_id=7)
mat7.set_density("sum")
mat7.add_element("Fe", 5.989900e-02)
mat7.add_element("Ni", 8.198400e-03)
mat7.add_element("Cr", 1.758200e-02)
mat7.add_element("Si", 6.424600e-04)
mat7.add_element("Mn", 8.758300e-04)
mat7.add_element("C", 6.009100e-05)
mat7.add_element("S", 2.250500e-05)
mat7.add_element("P", 3.495300e-05)
mat7.add_element("N", 1.717600e-04)

# Al2O3
mat8 = openmc.Material(material_id=8)
mat8.set_density("sum")
mat8.add_element("Al", 4.677800e-02)
mat8.add_nuclide("O16", 7.016600e-02)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8])

# ==============================================================================
# Geometry
# ==============================================================================

# Al6061 side plates: inner
surf11 = openmc.model.RectangularParallelepiped(-3.3, 3.3, -4, 4, -38.35, 38.35)
# Al6061 side plates: outer
surf12 = openmc.model.RectangularParallelepiped(-3.8, 3.8, -4, 4, -38.35, 38.35)
# Al-U3Si2 fuel
surf21 = openmc.model.RectangularParallelepiped(-3, 3, -0.02533, 0.02533, -31.15, 31.15)
# Al6061 clad: inner fuel plates
surf22 = openmc.model.RectangularParallelepiped(-3.525, 3.525, -0.0745, 0.0745, -33.55, 33.55)
# Water slot:  inner fuel plates
surf23 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# Al6061 clad: outer fuel plates
surf24 = openmc.model.RectangularParallelepiped(-3.525, 3.525, -0.0745, 0.0745, -38.35, 35.15)
# Water slot:  outer fuel plates
surf25 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -38.35, 35.15)
# Cd wire
surf26 = openmc.ZCylinder(surface_id=26, r=0.02425)
# Al6061 crossbar
surf29 = openmc.XCylinder(surface_id=29, x0=-3.3, y0=3.3, r=0.625)
# Al6061 nozzle: inner
surf30 = openmc.ZCylinder(surface_id=30, r=2.4895)
# Al6061 nozzle: outer
surf31 = openmc.ZCylinder(surface_id=31, r=3.0895)
# Al6061 side plates: inner
surf41 = openmc.model.RectangularParallelepiped(-3.3, 3.3, -4, 4, -38.35, 39.65)
# Al6061 side plates: outer
surf42 = openmc.model.RectangularParallelepiped(-3.8, 3.8, -4, 4, -38.35, 39.65)
# Al6061 internal guide plate
surf43 = openmc.model.RectangularParallelepiped(-3.525, 3.525, -0.065, 0.065, -33.55, 33.55)
# Al6061 external guide plate
surf44 = openmc.model.RectangularParallelepiped(-3.525, 3.525, -0.065, 0.065, -38.35, 39.65)
# In-Ag-Cd
surf50 = openmc.model.RectangularParallelepiped(-3.09, 3.09, -0.11, 0.11, 31.55, 94.95)
# edge
surf51 = openmc.XPlane(surface_id=51, x0=-2.98)
# edge
surf52 = openmc.XPlane(surface_id=52, x0=2.98)
# edge
surf53 = openmc.ZCylinder(surface_id=53, x0=tr, y0=-2.98, r=0.11)
# edge
surf54 = openmc.ZCylinder(surface_id=54, x0=tr, y0=2.98, r=0.11)
# indentation
# surf55: Unsupported surface type "s" with params ['0.1625', 'tr', '0', '0.1725', '63.25']
# SS304L cladding: inner
surf56 = openmc.model.RectangularParallelepiped(-3.175, 3.175, -0.145, 0.145, 31.55, 94.95)
# SS304L cladding: outer
surf57 = openmc.model.RectangularParallelepiped(-3.245, 3.245, -0.215, 0.215, 30.45, 117.15)
# BNCT
surf60 = openmc.model.RectangularParallelepiped(-38.55, 38.55, -199.0, -32.40, -40.75, 41.60)
# Al/Cd
surf61 = openmc.YPlane(surface_id=61, y0=-49.50)
# Cd/Al
surf62 = openmc.YPlane(surface_id=62, y0=-49.55)
# Al/Cd
surf63 = openmc.YPlane(surface_id=63, y0=-59.55)
# Cd/Al2O3
surf64 = openmc.YPlane(surface_id=64, y0=-59.70)
# Pure aluminum grid plate
surf90 = openmc.model.RectangularParallelepiped(-30.8, 30.8, -40.5, 40.5, -60.75, -40.75)
# Primary water hole
surf91 = openmc.ZCylinder(surface_id=91, r=3.0895)
# 1st (outer) plate and slot
surf101 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -38.35, 35.15)
# 2nd (inner) plate and slot
surf102 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 3rd (inner) plate and slot
surf103 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 4th (inner) plate and slot
surf104 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 5th (inner) plate and slot
surf105 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 6th (inner) plate and slot
surf106 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 7th (inner) plate and slot
surf107 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 8th (inner) plate and slot
surf108 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 9th (inner) plate and slot
surf109 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 11th (inner) plate and slot
surf111 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 12th (inner) plate and slot
surf112 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 13th (inner) plate and slot
surf113 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 14th (inner) plate and slot
surf114 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 15th (inner) plate and slot
surf115 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 16th (inner) plate and slot
surf116 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 17th (inner) plate and slot
surf117 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 18th (inner) plate and slot
surf118 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 35.15)
# 19th (outer) plate and slot
surf119 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -38.35, 35.15)
# Slot for Cd wire: outer fuel plate
surf121 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -38.35, 31.65)
# Slot for Cd wire: inner fuel plate
surf123 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: inner fuel plate
surf125 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: inner fuel plate
surf127 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: inner fuel plate
surf129 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: inner fuel plate
surf131 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: inner fuel plate
surf133 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: inner fuel plate
surf135 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: inner fuel plate
surf137 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 31.65)
# Slot for Cd wire: outer fuel plate
surf139 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -38.35, 31.65)
surf201 = openmc.ZCylinder(surface_id=201, x0=tr, y0=-15.4, r=1.1125)
surf206 = openmc.ZCylinder(surface_id=206, x0=tr, y0=-23.1, r=1.1125)
surf211 = openmc.ZCylinder(surface_id=211, x0=tr, y0=15.4, r=1.1125)
surf216 = openmc.ZCylinder(surface_id=216, x0=tr, y0=0.0, r=1.1125)
surf221 = openmc.ZCylinder(surface_id=221, x0=tr, y0=-15.4, r=1.1125)
surf226 = openmc.ZCylinder(surface_id=226, x0=tr, y0=23.1, r=1.1125)
surf231 = openmc.ZCylinder(surface_id=231, x0=tr, y0=7.7, r=1.1125)
surf236 = openmc.ZCylinder(surface_id=236, x0=tr, y0=-7.7, r=1.1125)
surf241 = openmc.ZCylinder(surface_id=241, x0=tr, y0=-23.1, r=1.1125)
surf246 = openmc.ZCylinder(surface_id=246, x0=tr, y0=15.4, r=1.1125)
surf251 = openmc.ZCylinder(surface_id=251, x0=tr, y0=0.0, r=1.1125)
surf256 = openmc.ZCylinder(surface_id=256, x0=tr, y0=-7.7, r=1.1125)
# 1st (inner) plate and slot
surf301 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 2nd (inner) plate and slot
surf302 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 3rd (inner) plate and slot
surf303 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 4th (inner) plate and slot
surf304 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 5th (inner) plate and slot
surf305 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 6th (inner) plate and slot
surf306 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 7th (inner) plate and slot
surf307 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 8th (inner) plate and slot
surf308 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 9th (inner) plate and slot
surf309 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 10th (inner) plate and slot
surf310 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 11th (inner) plate and slot
surf311 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 12th (inner) plate and slot
surf312 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 13th (inner) plate and slot
surf313 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# 14th (inner) plate and slot
surf314 = openmc.model.RectangularParallelepiped(-3.55, 3.55, -0.08, 0.08, -33.85, 39.65)
# Slot for Cd wire: outer fuel plate
surf321 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -38.35, 39.65)
# Slot for Cd wire: inner fuel plate
surf323 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 39.65)
# Slot for Cd wire: inner fuel plate
surf325 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 39.65)
# Slot for Cd wire: inner fuel plate
surf327 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 39.65)
# Slot for Cd wire: inner fuel plate
surf329 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 39.65)
# Slot for Cd wire: inner fuel plate
surf331 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 39.65)
# Slot for Cd wire: inner fuel plate
surf333 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 39.65)
# Slot for Cd wire: inner fuel plate
surf335 = openmc.model.RectangularParallelepiped(-0.03, 0.03, -0.025, 0.025, -33.85, 39.65)
# Al6061 external guide plate slot
surf341 = openmc.model.RectangularParallelepiped(-3.555, 3.555, -0.065, 0.065, -38.35, 39.65)
# Al6061 internal guide plate slot
surf342 = openmc.model.RectangularParallelepiped(-3.555, 3.555, -0.065, 0.065, -33.85, 39.65)
# Al6061 internal guide plate slot
surf343 = openmc.model.RectangularParallelepiped(-3.555, 3.555, -0.065, 0.065, -33.85, 39.65)
# Al6061 external guide plate slot
surf344 = openmc.model.RectangularParallelepiped(-3.555, 3.555, -0.065, 0.065, -38.35, 39.65)
# Control element, fully extracted: upper
surf401 = openmc.model.RectangularParallelepiped(-3.245, 3.245, -0.215, 0.215, 30.45, 117.15)
# Control element, fully extracted: lower
surf402 = openmc.model.RectangularParallelepiped(-3.245, 3.245, -0.215, 0.215, 30.45, 117.15)
# Control element 4: upper
surf403 = openmc.model.RectangularParallelepiped(-3.245, 3.245, -0.215, 0.215, 30.45, 117.15)
# Control element 4: lower
surf404 = openmc.model.RectangularParallelepiped(-3.245, 3.245, -0.215, 0.215, 30.45, 117.15)
surf901 = openmc.model.RectangularParallelepiped(-15.4, 23.1, -24.3, 16.2, -180, 800)
surf999 = openmc.ZCylinder(surface_id=999, r=120.0, boundary_type="vacuum")

# ------------------------------------------------------------------------------
# Universes
# ------------------------------------------------------------------------------

u1_cell0 = openmc.Cell(fill=mat1)
u1_cell0.region = -surf21
u1_cell1 = openmc.Cell(fill=mat2)
u1_cell1.region = +surf21 & -surf22
u1_cell2 = openmc.Cell(fill=mat3)
u1_cell2.region = +surf22 & -surf23
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

u2_cell0 = openmc.Cell(fill=mat1)
u2_cell0.region = -surf21
u2_cell1 = openmc.Cell(fill=mat2)
u2_cell1.region = +surf21 & -surf24
u2_cell2 = openmc.Cell(fill=mat3)
u2_cell2.region = +surf24 & -surf25
universe2 = openmc.Universe(universe_id=2, cells=[u2_cell0, u2_cell1, u2_cell2])

u3_cell0 = openmc.Cell(fill=mat4)
u3_cell0.region = -surf26
universe3 = openmc.Universe(universe_id=3, cells=[u3_cell0])

u4_cell0 = openmc.Cell(fill=mat2)
u4_cell0.region = +surf41 & -surf42 & +surf23 & +surf101 & +surf102 & +surf103 & +surf104 & +surf105 & +surf106 & +surf107 & +surf108 & +surf109
u4_cell1 = openmc.Cell(fill=mat2)
u4_cell1.region = +surf41 & +surf42 & +surf30 & -surf31
universe4 = openmc.Universe(universe_id=4, cells=[u4_cell0, u4_cell1])

u5_cell0 = openmc.Cell(fill=mat2)
u5_cell0.region = +surf41 & -surf42 & +surf301 & +surf302 & +surf303 & +surf304 & +surf305 & +surf306 & +surf307 & +surf308 & +surf309 & +surf310
u5_cell1 = openmc.Cell(fill=mat2)
u5_cell1.region = +surf41 & +surf42 & +surf30 & -surf31
universe5 = openmc.Universe(universe_id=5, cells=[u5_cell0, u5_cell1])

u6_cell0 = openmc.Cell(fill=mat5)
u6_cell0.region = -surf90 & +surf91
u6_cell1 = openmc.Cell(fill=mat3)
u6_cell1.region = -surf90 & -surf91
universe6 = openmc.Universe(universe_id=6, cells=[u6_cell0, u6_cell1])

u7_cell0 = openmc.Cell(fill=mat5)
u7_cell0.region = -surf90 & +surf91
universe7 = openmc.Universe(universe_id=7, cells=[u7_cell0])

# Lattice 8: 8x10 array
lattice8 = openmc.RectLattice(lattice_id=8)
lattice8.lower_left = [-30.8, -40.5]
lattice8.pitch = [7.700000, 8.100000]
lattice8.universes = [
    [universe6, universe6, universe6, universe6, universe6, universe6, universe6, universe6],
    [universe6, universe6, universe6, universe6, universe6, universe6, universe6, universe6],
    [universe6, universe6, universe6, universe7, universe7, universe7, universe6, universe6],
    [universe6, universe6, universe7, universe7, universe7, universe7, universe7, universe6],
    [universe6, universe6, universe7, universe7, universe7, universe7, universe7, universe6],
    [universe6, universe6, universe7, universe7, universe7, universe7, universe7, universe6],
    [universe6, universe6, universe6, universe7, universe7, universe6, universe6, universe6],
    [universe6, universe6, universe6, universe6, universe6, universe6, universe6, universe6],
    [universe6, universe6, universe6, universe6, universe6, universe6, universe6, universe6],
    [universe6, universe6, universe6, universe6, universe6, universe6, universe6, universe6],
]
universe8 = openmc.Universe(universe_id=8)
universe8.add_cell(openmc.Cell(fill=lattice8))

universe9 = openmc.Universe(universe_id=9, cells=[])

u10_cell0 = openmc.Cell(fill=mat2)
u10_cell0.region = -surf43
universe10 = openmc.Universe(universe_id=10, cells=[u10_cell0])

u11_cell0 = openmc.Cell(fill=mat2)
u11_cell0.region = -surf44
universe11 = openmc.Universe(universe_id=11, cells=[u11_cell0])

u12_cell0 = openmc.Cell(fill=mat6)
u12_cell0.region = -surf50 & +surf51 & -surf52 & +surf55 & -surf56
u12_cell1 = openmc.Cell(fill=mat6)
u12_cell1.region = -surf50 & -surf51 & -surf53 & -surf56
u12_cell2 = openmc.Cell(fill=mat6)
u12_cell2.region = -surf50 & +surf52 & -surf54 & -surf56
u12_cell3 = openmc.Cell(fill=mat7)
u12_cell3.region = +surf56 & -surf57
u12_cell4 = openmc.Cell(fill=mat3)
u12_cell4.region = +surf57 & -surf999
universe12 = openmc.Universe(universe_id=12, cells=[u12_cell0, u12_cell1, u12_cell2, u12_cell3, u12_cell4])

universe13 = openmc.Universe(universe_id=13, cells=[])

universe14 = openmc.Universe(universe_id=14, cells=[])

u15_cell0 = openmc.Cell(fill=mat3)
u15_cell0.region = -surf999
universe15 = openmc.Universe(universe_id=15, cells=[u15_cell0])

# Lattice 16: 5x5 array
lattice16 = openmc.RectLattice(lattice_id=16)
lattice16.lower_left = [-15.4, -24.3]
lattice16.pitch = [7.700000, 8.100000]
lattice16.universes = [
    [universe15, universe4, universe4, universe4, universe15],
    [universe4, universe13, universe4, universe13, universe4],
    [universe4, universe4, universe4, universe4, universe4],
    [universe4, universe13, universe4, universe14, universe4],
    [universe15, universe4, universe4, universe15, universe15],
]
universe16 = openmc.Universe(universe_id=16)
universe16.add_cell(openmc.Cell(fill=lattice16))

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Core
cell1 = openmc.Cell(cell_id=1, fill=universe16)
cell1.region = -surf999 & -surf901 & +surf90

# GP
cell2 = openmc.Cell(cell_id=2, fill=universe9)
cell2.region = -surf999 & -surf90 & +surf60

# Al
cell3 = openmc.Cell(cell_id=3, fill=mat5)
cell3.region = -surf999 & -surf60 & +surf61

# Cd
cell4 = openmc.Cell(cell_id=4, fill=mat4)
cell4.region = -surf999 & -surf60 & -surf61 & +surf62

# Al
cell5 = openmc.Cell(cell_id=5, fill=mat5)
cell5.region = -surf999 & -surf60 & -surf62 & +surf63

# Cd
cell6 = openmc.Cell(cell_id=6, fill=mat4)
cell6.region = -surf999 & -surf60 & -surf63 & +surf64

# Al2O3
cell7 = openmc.Cell(cell_id=7, fill=mat8)
cell7.region = -surf999 & -surf60 & -surf64

# H2O
cell8 = openmc.Cell(cell_id=8, fill=mat3)
cell8.region = -surf999 & +surf901 & +surf90 & +surf60

# H2O
cell12 = openmc.Cell(cell_id=12, fill=mat3)
cell12.region = +surf22 & -surf23

# H2O
cell16 = openmc.Cell(cell_id=16, fill=mat3)
cell16.region = +surf24 & -surf25

# H2O
cell24 = openmc.Cell(cell_id=24, fill=mat3)
cell24.region = -surf90 & -surf91

# Al
cell26 = openmc.Cell(cell_id=26, fill=mat5)
cell26.region = -surf90 & +surf91

# H2O
cell34 = openmc.Cell(cell_id=34, fill=mat3)
cell34.region = +surf57 & -surf999

# H2O
cell36 = openmc.Cell(cell_id=36, fill=mat3)
cell36.region = -surf999

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell12, cell16, cell24, cell26, cell34, cell36])
geometry = openmc.Geometry(root_universe)

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 1000
settings.batches = 5050
settings.inactive = 51
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
