"""
HMF072-3: ZEUS Config. No. 3 (Rev. 1)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# HEU: U1
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")

# HEU: U2 (case 2)
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")

# HEU: U3
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")

# HEU: U4
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")

# HEU: U5
mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")

# HEU: U5
mat6 = openmc.Material(material_id=6)
mat6.set_density("sum")

# HEU: L1
mat11 = openmc.Material(material_id=11)
mat11.set_density("sum")

# HEU: L2
mat12 = openmc.Material(material_id=12)
mat12.set_density("sum")

# HEU
mat13 = openmc.Material(material_id=13)
mat13.set_density("sum")

# HEU
mat14 = openmc.Material(material_id=14)
mat14.set_density("sum")

# HEU
mat15 = openmc.Material(material_id=15)
mat15.set_density("sum")

# Fe
mat21 = openmc.Material(material_id=21)
mat21.set_density("sum")

# Fe
mat22 = openmc.Material(material_id=22)
mat22.set_density("sum")

# Fe
mat23 = openmc.Material(material_id=23)
mat23.set_density("sum")

# Fe
mat24 = openmc.Material(material_id=24)
mat24.set_density("sum")

# Fe
mat25 = openmc.Material(material_id=25)
mat25.set_density("sum")

# Fe
mat26 = openmc.Material(material_id=26)
mat26.set_density("sum")

# Fe
mat27 = openmc.Material(material_id=27)
mat27.set_density("sum")

# Fe
mat28 = openmc.Material(material_id=28)
mat28.set_density("sum")

# Fe
mat29 = openmc.Material(material_id=29)
mat29.set_density("sum")

# Fe
mat30 = openmc.Material(material_id=30)
mat30.set_density("sum")

# Fe
mat31 = openmc.Material(material_id=31)
mat31.set_density("sum")

# Fe
mat32 = openmc.Material(material_id=32)
mat32.set_density("sum")

# Fe
mat41 = openmc.Material(material_id=41)
mat41.set_density("sum")

# Fe
mat42 = openmc.Material(material_id=42)
mat42.set_density("sum")

# Fe
mat43 = openmc.Material(material_id=43)
mat43.set_density("sum")

# Fe
mat44 = openmc.Material(material_id=44)
mat44.set_density("sum")

# Fe
mat45 = openmc.Material(material_id=45)
mat45.set_density("sum")

# Fe
mat46 = openmc.Material(material_id=46)
mat46.set_density("sum")

# Fe
mat47 = openmc.Material(material_id=47)
mat47.set_density("sum")

# Fe
mat48 = openmc.Material(material_id=48)
mat48.set_density("sum")

# Fe
mat49 = openmc.Material(material_id=49)
mat49.set_density("sum")

# Fe
mat50 = openmc.Material(material_id=50)
mat50.set_density("sum")

# CH2
mat60 = openmc.Material(material_id=60)
mat60.set_density("sum")
mat60.add_nuclide("H1", 8.071370e-02)
mat60.add_element("C", 4.035690e-02)
mat60.add_s_alpha_beta("c_H_in_CH2")

# Cu
mat61 = openmc.Material(material_id=61)
mat61.set_density("sum")

# Al
mat62 = openmc.Material(material_id=62)
mat62.set_density("sum")

# Al
mat63 = openmc.Material(material_id=63)
mat63.set_density("sum")

# C
mat64 = openmc.Material(material_id=64)
mat64.set_density("sum")

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat11, mat12, mat13, mat14, mat15, mat21, mat22, mat23, mat24, mat25, mat26, mat27, mat28, mat29, mat30, mat31, mat32, mat41, mat42, mat43, mat44, mat45, mat46, mat47, mat48, mat49, mat50, mat60, mat61, mat62, mat63, mat64])

# ==============================================================================
# Geometry
# ==============================================================================

# alignment tube inner radius
surf1 = openmc.ZCylinder(surface_id=1, r=2.54)
# alignment tube outer radius
surf2 = openmc.ZCylinder(surface_id=2, r=3.1496)
# lower reflector inner radius
surf3 = openmc.ZCylinder(surface_id=3, r=3.175)
# platen/adapter plate inner radius
surf4 = openmc.ZCylinder(surface_id=4, r=4.7625)
# inner HEU IR
surf5 = openmc.ZCylinder(surface_id=5, r=7.62)
# inner HEU OR/outer HEU IR
surf6 = openmc.ZCylinder(surface_id=6, r=19.05)
# outer HEU/OR
surf7 = openmc.ZCylinder(surface_id=7, r=26.67)
# corner reflector radius
surf8 = openmc.ZCylinder(surface_id=8, r=26.797)
surf11 = openmc.XPlane(surface_id=11, x0=-44.1452)
surf12 = openmc.XPlane(surface_id=12, x0=-27.940)
surf13 = openmc.XPlane(surface_id=13, x0=27.940)
surf14 = openmc.XPlane(surface_id=14, x0=44.1452)
surf16 = openmc.YPlane(surface_id=16, y0=-44.1452)
surf17 = openmc.YPlane(surface_id=17, y0=-27.940)
surf18 = openmc.YPlane(surface_id=18, y0=27.940)
surf19 = openmc.YPlane(surface_id=19, y0=44.1452)
# bottom of alignment tube (bottom of model)
surf20 = openmc.ZPlane(surface_id=20, z0=-3.2512)
# bottom of side & corner reflector
surf21 = openmc.ZPlane(surface_id=21, z0=2.540)
# top of side reflector (top of model)
surf22 = openmc.ZPlane(surface_id=22, z0=105.791)
# bottom of diaphragm/   top of lower corner reflector
surf23 = openmc.ZPlane(surface_id=23, z0=60.24880)
# top of diaphragm/bottom of upper corner reflector
surf24 = openmc.ZPlane(surface_id=24, z0=60.51296)
# top of upper corner reflector
surf25 = openmc.ZPlane(surface_id=25, z0=77.48016)
# bottom of platen/adapter plate
surf26 = openmc.ZPlane(surface_id=26, z0=25.603019)
# top of platen/adapter plate / start lower axial Cu
surf27 = openmc.ZPlane(surface_id=27, z0=31.953019)
# end lower axial Cu
surf28 = openmc.ZPlane(surface_id=28, z0=46.380219)
# bottom of upper axial void / end
surf29 = openmc.ZPlane(surface_id=29, z0=77.420304)
# end upper axial Cu
surf30 = openmc.ZPlane(surface_id=30, z0=91.907736)
# start lp1
surf51 = openmc.ZPlane(surface_id=51, z0=46.380219)
# end lp1   /start l1-l
surf52 = openmc.ZPlane(surface_id=52, z0=46.486391)
# end l1-l  /start l1-heu & l1-al spacer
surf53 = openmc.ZPlane(surface_id=53, z0=47.678923)
# end l1-heu & l1-al spacer/start l1-u
surf54 = openmc.ZPlane(surface_id=54, z0=47.978643)
# end l1-u  /start lp2
surf55 = openmc.ZPlane(surface_id=55, z0=49.161859)
# end lp2   /start l2-l
surf56 = openmc.ZPlane(surface_id=56, z0=49.268031)
# end l2-l  /start l2-heu
surf57 = openmc.ZPlane(surface_id=57, z0=50.451247)
# end l2-heu/start l2-u
surf58 = openmc.ZPlane(surface_id=58, z0=50.750967)
# end l2-u  /start lp3
surf59 = openmc.ZPlane(surface_id=59, z0=51.933866)
# end lp3   /start l3-l
surf60 = openmc.ZPlane(surface_id=60, z0=52.040038)
# end l3-l  /start l3-heu
surf61 = openmc.ZPlane(surface_id=61, z0=53.222727)
# end l3-heu/start l3-u
surf62 = openmc.ZPlane(surface_id=62, z0=53.522447)
# end l3-u  /start lp4
surf63 = openmc.ZPlane(surface_id=63, z0=54.707359)
# end lp4   /start l4-l
surf64 = openmc.ZPlane(surface_id=64, z0=54.813531)
# end l4-l  /start l4-heu
surf65 = openmc.ZPlane(surface_id=65, z0=55.997875)
# end l4-heu/start l4-u
surf66 = openmc.ZPlane(surface_id=66, z0=56.297595)
# end l4-u  /start lp5
surf67 = openmc.ZPlane(surface_id=67, z0=57.479966)
# end lp5   /start l5-l
surf68 = openmc.ZPlane(surface_id=68, z0=57.586138)
# end l5-l  /start l5-heu
surf69 = openmc.ZPlane(surface_id=69, z0=58.768720)
# end l5-heu/start l5-u
surf70 = openmc.ZPlane(surface_id=70, z0=59.068440)
# end l5-u
surf71 = openmc.ZPlane(surface_id=71, z0=60.248800)
# start up1
surf72 = openmc.ZPlane(surface_id=72, z0=60.512960)
# end up1   /start u1-l
surf73 = openmc.ZPlane(surface_id=73, z0=60.619132)
# end u1-l  /start u1-heu
surf74 = openmc.ZPlane(surface_id=74, z0=61.814624)
# end u1-heu/start u1-u
surf75 = openmc.ZPlane(surface_id=75, z0=62.114344)
# end u1-u  /start up2
surf76 = openmc.ZPlane(surface_id=76, z0=63.315340)
# end up2   /start u2-l
surf77 = openmc.ZPlane(surface_id=77, z0=63.421512)
# end u2-l  /start u2-heu
surf78 = openmc.ZPlane(surface_id=78, z0=64.612348)
# end u2-heu/start u2-u
surf79 = openmc.ZPlane(surface_id=79, z0=64.912068)
# end u2-u  /start up3
surf80 = openmc.ZPlane(surface_id=80, z0=66.105444)
# end up3   /start u3-l
surf81 = openmc.ZPlane(surface_id=81, z0=66.211616)
# end u3-l  /start u3-heu
surf82 = openmc.ZPlane(surface_id=82, z0=67.401608)
# end u3-heu/start u3-u
surf83 = openmc.ZPlane(surface_id=83, z0=67.701328)
# end u3-u  /start up4
surf84 = openmc.ZPlane(surface_id=84, z0=68.893860)
# end up4   /start u4-l
surf85 = openmc.ZPlane(surface_id=85, z0=69.000032)
# end u4-l  /start u4-heu
surf86 = openmc.ZPlane(surface_id=86, z0=70.191292)
# end u4-heu/start u4-u
surf87 = openmc.ZPlane(surface_id=87, z0=70.491012)
# end u4-u  /start up5
surf88 = openmc.ZPlane(surface_id=88, z0=71.717408)
# end up5   /start u5-l
surf89 = openmc.ZPlane(surface_id=89, z0=71.823580)
# end u5-l  /start u5-heu
surf90 = openmc.ZPlane(surface_id=90, z0=73.025424)
# end u5-heu/start u5-u
surf91 = openmc.ZPlane(surface_id=91, z0=73.325144)
# end u5-u  /start up6
surf92 = openmc.ZPlane(surface_id=92, z0=74.532068)
# end up6   /start u6-l
surf93 = openmc.ZPlane(surface_id=93, z0=74.638240)
# end u6-l  /start u6-heu
surf94 = openmc.ZPlane(surface_id=94, z0=75.826960)
# end u6-heu/start u6-u
surf95 = openmc.ZPlane(surface_id=95, z0=76.126680)
# end u6-u  /start up7
surf96 = openmc.ZPlane(surface_id=96, z0=77.314132)
# end up7
surf97 = openmc.ZPlane(surface_id=97, z0=77.420304)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# MAT6
cell1 = openmc.Cell(cell_id=1, fill=mat6)
cell1.region = -surf95 & +surf94 & -surf7

# MAT5
cell2 = openmc.Cell(cell_id=2, fill=mat5)
cell2.region = -surf91 & +surf90 & -surf7

# MAT4
cell3 = openmc.Cell(cell_id=3, fill=mat4)
cell3.region = -surf87 & +surf86 & -surf7

# MAT3
cell4 = openmc.Cell(cell_id=4, fill=mat3)
cell4.region = -surf83 & +surf82 & -surf7

# MAT2
cell5 = openmc.Cell(cell_id=5, fill=mat2)
cell5.region = -surf79 & +surf78 & -surf7

# MAT1
cell6 = openmc.Cell(cell_id=6, fill=mat1)
cell6.region = -surf75 & +surf74 & -surf7

# MAT15
cell7 = openmc.Cell(cell_id=7, fill=mat15)
cell7.region = -surf70 & +surf69 & -surf7 & +surf3

# MAT14
cell8 = openmc.Cell(cell_id=8, fill=mat14)
cell8.region = -surf66 & +surf65 & -surf7 & +surf3

# MAT13
cell9 = openmc.Cell(cell_id=9, fill=mat13)
cell9.region = -surf62 & +surf61 & -surf7 & +surf3

# MAT12
cell10 = openmc.Cell(cell_id=10, fill=mat12)
cell10.region = -surf58 & +surf57 & -surf7 & +surf3

# MAT11
cell11 = openmc.Cell(cell_id=11, fill=mat11)
cell11.region = -surf54 & +surf53 & -surf7 & +surf5

# MAT63
cell12 = openmc.Cell(cell_id=12, fill=mat63)
cell12.region = -surf54 & +surf53 & -surf5 & +surf3

# MAT31
cell13 = openmc.Cell(cell_id=13, fill=mat31)
cell13.region = -surf96 & +surf95 & -surf7

# MAT32
cell14 = openmc.Cell(cell_id=14, fill=mat32)
cell14.region = -surf94 & +surf93 & -surf7

# MAT29
cell15 = openmc.Cell(cell_id=15, fill=mat29)
cell15.region = -surf92 & +surf91 & -surf7

# MAT30
cell16 = openmc.Cell(cell_id=16, fill=mat30)
cell16.region = -surf90 & +surf89 & -surf7

# MAT27
cell17 = openmc.Cell(cell_id=17, fill=mat27)
cell17.region = -surf88 & +surf87 & -surf7

# MAT28
cell18 = openmc.Cell(cell_id=18, fill=mat28)
cell18.region = -surf86 & +surf85 & -surf7

# MAT25
cell19 = openmc.Cell(cell_id=19, fill=mat25)
cell19.region = -surf84 & +surf83 & -surf7

# MAT26
cell20 = openmc.Cell(cell_id=20, fill=mat26)
cell20.region = -surf82 & +surf81 & -surf7

# MAT23
cell21 = openmc.Cell(cell_id=21, fill=mat23)
cell21.region = -surf80 & +surf79 & -surf7

# MAT24
cell22 = openmc.Cell(cell_id=22, fill=mat24)
cell22.region = -surf78 & +surf77 & -surf7

# MAT21
cell23 = openmc.Cell(cell_id=23, fill=mat21)
cell23.region = -surf76 & +surf75 & -surf7

# MAT22
cell24 = openmc.Cell(cell_id=24, fill=mat22)
cell24.region = -surf74 & +surf73 & -surf7

# MAT49
cell25 = openmc.Cell(cell_id=25, fill=mat49)
cell25.region = -surf71 & +surf70 & -surf7 & +surf3

# MAT50
cell26 = openmc.Cell(cell_id=26, fill=mat50)
cell26.region = -surf69 & +surf68 & -surf7 & +surf3

# MAT47
cell27 = openmc.Cell(cell_id=27, fill=mat47)
cell27.region = -surf67 & +surf66 & -surf7 & +surf3

# MAT48
cell28 = openmc.Cell(cell_id=28, fill=mat48)
cell28.region = -surf65 & +surf64 & -surf7 & +surf3

# MAT45
cell29 = openmc.Cell(cell_id=29, fill=mat45)
cell29.region = -surf63 & +surf62 & -surf7 & +surf3

# MAT46
cell30 = openmc.Cell(cell_id=30, fill=mat46)
cell30.region = -surf61 & +surf60 & -surf7 & +surf3

# MAT43
cell31 = openmc.Cell(cell_id=31, fill=mat43)
cell31.region = -surf59 & +surf58 & -surf7 & +surf3

# MAT44
cell32 = openmc.Cell(cell_id=32, fill=mat44)
cell32.region = -surf57 & +surf56 & -surf7 & +surf3

# MAT41
cell33 = openmc.Cell(cell_id=33, fill=mat41)
cell33.region = -surf55 & +surf54 & -surf7 & +surf3

# MAT42
cell34 = openmc.Cell(cell_id=34, fill=mat42)
cell34.region = -surf53 & +surf52 & -surf7 & +surf3

# POLY
cell35 = openmc.Cell(cell_id=35, fill=mat60)
cell35.region = -surf97 & +surf96 & -surf7 & -surf93 & +surf97 & -surf7

# MAT61
cell36 = openmc.Cell(cell_id=36, fill=mat61)
cell36.region = +surf11 & -surf12 & +surf16 & -surf19 & +surf21 & -surf22 & +surf12 & -surf13 & +surf16 & -surf17 & +surf21 & -surf22

# MAT61
cell37 = openmc.Cell(cell_id=37, fill=mat61)
cell37.region = +surf12 & -surf13 & +surf17 & -surf18 & +surf21 & -surf23 & +surf8

# MAT61
cell38 = openmc.Cell(cell_id=38, fill=mat61)
cell38.region = +surf12 & -surf13 & +surf17 & -surf18 & +surf24 & -surf25 & +surf8

# MAT61
cell39 = openmc.Cell(cell_id=39, fill=mat61)
cell39.region = +surf12 & -surf13 & +surf17 & -surf18 & -surf30 & +surf25

# MAT64
cell40 = openmc.Cell(cell_id=40, fill=mat64)
cell40.region = +surf12 & -surf13 & +surf17 & -surf18 & -surf24 & +surf23

# MAT62
cell41 = openmc.Cell(cell_id=41, fill=mat62)
cell41.region = -surf23 & +surf20 & -surf2 & +surf1

# MAT61
cell42 = openmc.Cell(cell_id=42, fill=mat61)
cell42.region = -surf28 & +surf27 & -surf7 & +surf3

# MAT62
cell43 = openmc.Cell(cell_id=43, fill=mat62)
cell43.region = -surf27 & +surf26 & -surf7 & +surf4

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24, cell25, cell26, cell27, cell28, cell29, cell30, cell31, cell32, cell33, cell34, cell35, cell36, cell37, cell38, cell39, cell40, cell41, cell42, cell43])
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
source.space = openmc.stats.Box((-7.0, -7.0, 46.83), (7.0, 7.0, 76.98))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
