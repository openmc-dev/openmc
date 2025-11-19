"""
HMF072-2: Zeus Iron Core Benchmark (No Al Shim)
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1, name="")
mat1.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat2 = openmc.Material(material_id=2, name="")
mat2.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat3 = openmc.Material(material_id=3, name="")
mat3.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat10 = openmc.Material(material_id=10, name="")
mat10.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat11 = openmc.Material(material_id=11, name="")
mat11.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat20 = openmc.Material(material_id=20, name="")
mat20.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat500 = openmc.Material(material_id=500, name="")
mat500.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat502 = openmc.Material(material_id=502, name="")
mat502.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat504 = openmc.Material(material_id=504, name="")
mat504.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat506 = openmc.Material(material_id=506, name="")
mat506.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat508 = openmc.Material(material_id=508, name="")
mat508.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat510 = openmc.Material(material_id=510, name="")
mat510.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat513 = openmc.Material(material_id=513, name="")
mat513.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat515 = openmc.Material(material_id=515, name="")
mat515.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat517 = openmc.Material(material_id=517, name="")
mat517.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat519 = openmc.Material(material_id=519, name="")
mat519.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat521 = openmc.Material(material_id=521, name="")
mat521.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat523 = openmc.Material(material_id=523, name="")
mat523.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat525 = openmc.Material(material_id=525, name="")
mat525.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat527 = openmc.Material(material_id=527, name="")
mat527.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat529 = openmc.Material(material_id=529, name="")
mat529.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat531 = openmc.Material(material_id=531, name="")
mat531.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4001 = openmc.Material(material_id=4001, name="")
mat4001.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4005 = openmc.Material(material_id=4005, name="")
mat4005.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4009 = openmc.Material(material_id=4009, name="")
mat4009.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4013 = openmc.Material(material_id=4013, name="")
mat4013.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4017 = openmc.Material(material_id=4017, name="")
mat4017.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4021 = openmc.Material(material_id=4021, name="")
mat4021.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4025 = openmc.Material(material_id=4025, name="")
mat4025.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4029 = openmc.Material(material_id=4029, name="")
mat4029.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4033 = openmc.Material(material_id=4033, name="")
mat4033.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4037 = openmc.Material(material_id=4037, name="")
mat4037.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4041 = openmc.Material(material_id=4041, name="")
mat4041.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4045 = openmc.Material(material_id=4045, name="")
mat4045.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4049 = openmc.Material(material_id=4049, name="")
mat4049.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4053 = openmc.Material(material_id=4053, name="")
mat4053.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4057 = openmc.Material(material_id=4057, name="")
mat4057.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4061 = openmc.Material(material_id=4061, name="")
mat4061.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4065 = openmc.Material(material_id=4065, name="")
mat4065.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4069 = openmc.Material(material_id=4069, name="")
mat4069.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4073 = openmc.Material(material_id=4073, name="")
mat4073.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4077 = openmc.Material(material_id=4077, name="")
mat4077.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4081 = openmc.Material(material_id=4081, name="")
mat4081.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4085 = openmc.Material(material_id=4085, name="")
mat4085.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4089 = openmc.Material(material_id=4089, name="")
mat4089.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4093 = openmc.Material(material_id=4093, name="")
mat4093.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4097 = openmc.Material(material_id=4097, name="")
mat4097.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4101 = openmc.Material(material_id=4101, name="")
mat4101.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4105 = openmc.Material(material_id=4105, name="")
mat4105.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4109 = openmc.Material(material_id=4109, name="")
mat4109.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4113 = openmc.Material(material_id=4113, name="")
mat4113.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4117 = openmc.Material(material_id=4117, name="")
mat4117.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4121 = openmc.Material(material_id=4121, name="")
mat4121.set_density("atom/b-cm", 1.0)  # TODO: Verify density

mat4125 = openmc.Material(material_id=4125, name="")
mat4125.set_density("atom/b-cm", 1.0)  # TODO: Verify density

materials = openmc.Materials([mat1, mat2, mat3, mat10, mat11, mat20, mat500, mat502, mat504, mat506, mat508, mat510, mat513, mat515, mat517, mat519, mat521, mat523, mat525, mat527, mat529, mat531, mat4001, mat4005, mat4009, mat4013, mat4017, mat4021, mat4025, mat4029, mat4033, mat4037, mat4041, mat4045, mat4049, mat4053, mat4057, mat4061, mat4065, mat4069, mat4073, mat4077, mat4081, mat4085, mat4089, mat4093, mat4097, mat4101, mat4105, mat4109, mat4113, mat4117, mat4121, mat4125])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Top of Top Plate
surf2 = openmc.ZPlane(surface_id=2, z0=2.54)

# 126.4412
surf8 = openmc.ZPlane(surface_id=8, z0=105.7910)

# Hole in Top Plate and Side Reflector
surf13 = openmc.ZCylinder(surface_id=13, r=26.797)

# Left Bound of Side Reflector
surf14 = openmc.XPlane(surface_id=14, x0=-44.1452)

# Right Bound of Side Reflector
surf15 = openmc.XPlane(surface_id=15, x0=44.1452)

# Left End Side Reflectors
surf16 = openmc.XPlane(surface_id=16, x0=-27.94)

# Right End Side Relfectors
surf17 = openmc.XPlane(surface_id=17, x0=27.94)

# Front Bound of Side Reflector
surf18 = openmc.YPlane(surface_id=18, y0=-44.1452)

# Back Bound of Side Reflector
surf19 = openmc.YPlane(surface_id=19, y0=44.1452)

# Front End Side Reflectors
surf20 = openmc.YPlane(surface_id=20, y0=-27.94)

# Back End Side Reflectors
surf21 = openmc.YPlane(surface_id=21, y0=27.94)

# Edge of Corner Reflector
surf22 = openmc.XPlane(surface_id=22, x0=0.00)

# Edge of Corner Reflector
surf23 = openmc.YPlane(surface_id=23, y0=0.00)

# Top of Lower Corner Reflector Fourth Row, Bottom Diaphragm
surf27 = openmc.ZPlane(surface_id=27, z0=60.2488)

# Top of Diaphragm
surf28 = openmc.ZPlane(surface_id=28, z0=60.51296)

# Top of Upper Corner Reflector Fifth Row, Bottom Upper Reflector
surf33 = openmc.ZPlane(surface_id=33, z0=82.26806)

# Top of Upper Reflector
surf34 = openmc.ZPlane(surface_id=34, z0=96.69526)

# Core Outside Radius (21" Diameter)
surf36 = openmc.ZCylinder(surface_id=36, r=26.67)

# Top of Iron Plate U1-4
surf40 = openmc.ZPlane(surface_id=40, z0=61.708453)

# Top of HEU Plate U1
surf41 = openmc.ZPlane(surface_id=41, z0=62.008173)

# Top of Iron Plate U1-8
surf45 = openmc.ZPlane(surface_id=45, z0=63.209170)

# Top of Iron Plate U2-4
surf49 = openmc.ZPlane(surface_id=49, z0=64.400007)

# Top of HEU Plate U2
surf50 = openmc.ZPlane(surface_id=50, z0=64.699727)

# Top of Iron Plate U2-8
surf54 = openmc.ZPlane(surface_id=54, z0=65.893103)

# Top of Iron Plate U3-4
surf58 = openmc.ZPlane(surface_id=58, z0=67.083093)

# Top of HEU Plate U3
surf59 = openmc.ZPlane(surface_id=59, z0=67.382813)

# Top of Iron Plate U3-8
surf63 = openmc.ZPlane(surface_id=63, z0=68.575343)

# Top of Iron Plate U4-4
surf67 = openmc.ZPlane(surface_id=67, z0=69.766603)

# Top of HEU Plate U4
surf68 = openmc.ZPlane(surface_id=68, z0=70.066323)

# Top of Iron Plate U4-8
surf72 = openmc.ZPlane(surface_id=72, z0=71.292720)

# Top of Iron Plate U5-4
surf76 = openmc.ZPlane(surface_id=76, z0=72.494563)

# Top of HEU Plate U5
surf77 = openmc.ZPlane(surface_id=77, z0=72.794283)

# Top of Iron Plate U5-8
surf81 = openmc.ZPlane(surface_id=81, z0=74.001207)

# Top of Iron Plate U6-4
surf85 = openmc.ZPlane(surface_id=85, z0=75.189927)

# Top of HEU Plate U6
surf86 = openmc.ZPlane(surface_id=86, z0=75.489647)

# Top of Iron Plate U6-8
surf90 = openmc.ZPlane(surface_id=90, z0=76.677097)

# Top of Iron Plate U7-4
surf94 = openmc.ZPlane(surface_id=94, z0=77.891640)

# Top of HEU Plate U7
surf95 = openmc.ZPlane(surface_id=95, z0=78.191360)

# Top of Iron Plate U7-8
surf99 = openmc.ZPlane(surface_id=99, z0=79.377963)

# Top of Iron Plate U8-4
surf103 = openmc.ZPlane(surface_id=103, z0=80.564990)

# Top of HEU Plate U8
surf104 = openmc.ZPlane(surface_id=104, z0=80.864710)

# Top of Iron Plate U8-8
surf108 = openmc.ZPlane(surface_id=108, z0=82.064437)

# 2.5" Diameter
surf110 = openmc.ZCylinder(surface_id=110, r=3.175)

# 58.754433
surf113 = openmc.ZPlane(surface_id=113, z0=59.054153)

# 58.454713
surf114 = openmc.ZPlane(surface_id=114, z0=58.754433)

# 57.271073
surf118 = openmc.ZPlane(surface_id=118, z0=57.570793)

# 56.084470
surf122 = openmc.ZPlane(surface_id=122, z0=56.384190)

# 55.784750
surf123 = openmc.ZPlane(surface_id=123, z0=56.084470)

# 54.601110
surf127 = openmc.ZPlane(surface_id=127, z0=54.900830)

# 53.420433
surf131 = openmc.ZPlane(surface_id=131, z0=53.720153)

# 53.120713
surf132 = openmc.ZPlane(surface_id=132, z0=53.420433)

# 51.937497
surf136 = openmc.ZPlane(surface_id=136, z0=52.237217)

# 50.758090
surf140 = openmc.ZPlane(surface_id=140, z0=51.057810)

# 50.458370
surf141 = openmc.ZPlane(surface_id=141, z0=50.758090)

# 49.274730
surf145 = openmc.ZPlane(surface_id=145, z0=49.574450)

# 48.092783
surf149 = openmc.ZPlane(surface_id=149, z0=48.392503)

# 47.793063
surf150 = openmc.ZPlane(surface_id=150, z0=48.092783)

# 46.608153
surf154 = openmc.ZPlane(surface_id=154, z0=46.907873)

# 45.423243
surf158 = openmc.ZPlane(surface_id=158, z0=45.722963)

# 45.123523
surf159 = openmc.ZPlane(surface_id=159, z0=45.423243)

# 43.941577
surf163 = openmc.ZPlane(surface_id=163, z0=44.241297)

# 42.758360
surf167 = openmc.ZPlane(surface_id=167, z0=43.058080)

# 42.458640
surf168 = openmc.ZPlane(surface_id=168, z0=42.758360)

# 41.275423
surf172 = openmc.ZPlane(surface_id=172, z0=41.575143)

# 40.092207
surf176 = openmc.ZPlane(surface_id=176, z0=40.391927)

# 39.792487
surf177 = openmc.ZPlane(surface_id=177, z0=40.092207)

# 38.599957
surf181 = openmc.ZPlane(surface_id=181, z0=38.899677)

# 24.172757
surf182 = openmc.ZPlane(surface_id=182, z0=24.472477)

# 20.362757
surf183 = openmc.ZPlane(surface_id=183, z0=20.662477)

# 17.822757
surf184 = openmc.ZPlane(surface_id=184, z0=18.122477)

# 3.75" Diameter Hole in Platen Adapter/Platen
surf185 = openmc.ZCylinder(surface_id=185, r=4.7625)

# Outside of Aligment Tube (2.48" Diameter)
surf186 = openmc.ZCylinder(surface_id=186, r=3.1496)

# Inside of Alignment Tube (2" Diameter)
surf187 = openmc.ZCylinder(surface_id=187, r=2.54)

# Bottom of Alignment Tube
surf188 = openmc.ZPlane(surface_id=188, z0=-3.2512)

# IR Outer Fuel Ring (15" Diameter)
surf410 = openmc.ZCylinder(surface_id=410, r=19.05)

# IR Inner Fuel Ring (6" Diameter)
surf860 = openmc.ZCylinder(surface_id=860, r=7.62)

# 59.94908
surf1080 = openmc.ZPlane(surface_id=1080, z0=59.054153)


# Cell: VOID
cell0 = openmc.Cell(cell_id=0, fill=None, name="VOID")  # mat0 undefined, using void
cell0.region = +surf103 & -surf104 & -surf860

# Cell: MAT1
cell1 = openmc.Cell(cell_id=1, fill=mat1, name="MAT1")
cell1.region = -surf36 & -surf182 & +surf183 & +surf185

# Cell: MAT1
cell2 = openmc.Cell(cell_id=2, fill=mat1, name="MAT1")
cell2.region = -surf36 & -surf183 & +surf184 & +surf185

# Cell: MAT1
cell3 = openmc.Cell(cell_id=3, fill=mat1, name="MAT1")
cell3.region = -surf186 & +surf187 & -surf27 & +surf188

# Cell: MAT2
cell4 = openmc.Cell(cell_id=4, fill=mat2, name="MAT2")
cell4.region = +surf14 & -surf16 & +surf18 & -surf21 & +surf2 & -surf8

# Cell: MAT2
cell5 = openmc.Cell(cell_id=5, fill=mat2, name="MAT2")
cell5.region = +surf16 & -surf15 & +surf18 & -surf20 & +surf2 & -surf8

# Cell: MAT2
cell6 = openmc.Cell(cell_id=6, fill=mat2, name="MAT2")
cell6.region = +surf17 & -surf15 & +surf20 & -surf19 & +surf2 & -surf8

# Cell: MAT2
cell7 = openmc.Cell(cell_id=7, fill=mat2, name="MAT2")
cell7.region = +surf14 & -surf17 & +surf21 & -surf19 & +surf2 & -surf8

# Cell: MAT3
cell8 = openmc.Cell(cell_id=8, fill=mat3, name="MAT3")
cell8.region = +surf16 & -surf17 & +surf20 & -surf21 & +surf27 & -surf28

# Cell: MAT10
cell9 = openmc.Cell(cell_id=9, fill=mat10, name="MAT10")
cell9.region = -surf113 & +surf114 & -surf860 & +surf110

# Cell: MAT20
cell10 = openmc.Cell(cell_id=10, fill=mat20, name="MAT20")
cell10.region = +surf16 & -surf22 & +surf20 & -surf23 & +surf2 & -surf27 & +surf13

# Cell: MAT20
cell11 = openmc.Cell(cell_id=11, fill=mat20, name="MAT20")
cell11.region = +surf22 & -surf17 & +surf20 & -surf23 & +surf2 & -surf27 & +surf13

# Cell: MAT20
cell12 = openmc.Cell(cell_id=12, fill=mat20, name="MAT20")
cell12.region = +surf22 & -surf17 & +surf23 & -surf21 & +surf2 & -surf27 & +surf13

# Cell: MAT20
cell13 = openmc.Cell(cell_id=13, fill=mat20, name="MAT20")
cell13.region = +surf16 & -surf22 & +surf23 & -surf21 & +surf2 & -surf27 & +surf13

# Cell: MAT20
cell14 = openmc.Cell(cell_id=14, fill=mat20, name="MAT20")
cell14.region = +surf16 & -surf22 & +surf20 & -surf23 & +surf28 & -surf33 & +surf13

# Cell: MAT20
cell15 = openmc.Cell(cell_id=15, fill=mat20, name="MAT20")
cell15.region = +surf22 & -surf17 & +surf20 & -surf23 & +surf28 & -surf33 & +surf13

# Cell: MAT20
cell16 = openmc.Cell(cell_id=16, fill=mat20, name="MAT20")
cell16.region = +surf22 & -surf17 & +surf23 & -surf21 & +surf28 & -surf33 & +surf13

# Cell: MAT20
cell17 = openmc.Cell(cell_id=17, fill=mat20, name="MAT20")
cell17.region = +surf16 & -surf22 & +surf23 & -surf21 & +surf28 & -surf33 & +surf13

# Cell: MAT20
cell18 = openmc.Cell(cell_id=18, fill=mat20, name="MAT20")
cell18.region = +surf16 & -surf17 & +surf20 & -surf21 & +surf33 & -surf34

# Cell: MAT20
cell19 = openmc.Cell(cell_id=19, fill=mat20, name="MAT20")
cell19.region = -surf36 & -surf181 & +surf182 & +surf110

# Cell: MAT500
cell20 = openmc.Cell(cell_id=20, fill=mat500, name="MAT500")
cell20.region = -surf36 & +surf40 & -surf41

# Cell: MAT502
cell21 = openmc.Cell(cell_id=21, fill=mat502, name="MAT502")
cell21.region = -surf36 & +surf49 & -surf50

# Cell: MAT504
cell22 = openmc.Cell(cell_id=22, fill=mat504, name="MAT504")
cell22.region = -surf36 & +surf58 & -surf59

# Cell: MAT506
cell23 = openmc.Cell(cell_id=23, fill=mat506, name="MAT506")
cell23.region = -surf36 & +surf67 & -surf68

# Cell: MAT508
cell24 = openmc.Cell(cell_id=24, fill=mat508, name="MAT508")
cell24.region = -surf36 & +surf76 & -surf77

# Cell: MAT510
cell25 = openmc.Cell(cell_id=25, fill=mat510, name="MAT510")
cell25.region = -surf36 & +surf85 & -surf86

# Cell: MAT513
cell26 = openmc.Cell(cell_id=26, fill=mat513, name="MAT513")
cell26.region = -surf36 & +surf94 & -surf95

# Cell: MAT515
cell27 = openmc.Cell(cell_id=27, fill=mat515, name="MAT515")
cell27.region = -surf36 & +surf103 & -surf104 & +surf860

# Cell: MAT517
cell28 = openmc.Cell(cell_id=28, fill=mat517, name="MAT517")
cell28.region = -surf36 & -surf113 & +surf114 & +surf860

# Cell: MAT519
cell29 = openmc.Cell(cell_id=29, fill=mat519, name="MAT519")
cell29.region = -surf36 & -surf122 & +surf123 & +surf110

# Cell: MAT521
cell30 = openmc.Cell(cell_id=30, fill=mat521, name="MAT521")
cell30.region = -surf36 & -surf131 & +surf132 & +surf110

# Cell: MAT523
cell31 = openmc.Cell(cell_id=31, fill=mat523, name="MAT523")
cell31.region = -surf36 & -surf140 & +surf141 & +surf110

# Cell: MAT525
cell32 = openmc.Cell(cell_id=32, fill=mat525, name="MAT525")
cell32.region = -surf36 & -surf149 & +surf150 & +surf110

# Cell: MAT527
cell33 = openmc.Cell(cell_id=33, fill=mat527, name="MAT527")
cell33.region = -surf36 & -surf158 & +surf159 & +surf110

# Cell: MAT529
cell34 = openmc.Cell(cell_id=34, fill=mat529, name="MAT529")
cell34.region = -surf36 & -surf167 & +surf168 & +surf110

# Cell: MAT531
cell35 = openmc.Cell(cell_id=35, fill=mat531, name="MAT531")
cell35.region = -surf36 & -surf176 & +surf177 & +surf110

# Cell: MAT4001
cell36 = openmc.Cell(cell_id=36, fill=mat4001, name="MAT4001")
cell36.region = -surf36 & +surf28 & -surf40

# Cell: MAT4005
cell37 = openmc.Cell(cell_id=37, fill=mat4005, name="MAT4005")
cell37.region = -surf36 & +surf41 & -surf45

# Cell: MAT4009
cell38 = openmc.Cell(cell_id=38, fill=mat4009, name="MAT4009")
cell38.region = -surf36 & +surf45 & -surf49

# Cell: MAT4013
cell39 = openmc.Cell(cell_id=39, fill=mat4013, name="MAT4013")
cell39.region = -surf36 & +surf50 & -surf54

# Cell: MAT4017
cell40 = openmc.Cell(cell_id=40, fill=mat4017, name="MAT4017")
cell40.region = -surf36 & +surf54 & -surf58

# Cell: MAT4021
cell41 = openmc.Cell(cell_id=41, fill=mat4021, name="MAT4021")
cell41.region = -surf36 & +surf59 & -surf63

# Cell: MAT4025
cell42 = openmc.Cell(cell_id=42, fill=mat4025, name="MAT4025")
cell42.region = -surf36 & +surf63 & -surf67

# Cell: MAT4029
cell43 = openmc.Cell(cell_id=43, fill=mat4029, name="MAT4029")
cell43.region = -surf36 & +surf68 & -surf72

# Cell: MAT4033
cell44 = openmc.Cell(cell_id=44, fill=mat4033, name="MAT4033")
cell44.region = -surf36 & +surf72 & -surf76

# Cell: MAT4037
cell45 = openmc.Cell(cell_id=45, fill=mat4037, name="MAT4037")
cell45.region = -surf36 & +surf77 & -surf81

# Cell: MAT4041
cell46 = openmc.Cell(cell_id=46, fill=mat4041, name="MAT4041")
cell46.region = -surf36 & +surf81 & -surf85

# Cell: MAT4045
cell47 = openmc.Cell(cell_id=47, fill=mat4045, name="MAT4045")
cell47.region = -surf36 & +surf86 & -surf90

# Cell: MAT4049
cell48 = openmc.Cell(cell_id=48, fill=mat4049, name="MAT4049")
cell48.region = -surf36 & +surf90 & -surf94

# Cell: MAT4053
cell49 = openmc.Cell(cell_id=49, fill=mat4053, name="MAT4053")
cell49.region = -surf36 & +surf95 & -surf99

# Cell: MAT4057
cell50 = openmc.Cell(cell_id=50, fill=mat4057, name="MAT4057")
cell50.region = -surf36 & +surf99 & -surf103

# Cell: MAT4061
cell51 = openmc.Cell(cell_id=51, fill=mat4061, name="MAT4061")
cell51.region = -surf36 & +surf104 & -surf108

# Cell: MAT4065
cell52 = openmc.Cell(cell_id=52, fill=mat4065, name="MAT4065")
cell52.region = -surf36 & -surf1080 & +surf113 & +surf110

# Cell: MAT4069
cell53 = openmc.Cell(cell_id=53, fill=mat4069, name="MAT4069")
cell53.region = -surf36 & -surf114 & +surf118 & +surf110

# Cell: MAT4073
cell54 = openmc.Cell(cell_id=54, fill=mat4073, name="MAT4073")
cell54.region = -surf36 & -surf118 & +surf122 & +surf110

# Cell: MAT4077
cell55 = openmc.Cell(cell_id=55, fill=mat4077, name="MAT4077")
cell55.region = -surf36 & -surf123 & +surf127 & +surf110

# Cell: MAT4081
cell56 = openmc.Cell(cell_id=56, fill=mat4081, name="MAT4081")
cell56.region = -surf36 & -surf127 & +surf131 & +surf110

# Cell: MAT4085
cell57 = openmc.Cell(cell_id=57, fill=mat4085, name="MAT4085")
cell57.region = -surf36 & -surf132 & +surf136 & +surf110

# Cell: MAT4089
cell58 = openmc.Cell(cell_id=58, fill=mat4089, name="MAT4089")
cell58.region = -surf36 & -surf136 & +surf140 & +surf110

# Cell: MAT4093
cell59 = openmc.Cell(cell_id=59, fill=mat4093, name="MAT4093")
cell59.region = -surf36 & -surf141 & +surf145 & +surf110

# Cell: MAT4097
cell60 = openmc.Cell(cell_id=60, fill=mat4097, name="MAT4097")
cell60.region = -surf36 & -surf145 & +surf149 & +surf110

# Cell: MAT4101
cell61 = openmc.Cell(cell_id=61, fill=mat4101, name="MAT4101")
cell61.region = -surf36 & -surf150 & +surf154 & +surf110

# Cell: MAT4105
cell62 = openmc.Cell(cell_id=62, fill=mat4105, name="MAT4105")
cell62.region = -surf36 & -surf154 & +surf158 & +surf110

# Cell: MAT4109
cell63 = openmc.Cell(cell_id=63, fill=mat4109, name="MAT4109")
cell63.region = -surf36 & -surf159 & +surf163 & +surf110

# Cell: MAT4113
cell64 = openmc.Cell(cell_id=64, fill=mat4113, name="MAT4113")
cell64.region = -surf36 & -surf163 & +surf167 & +surf110

# Cell: MAT4117
cell65 = openmc.Cell(cell_id=65, fill=mat4117, name="MAT4117")
cell65.region = -surf36 & -surf168 & +surf172 & +surf110

# Cell: MAT4121
cell66 = openmc.Cell(cell_id=66, fill=mat4121, name="MAT4121")
cell66.region = -surf36 & -surf172 & +surf176 & +surf110

# Cell: MAT4125
cell67 = openmc.Cell(cell_id=67, fill=mat4125, name="MAT4125")
cell67.region = -surf36 & -surf177 & +surf181 & +surf110

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary (6 planes)
# TODO: Adjust dimensions to encompass your entire geometry
boundary_xmin = openmc.XPlane(surface_id=10000, x0=-200, boundary_type="vacuum")
boundary_xmax = openmc.XPlane(surface_id=10001, x0=200, boundary_type="vacuum")
boundary_ymin = openmc.YPlane(surface_id=10002, y0=-200, boundary_type="vacuum")
boundary_ymax = openmc.YPlane(surface_id=10003, y0=200, boundary_type="vacuum")
boundary_zmin = openmc.ZPlane(surface_id=10004, z0=-200, boundary_type="vacuum")
boundary_zmax = openmc.ZPlane(surface_id=10005, z0=200, boundary_type="vacuum")

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
outer_region = outer_region & ~cell25.region
outer_region = outer_region & ~cell26.region
outer_region = outer_region & ~cell27.region
outer_region = outer_region & ~cell28.region
outer_region = outer_region & ~cell29.region
outer_region = outer_region & ~cell30.region
outer_region = outer_region & ~cell31.region
outer_region = outer_region & ~cell32.region
outer_region = outer_region & ~cell33.region
outer_region = outer_region & ~cell34.region
outer_region = outer_region & ~cell35.region
outer_region = outer_region & ~cell36.region
outer_region = outer_region & ~cell37.region
outer_region = outer_region & ~cell38.region
outer_region = outer_region & ~cell39.region
outer_region = outer_region & ~cell40.region
outer_region = outer_region & ~cell41.region
outer_region = outer_region & ~cell42.region
outer_region = outer_region & ~cell43.region
outer_region = outer_region & ~cell44.region
outer_region = outer_region & ~cell45.region
outer_region = outer_region & ~cell46.region
outer_region = outer_region & ~cell47.region
outer_region = outer_region & ~cell48.region
outer_region = outer_region & ~cell49.region
outer_region = outer_region & ~cell50.region
outer_region = outer_region & ~cell51.region
outer_region = outer_region & ~cell52.region
outer_region = outer_region & ~cell53.region
outer_region = outer_region & ~cell54.region
outer_region = outer_region & ~cell55.region
outer_region = outer_region & ~cell56.region
outer_region = outer_region & ~cell57.region
outer_region = outer_region & ~cell58.region
outer_region = outer_region & ~cell59.region
outer_region = outer_region & ~cell60.region
outer_region = outer_region & ~cell61.region
outer_region = outer_region & ~cell62.region
outer_region = outer_region & ~cell63.region
outer_region = outer_region & ~cell64.region
outer_region = outer_region & ~cell65.region
outer_region = outer_region & ~cell66.region
outer_region = outer_region & ~cell67.region
outer_cell = openmc.Cell(cell_id=68, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24, cell25, cell26, cell27, cell28, cell29, cell30, cell31, cell32, cell33, cell34, cell35, cell36, cell37, cell38, cell39, cell40, cell41, cell42, cell43, cell44, cell45, cell46, cell47, cell48, cell49, cell50, cell51, cell52, cell53, cell54, cell55, cell56, cell57, cell58, cell59, cell60, cell61, cell62, cell63, cell64, cell65, cell66, cell67, outer_cell])
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
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
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
