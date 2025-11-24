"""
HMF072-2: Zeus Iron Core Benchmark (No Al Shim)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")

mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")

mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")

mat10 = openmc.Material(material_id=10)
mat10.set_density("sum")

mat11 = openmc.Material(material_id=11)
mat11.set_density("sum")

mat20 = openmc.Material(material_id=20)
mat20.set_density("sum")

mat500 = openmc.Material(material_id=500)
mat500.set_density("sum")

mat502 = openmc.Material(material_id=502)
mat502.set_density("sum")

mat504 = openmc.Material(material_id=504)
mat504.set_density("sum")

mat506 = openmc.Material(material_id=506)
mat506.set_density("sum")

mat508 = openmc.Material(material_id=508)
mat508.set_density("sum")

mat510 = openmc.Material(material_id=510)
mat510.set_density("sum")

mat513 = openmc.Material(material_id=513)
mat513.set_density("sum")

mat515 = openmc.Material(material_id=515)
mat515.set_density("sum")

mat517 = openmc.Material(material_id=517)
mat517.set_density("sum")

mat519 = openmc.Material(material_id=519)
mat519.set_density("sum")

mat521 = openmc.Material(material_id=521)
mat521.set_density("sum")

mat523 = openmc.Material(material_id=523)
mat523.set_density("sum")

mat525 = openmc.Material(material_id=525)
mat525.set_density("sum")

mat527 = openmc.Material(material_id=527)
mat527.set_density("sum")

mat529 = openmc.Material(material_id=529)
mat529.set_density("sum")

mat531 = openmc.Material(material_id=531)
mat531.set_density("sum")

mat4001 = openmc.Material(material_id=4001)
mat4001.set_density("sum")

mat4005 = openmc.Material(material_id=4005)
mat4005.set_density("sum")

mat4009 = openmc.Material(material_id=4009)
mat4009.set_density("sum")

mat4013 = openmc.Material(material_id=4013)
mat4013.set_density("sum")

mat4017 = openmc.Material(material_id=4017)
mat4017.set_density("sum")

mat4021 = openmc.Material(material_id=4021)
mat4021.set_density("sum")

mat4025 = openmc.Material(material_id=4025)
mat4025.set_density("sum")

mat4029 = openmc.Material(material_id=4029)
mat4029.set_density("sum")

mat4033 = openmc.Material(material_id=4033)
mat4033.set_density("sum")

mat4037 = openmc.Material(material_id=4037)
mat4037.set_density("sum")

mat4041 = openmc.Material(material_id=4041)
mat4041.set_density("sum")

mat4045 = openmc.Material(material_id=4045)
mat4045.set_density("sum")

mat4049 = openmc.Material(material_id=4049)
mat4049.set_density("sum")

mat4053 = openmc.Material(material_id=4053)
mat4053.set_density("sum")

mat4057 = openmc.Material(material_id=4057)
mat4057.set_density("sum")

mat4061 = openmc.Material(material_id=4061)
mat4061.set_density("sum")

mat4065 = openmc.Material(material_id=4065)
mat4065.set_density("sum")

mat4069 = openmc.Material(material_id=4069)
mat4069.set_density("sum")

mat4073 = openmc.Material(material_id=4073)
mat4073.set_density("sum")

mat4077 = openmc.Material(material_id=4077)
mat4077.set_density("sum")

mat4081 = openmc.Material(material_id=4081)
mat4081.set_density("sum")

mat4085 = openmc.Material(material_id=4085)
mat4085.set_density("sum")

mat4089 = openmc.Material(material_id=4089)
mat4089.set_density("sum")

mat4093 = openmc.Material(material_id=4093)
mat4093.set_density("sum")

mat4097 = openmc.Material(material_id=4097)
mat4097.set_density("sum")

mat4101 = openmc.Material(material_id=4101)
mat4101.set_density("sum")

mat4105 = openmc.Material(material_id=4105)
mat4105.set_density("sum")

mat4109 = openmc.Material(material_id=4109)
mat4109.set_density("sum")

mat4113 = openmc.Material(material_id=4113)
mat4113.set_density("sum")

mat4117 = openmc.Material(material_id=4117)
mat4117.set_density("sum")

mat4121 = openmc.Material(material_id=4121)
mat4121.set_density("sum")

mat4125 = openmc.Material(material_id=4125)
mat4125.set_density("sum")

materials = openmc.Materials([mat1, mat2, mat3, mat10, mat11, mat20, mat500, mat502, mat504, mat506, mat508, mat510, mat513, mat515, mat517, mat519, mat521, mat523, mat525, mat527, mat529, mat531, mat4001, mat4005, mat4009, mat4013, mat4017, mat4021, mat4025, mat4029, mat4033, mat4037, mat4041, mat4045, mat4049, mat4053, mat4057, mat4061, mat4065, mat4069, mat4073, mat4077, mat4081, mat4085, mat4089, mat4093, mat4097, mat4101, mat4105, mat4109, mat4113, mat4117, mat4121, mat4125])

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

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# MAT2
cell1 = openmc.Cell(cell_id=1, fill=mat2)
cell1.region = +surf14 & -surf16 & +surf18 & -surf21 & +surf2 & -surf8

# MAT2
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = +surf16 & -surf15 & +surf18 & -surf20 & +surf2 & -surf8

# MAT2
cell3 = openmc.Cell(cell_id=3, fill=mat2)
cell3.region = +surf17 & -surf15 & +surf20 & -surf19 & +surf2 & -surf8

# MAT2
cell4 = openmc.Cell(cell_id=4, fill=mat2)
cell4.region = +surf14 & -surf17 & +surf21 & -surf19 & +surf2 & -surf8

# MAT20
cell5 = openmc.Cell(cell_id=5, fill=mat20)
cell5.region = +surf16 & -surf22 & +surf20 & -surf23 & +surf2 & -surf27 & +surf13

# MAT20
cell6 = openmc.Cell(cell_id=6, fill=mat20)
cell6.region = +surf22 & -surf17 & +surf20 & -surf23 & +surf2 & -surf27 & +surf13

# MAT20
cell7 = openmc.Cell(cell_id=7, fill=mat20)
cell7.region = +surf22 & -surf17 & +surf23 & -surf21 & +surf2 & -surf27 & +surf13

# MAT20
cell8 = openmc.Cell(cell_id=8, fill=mat20)
cell8.region = +surf16 & -surf22 & +surf23 & -surf21 & +surf2 & -surf27 & +surf13

# MAT3
cell9 = openmc.Cell(cell_id=9, fill=mat3)
cell9.region = +surf16 & -surf17 & +surf20 & -surf21 & +surf27 & -surf28

# MAT20
cell10 = openmc.Cell(cell_id=10, fill=mat20)
cell10.region = +surf16 & -surf22 & +surf20 & -surf23 & +surf28 & -surf33 & +surf13

# MAT20
cell11 = openmc.Cell(cell_id=11, fill=mat20)
cell11.region = +surf22 & -surf17 & +surf20 & -surf23 & +surf28 & -surf33 & +surf13

# MAT20
cell12 = openmc.Cell(cell_id=12, fill=mat20)
cell12.region = +surf22 & -surf17 & +surf23 & -surf21 & +surf28 & -surf33 & +surf13

# MAT20
cell13 = openmc.Cell(cell_id=13, fill=mat20)
cell13.region = +surf16 & -surf22 & +surf23 & -surf21 & +surf28 & -surf33 & +surf13

# MAT20
cell14 = openmc.Cell(cell_id=14, fill=mat20)
cell14.region = +surf16 & -surf17 & +surf20 & -surf21 & +surf33 & -surf34

# MAT4001
cell15 = openmc.Cell(cell_id=15, fill=mat4001)
cell15.region = -surf36 & +surf28 & -surf40

# MAT500
cell16 = openmc.Cell(cell_id=16, fill=mat500)
cell16.region = -surf36 & +surf40 & -surf41

# MAT4005
cell17 = openmc.Cell(cell_id=17, fill=mat4005)
cell17.region = -surf36 & +surf41 & -surf45

# MAT4009
cell18 = openmc.Cell(cell_id=18, fill=mat4009)
cell18.region = -surf36 & +surf45 & -surf49

# MAT502
cell19 = openmc.Cell(cell_id=19, fill=mat502)
cell19.region = -surf36 & +surf49 & -surf50

# MAT4013
cell20 = openmc.Cell(cell_id=20, fill=mat4013)
cell20.region = -surf36 & +surf50 & -surf54

# MAT4017
cell21 = openmc.Cell(cell_id=21, fill=mat4017)
cell21.region = -surf36 & +surf54 & -surf58

# MAT504
cell22 = openmc.Cell(cell_id=22, fill=mat504)
cell22.region = -surf36 & +surf58 & -surf59

# MAT4021
cell23 = openmc.Cell(cell_id=23, fill=mat4021)
cell23.region = -surf36 & +surf59 & -surf63

# MAT4025
cell24 = openmc.Cell(cell_id=24, fill=mat4025)
cell24.region = -surf36 & +surf63 & -surf67

# MAT506
cell25 = openmc.Cell(cell_id=25, fill=mat506)
cell25.region = -surf36 & +surf67 & -surf68

# MAT4029
cell26 = openmc.Cell(cell_id=26, fill=mat4029)
cell26.region = -surf36 & +surf68 & -surf72

# MAT4033
cell27 = openmc.Cell(cell_id=27, fill=mat4033)
cell27.region = -surf36 & +surf72 & -surf76

# MAT508
cell28 = openmc.Cell(cell_id=28, fill=mat508)
cell28.region = -surf36 & +surf76 & -surf77

# MAT4037
cell29 = openmc.Cell(cell_id=29, fill=mat4037)
cell29.region = -surf36 & +surf77 & -surf81

# MAT4041
cell30 = openmc.Cell(cell_id=30, fill=mat4041)
cell30.region = -surf36 & +surf81 & -surf85

# MAT510
cell31 = openmc.Cell(cell_id=31, fill=mat510)
cell31.region = -surf36 & +surf85 & -surf86

# MAT4045
cell32 = openmc.Cell(cell_id=32, fill=mat4045)
cell32.region = -surf36 & +surf86 & -surf90

# MAT4049
cell33 = openmc.Cell(cell_id=33, fill=mat4049)
cell33.region = -surf36 & +surf90 & -surf94

# MAT513
cell34 = openmc.Cell(cell_id=34, fill=mat513)
cell34.region = -surf36 & +surf94 & -surf95

# MAT4053
cell35 = openmc.Cell(cell_id=35, fill=mat4053)
cell35.region = -surf36 & +surf95 & -surf99

# MAT4057
cell36 = openmc.Cell(cell_id=36, fill=mat4057)
cell36.region = -surf36 & +surf99 & -surf103

# MAT515
cell37 = openmc.Cell(cell_id=37, fill=mat515)
cell37.region = -surf36 & +surf103 & -surf104 & +surf860

# VOID
cell38 = openmc.Cell(cell_id=38)
cell38.region = +surf103 & -surf104 & -surf860

# MAT4061
cell39 = openmc.Cell(cell_id=39, fill=mat4061)
cell39.region = -surf36 & +surf104 & -surf108

# MAT4065
cell40 = openmc.Cell(cell_id=40, fill=mat4065)
cell40.region = -surf36 & -surf1080 & +surf113 & +surf110

# MAT517
cell41 = openmc.Cell(cell_id=41, fill=mat517)
cell41.region = -surf36 & -surf113 & +surf114 & +surf860

# MAT10
cell42 = openmc.Cell(cell_id=42, fill=mat10)
cell42.region = -surf113 & +surf114 & -surf860 & +surf110

# MAT4069
cell43 = openmc.Cell(cell_id=43, fill=mat4069)
cell43.region = -surf36 & -surf114 & +surf118 & +surf110

# MAT4073
cell44 = openmc.Cell(cell_id=44, fill=mat4073)
cell44.region = -surf36 & -surf118 & +surf122 & +surf110

# MAT519
cell45 = openmc.Cell(cell_id=45, fill=mat519)
cell45.region = -surf36 & -surf122 & +surf123 & +surf110

# MAT4077
cell46 = openmc.Cell(cell_id=46, fill=mat4077)
cell46.region = -surf36 & -surf123 & +surf127 & +surf110

# MAT4081
cell47 = openmc.Cell(cell_id=47, fill=mat4081)
cell47.region = -surf36 & -surf127 & +surf131 & +surf110

# MAT521
cell48 = openmc.Cell(cell_id=48, fill=mat521)
cell48.region = -surf36 & -surf131 & +surf132 & +surf110

# MAT4085
cell49 = openmc.Cell(cell_id=49, fill=mat4085)
cell49.region = -surf36 & -surf132 & +surf136 & +surf110

# MAT4089
cell50 = openmc.Cell(cell_id=50, fill=mat4089)
cell50.region = -surf36 & -surf136 & +surf140 & +surf110

# MAT523
cell51 = openmc.Cell(cell_id=51, fill=mat523)
cell51.region = -surf36 & -surf140 & +surf141 & +surf110

# MAT4093
cell52 = openmc.Cell(cell_id=52, fill=mat4093)
cell52.region = -surf36 & -surf141 & +surf145 & +surf110

# MAT4097
cell53 = openmc.Cell(cell_id=53, fill=mat4097)
cell53.region = -surf36 & -surf145 & +surf149 & +surf110

# MAT525
cell54 = openmc.Cell(cell_id=54, fill=mat525)
cell54.region = -surf36 & -surf149 & +surf150 & +surf110

# MAT4101
cell55 = openmc.Cell(cell_id=55, fill=mat4101)
cell55.region = -surf36 & -surf150 & +surf154 & +surf110

# MAT4105
cell56 = openmc.Cell(cell_id=56, fill=mat4105)
cell56.region = -surf36 & -surf154 & +surf158 & +surf110

# MAT527
cell57 = openmc.Cell(cell_id=57, fill=mat527)
cell57.region = -surf36 & -surf158 & +surf159 & +surf110

# MAT4109
cell58 = openmc.Cell(cell_id=58, fill=mat4109)
cell58.region = -surf36 & -surf159 & +surf163 & +surf110

# MAT4113
cell59 = openmc.Cell(cell_id=59, fill=mat4113)
cell59.region = -surf36 & -surf163 & +surf167 & +surf110

# MAT529
cell60 = openmc.Cell(cell_id=60, fill=mat529)
cell60.region = -surf36 & -surf167 & +surf168 & +surf110

# MAT4117
cell61 = openmc.Cell(cell_id=61, fill=mat4117)
cell61.region = -surf36 & -surf168 & +surf172 & +surf110

# MAT4121
cell62 = openmc.Cell(cell_id=62, fill=mat4121)
cell62.region = -surf36 & -surf172 & +surf176 & +surf110

# MAT531
cell63 = openmc.Cell(cell_id=63, fill=mat531)
cell63.region = -surf36 & -surf176 & +surf177 & +surf110

# MAT4125
cell64 = openmc.Cell(cell_id=64, fill=mat4125)
cell64.region = -surf36 & -surf177 & +surf181 & +surf110

# MAT20
cell65 = openmc.Cell(cell_id=65, fill=mat20)
cell65.region = -surf36 & -surf181 & +surf182 & +surf110

# MAT1
cell66 = openmc.Cell(cell_id=66, fill=mat1)
cell66.region = -surf36 & -surf182 & +surf183 & +surf185

# MAT1
cell67 = openmc.Cell(cell_id=67, fill=mat1)
cell67.region = -surf36 & -surf183 & +surf184 & +surf185

# MAT1
cell68 = openmc.Cell(cell_id=68, fill=mat1)
cell68.region = -surf186 & +surf187 & -surf27 & +surf188

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24, cell25, cell26, cell27, cell28, cell29, cell30, cell31, cell32, cell33, cell34, cell35, cell36, cell37, cell38, cell39, cell40, cell41, cell42, cell43, cell44, cell45, cell46, cell47, cell48, cell49, cell50, cell51, cell52, cell53, cell54, cell55, cell56, cell57, cell58, cell59, cell60, cell61, cell62, cell63, cell64, cell65, cell66, cell67, cell68])
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
source.space = openmc.stats.Box((-14.0, -12.26, 39.22), (14.0, 12.26, 81.75))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
