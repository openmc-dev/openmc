"""
IEU-MET-FAST-007 (Rev. 2): Big Ten Detailed Benchmark Model (translated)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Highly Enriched Uranium (93 wt.%)
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 4.981400e-04)
mat1.add_nuclide("U235", 4.503400e-02)
mat1.add_nuclide("U236", 1.323600e-04)
mat1.add_nuclide("U238", 2.605600e-03)

# Intermediate Enriched Uranium (10 wt.%)
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("U234", 2.476100e-05)
mat2.add_nuclide("U235", 4.846100e-03)
mat2.add_nuclide("U236", 4.334800e-05)
mat2.add_nuclide("U238", 4.269500e-02)

# Natural Uranium
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("U234", 2.651800e-06)
mat3.add_nuclide("U235", 3.470100e-04)
mat3.add_nuclide("U238", 4.784600e-02)

# Depleted Uranium
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_nuclide("U234", 2.867200e-07)
mat4.add_nuclide("U235", 1.005800e-04)
mat4.add_nuclide("U236", 1.146800e-06)
mat4.add_nuclide("U238", 4.767700e-02)

# Natural Uranium (P66)
mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")
mat5.add_nuclide("U234", 2.645800e-06)
mat5.add_nuclide("U235", 3.462300e-04)
mat5.add_nuclide("U238", 4.773800e-02)

# Steel 347
mat6 = openmc.Material(material_id=6)
mat6.set_density("sum")
mat6.add_element("Fe", 5.779800e-02)
mat6.add_element("Cr", 1.667800e-02)
mat6.add_element("Ni", 9.029600e-03)
mat6.add_element("Mn", 1.753900e-03)
mat6.add_element("Si", 1.715400e-03)
mat6.add_element("Nb", 5.185500e-04)

# Steel 304
mat7 = openmc.Material(material_id=7)
mat7.set_density("sum")
mat7.add_element("Fe", 5.930800e-02)
mat7.add_element("Cr", 1.760400e-02)
mat7.add_element("Ni", 7.593000e-03)
mat7.add_element("Mn", 1.753900e-03)
mat7.add_element("Si", 1.715400e-03)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7])

# ==============================================================================
# Geometry
# ==============================================================================


# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# NatU
cell1 = openmc.Cell(cell_id=1, fill=mat3)
cell1.region = +surf103 & -surf108 & +surf4 & -surf3

# HEU
cell2 = openmc.Cell(cell_id=2, fill=mat1)
cell2.region = +surf103 & -surf108 & +surf5 & -surf4

# NatU
cell3 = openmc.Cell(cell_id=3, fill=mat3)
cell3.region = +surf103 & -surf108 & +surf6 & -surf5

# HEU
cell4 = openmc.Cell(cell_id=4, fill=mat1)
cell4.region = +surf103 & -surf108 & +surf7 & -surf6

# NatU
cell5 = openmc.Cell(cell_id=5, fill=mat3)
cell5.region = +surf103 & -surf108 & +surf8 & -surf7

# HEU
cell6 = openmc.Cell(cell_id=6, fill=mat1)
cell6.region = +surf103 & -surf108 & +surf9 & -surf8

# NatU
cell7 = openmc.Cell(cell_id=7, fill=mat3)
cell7.region = +surf103 & -surf108 & +surf10 & -surf9

# HEU
cell8 = openmc.Cell(cell_id=8, fill=mat1)
cell8.region = +surf103 & -surf108 & +surf11 & -surf10

# NatU
cell9 = openmc.Cell(cell_id=9, fill=mat3)
cell9.region = +surf103 & -surf108 & +surf12 & -surf11

# HEU
cell10 = openmc.Cell(cell_id=10, fill=mat1)
cell10.region = +surf103 & -surf108 & +surf13 & -surf12

# NatU
cell11 = openmc.Cell(cell_id=11, fill=mat3)
cell11.region = +surf107 & -surf108 & +surf15 & -surf13

# HEU
cell12 = openmc.Cell(cell_id=12, fill=mat1)
cell12.region = +surf107 & -surf108 & +surf16 & -surf15

# D38
cell13 = openmc.Cell(cell_id=13, fill=mat4)
cell13.region = +surf106 & -surf108 & +surf18 & -surf16

# HEU
cell14 = openmc.Cell(cell_id=14, fill=mat1)
cell14.region = +surf107 & -surf108 & +surf20 & -surf19

# NatU
cell15 = openmc.Cell(cell_id=15, fill=mat3)
cell15.region = +surf107 & -surf108 & +surf21 & -surf20

# NatU
cell16 = openmc.Cell(cell_id=16, fill=mat3)
cell16.region = +surf107 & -surf108 & +surf23 & -surf22

# HEU
cell17 = openmc.Cell(cell_id=17, fill=mat1)
cell17.region = +surf107 & -surf108 & +surf24 & -surf23

# NatU
cell18 = openmc.Cell(cell_id=18, fill=mat3)
cell18.region = +surf107 & -surf108 & +surf25 & -surf24

# NatU
cell19 = openmc.Cell(cell_id=19, fill=mat3)
cell19.region = +surf107 & -surf108 & +surf27 & -surf26

# HEU
cell20 = openmc.Cell(cell_id=20, fill=mat1)
cell20.region = +surf107 & -surf108 & +surf28 & -surf27

# NatU
cell21 = openmc.Cell(cell_id=21, fill=mat3)
cell21.region = +surf107 & -surf108 & +surf29 & -surf28

# NatU
cell22 = openmc.Cell(cell_id=22, fill=mat3)
cell22.region = +surf107 & -surf108 & +surf31 & -surf30

# HEU
cell23 = openmc.Cell(cell_id=23, fill=mat1)
cell23.region = +surf107 & -surf108 & +surf32 & -surf31

# NatU
cell24 = openmc.Cell(cell_id=24, fill=mat3)
cell24.region = +surf107 & -surf108 & +surf33 & -surf32

# NatU
cell25 = openmc.Cell(cell_id=25, fill=mat3)
cell25.region = +surf107 & -surf108 & +surf35 & -surf34

# HEU
cell26 = openmc.Cell(cell_id=26, fill=mat1)
cell26.region = +surf107 & -surf108 & +surf36 & -surf35

# NatU
cell27 = openmc.Cell(cell_id=27, fill=mat3)
cell27.region = +surf107 & -surf108 & +surf37 & -surf36

# NatU
cell28 = openmc.Cell(cell_id=28, fill=mat3)
cell28.region = +surf107 & -surf108 & +surf39 & -surf38

# HEU
cell29 = openmc.Cell(cell_id=29, fill=mat1)
cell29.region = +surf107 & -surf108 & +surf40 & -surf39

# NatU
cell30 = openmc.Cell(cell_id=30, fill=mat3)
cell30.region = +surf107 & -surf108 & +surf41 & -surf40

# NatU
cell31 = openmc.Cell(cell_id=31, fill=mat3)
cell31.region = +surf107 & -surf108 & +surf43 & -surf42

# HEU
cell32 = openmc.Cell(cell_id=32, fill=mat1)
cell32.region = +surf107 & -surf108 & +surf44 & -surf43

# NatU
cell33 = openmc.Cell(cell_id=33, fill=mat3)
cell33.region = +surf107 & -surf108 & +surf45 & -surf44

# NatU
cell34 = openmc.Cell(cell_id=34, fill=mat3)
cell34.region = +surf107 & -surf108 & +surf47 & -surf46

# HEU
cell35 = openmc.Cell(cell_id=35, fill=mat1)
cell35.region = +surf105 & -surf108 & +surf49 & -surf48

# NatU
cell36 = openmc.Cell(cell_id=36, fill=mat3)
cell36.region = +surf105 & -surf108 & +surf50 & -surf49

# NatU
cell37 = openmc.Cell(cell_id=37, fill=mat3)
cell37.region = +surf105 & -surf108 & +surf52 & -surf51

# HEU
cell38 = openmc.Cell(cell_id=38, fill=mat1)
cell38.region = +surf105 & -surf108 & +surf53 & -surf52

# NatU
cell39 = openmc.Cell(cell_id=39, fill=mat3)
cell39.region = +surf105 & -surf108 & +surf54 & -surf53

# NatU
cell40 = openmc.Cell(cell_id=40, fill=mat3)
cell40.region = +surf105 & -surf108 & +surf56 & -surf55

# HEU
cell41 = openmc.Cell(cell_id=41, fill=mat1)
cell41.region = +surf105 & -surf108 & +surf57 & -surf56

# NatU
cell42 = openmc.Cell(cell_id=42, fill=mat3)
cell42.region = +surf105 & -surf108 & +surf58 & -surf57

# NatU
cell43 = openmc.Cell(cell_id=43, fill=mat3)
cell43.region = +surf105 & -surf108 & +surf60 & -surf59

# HEU
cell44 = openmc.Cell(cell_id=44, fill=mat1)
cell44.region = +surf105 & -surf108 & +surf61 & -surf60

# NatU
cell45 = openmc.Cell(cell_id=45, fill=mat3)
cell45.region = +surf105 & -surf108 & +surf62 & -surf61

# NatU
cell46 = openmc.Cell(cell_id=46, fill=mat3)
cell46.region = +surf105 & -surf108 & +surf64 & -surf63

# HEU
cell47 = openmc.Cell(cell_id=47, fill=mat1)
cell47.region = +surf105 & -surf108 & +surf65 & -surf64

# NatU
cell48 = openmc.Cell(cell_id=48, fill=mat3)
cell48.region = +surf105 & -surf108 & +surf66 & -surf65

# NatU
cell49 = openmc.Cell(cell_id=49, fill=mat3)
cell49.region = +surf105 & -surf108 & +surf68 & -surf67

# HEU
cell50 = openmc.Cell(cell_id=50, fill=mat1)
cell50.region = +surf105 & -surf108 & +surf69 & -surf68

# NatUp6
cell51 = openmc.Cell(cell_id=51, fill=mat5)
cell51.region = +surf105 & -surf108 & +surf71 & -surf69

# SS304
cell52 = openmc.Cell(cell_id=52, fill=mat7)
cell52.region = +surf101 & -surf102 & +surf17 & -surf1

# U10
cell53 = openmc.Cell(cell_id=53, fill=mat2)
cell53.region = -surf105 & +surf70 & -surf48

# U10
cell54 = openmc.Cell(cell_id=54, fill=mat2)
cell54.region = -surf107 & +surf48 & -surf18

# U10
cell55 = openmc.Cell(cell_id=55, fill=mat2)
cell55.region = -surf106 & +surf18 & -surf17

# U10
cell56 = openmc.Cell(cell_id=56, fill=mat2)
cell56.region = +surf103 & -surf106 & +surf17 & -surf16

# U10
cell57 = openmc.Cell(cell_id=57, fill=mat2)
cell57.region = +surf103 & -surf107 & +surf16 & -surf14

# U10
cell58 = openmc.Cell(cell_id=58, fill=mat2)
cell58.region = -surf100 & +surf17 & -surf1

# U10
cell59 = openmc.Cell(cell_id=59, fill=mat2)
cell59.region = +surf102 & -surf103 & +surf17 & -surf3

# U10
cell60 = openmc.Cell(cell_id=60, fill=mat2)
cell60.region = +surf102 & -surf104 & +surf3 & -surf2

# D38
cell61 = openmc.Cell(cell_id=61, fill=mat4)
cell61.region = -surf105 & +surf72 & -surf70

# D38
cell62 = openmc.Cell(cell_id=62, fill=mat4)
cell62.region = +surf105 & -surf108 & +surf72 & -surf71

# D38
cell63 = openmc.Cell(cell_id=63, fill=mat4)
cell63.region = +surf108 & -surf109 & +surf72 & -surf1 & +surf112 & +surf115 & +surf118 & +surf121 & +surf124 & +surf127

# D38
cell64 = openmc.Cell(cell_id=64, fill=mat4)
cell64.region = +surf105 & -surf108 & +surf3 & -surf1

# D38
cell65 = openmc.Cell(cell_id=65, fill=mat4)
cell65.region = +surf102 & -surf105 & +surf2 & -surf1

# D38
cell66 = openmc.Cell(cell_id=66, fill=mat4)
cell66.region = +surf104 & -surf105 & +surf3 & -surf2

# D38
cell67 = openmc.Cell(cell_id=67, fill=mat4)
cell67.region = -surf110 & +surf72 & -surf1

# SS347
cell68 = openmc.Cell(cell_id=68, fill=mat6)
cell68.region = +surf111 & -surf112 & +surf72 & -surf1

# D38
cell69 = openmc.Cell(cell_id=69, fill=mat4)
cell69.region = -surf113 & +surf72 & -surf1

# SS347
cell70 = openmc.Cell(cell_id=70, fill=mat6)
cell70.region = +surf114 & -surf115 & +surf72 & -surf1

# D38
cell71 = openmc.Cell(cell_id=71, fill=mat4)
cell71.region = -surf116 & +surf72 & -surf1

# SS347
cell72 = openmc.Cell(cell_id=72, fill=mat6)
cell72.region = +surf117 & -surf118 & +surf72 & -surf1

# D38
cell73 = openmc.Cell(cell_id=73, fill=mat4)
cell73.region = -surf119 & +surf72 & -surf1

# SS347
cell74 = openmc.Cell(cell_id=74, fill=mat6)
cell74.region = +surf120 & -surf121 & +surf72 & -surf1

# D38
cell75 = openmc.Cell(cell_id=75, fill=mat4)
cell75.region = -surf122 & +surf72 & -surf1

# SS347
cell76 = openmc.Cell(cell_id=76, fill=mat6)
cell76.region = +surf123 & -surf124 & +surf72 & -surf1

# D38
cell77 = openmc.Cell(cell_id=77, fill=mat4)
cell77.region = -surf125 & +surf72 & -surf1

# SS347
cell78 = openmc.Cell(cell_id=78, fill=mat6)
cell78.region = +surf126 & -surf127 & +surf72 & -surf1

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24, cell25, cell26, cell27, cell28, cell29, cell30, cell31, cell32, cell33, cell34, cell35, cell36, cell37, cell38, cell39, cell40, cell41, cell42, cell43, cell44, cell45, cell46, cell47, cell48, cell49, cell50, cell51, cell52, cell53, cell54, cell55, cell56, cell57, cell58, cell59, cell60, cell61, cell62, cell63, cell64, cell65, cell66, cell67, cell68, cell69, cell70, cell71, cell72, cell73, cell74, cell75, cell76, cell77, cell78])
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
source.space = openmc.stats.Box((-20.7, -20.7, -40.6), (20.7, 20.7, 33.3))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
