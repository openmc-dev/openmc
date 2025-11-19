"""
IEU-MET-FAST-007 (Rev. 2): Big Ten Detailed Benchmark Model (translated)
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Highly Enriched Uranium (93 wt.%)
mat1 = openmc.Material(material_id=1, name="Highly Enriched Uranium (93 wt.%)")
mat1.set_density("atom/b-cm", 4.827010e-02)
mat1.add_nuclide("U234", 4.9814e-4)
mat1.add_nuclide("U235", 4.5034e-2)
mat1.add_nuclide("U236", 1.3236e-4)
mat1.add_nuclide("U238", 2.6056e-3)

# Intermediate Enriched Uranium (10 wt.%)
mat2 = openmc.Material(material_id=2, name="Intermediate Enriched Uranium (10 wt.%)")
mat2.set_density("atom/b-cm", 4.760921e-02)
mat2.add_nuclide("U234", 2.4761e-5)
mat2.add_nuclide("U235", 4.8461e-3)
mat2.add_nuclide("U236", 4.3348e-5)
mat2.add_nuclide("U238", 4.2695e-2)

# Natural Uranium
mat3 = openmc.Material(material_id=3, name="Natural Uranium")
mat3.set_density("atom/b-cm", 4.819566e-02)
mat3.add_nuclide("U234", 2.6518e-6)
mat3.add_nuclide("U235", 3.4701e-4)
mat3.add_nuclide("U238", 4.7846e-2)

# Depleted Uranium
mat4 = openmc.Material(material_id=4, name="Depleted Uranium")
mat4.set_density("atom/b-cm", 4.777901e-02)
mat4.add_nuclide("U234", 2.8672e-7)
mat4.add_nuclide("U235", 1.0058e-4)
mat4.add_nuclide("U236", 1.1468e-6)
mat4.add_nuclide("U238", 4.7677e-2)

# Natural Uranium (P66)
mat5 = openmc.Material(material_id=5, name="Natural Uranium (P66)")
mat5.set_density("atom/b-cm", 4.808688e-02)
mat5.add_nuclide("U234", 2.6458e-6)
mat5.add_nuclide("U235", 3.4623e-4)
mat5.add_nuclide("U238", 4.7738e-2)

# Steel 347
mat6 = openmc.Material(material_id=6, name="Steel 347")
mat6.set_density("atom/b-cm", 8.525950e-02)
mat6.add_element("Fe", 5.7798e-2)
mat6.add_element("Cr", 1.6678e-2)
mat6.add_element("Ni", 9.0296e-3)
mat6.add_element("Mn", 1.7539e-3)

# Steel 304
mat7 = openmc.Material(material_id=7, name="Steel 304")
mat7.set_density("atom/b-cm", 8.625890e-02)
mat7.add_element("Fe", 5.9308e-2)
mat7.add_element("Cr", 1.7604e-2)
mat7.add_element("Ni", 7.5930e-3)
mat7.add_element("Mn", 1.7539e-3)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================


# Cell: HEU
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="HEU")
cell0.region = +surf103 & -surf108 & +surf5 & -surf4

# Cell: HEU
cell1 = openmc.Cell(cell_id=1, fill=mat1, name="HEU")
cell1.region = +surf103 & -surf108 & +surf7 & -surf6

# Cell: HEU
cell2 = openmc.Cell(cell_id=2, fill=mat1, name="HEU")
cell2.region = +surf103 & -surf108 & +surf9 & -surf8

# Cell: HEU
cell3 = openmc.Cell(cell_id=3, fill=mat1, name="HEU")
cell3.region = +surf103 & -surf108 & +surf11 & -surf10

# Cell: HEU
cell4 = openmc.Cell(cell_id=4, fill=mat1, name="HEU")
cell4.region = +surf103 & -surf108 & +surf13 & -surf12

# Cell: HEU
cell5 = openmc.Cell(cell_id=5, fill=mat1, name="HEU")
cell5.region = +surf107 & -surf108 & +surf16 & -surf15

# Cell: HEU
cell6 = openmc.Cell(cell_id=6, fill=mat1, name="HEU")
cell6.region = +surf107 & -surf108 & +surf20 & -surf19

# Cell: HEU
cell7 = openmc.Cell(cell_id=7, fill=mat1, name="HEU")
cell7.region = +surf107 & -surf108 & +surf24 & -surf23

# Cell: HEU
cell8 = openmc.Cell(cell_id=8, fill=mat1, name="HEU")
cell8.region = +surf107 & -surf108 & +surf28 & -surf27

# Cell: HEU
cell9 = openmc.Cell(cell_id=9, fill=mat1, name="HEU")
cell9.region = +surf107 & -surf108 & +surf32 & -surf31

# Cell: HEU
cell10 = openmc.Cell(cell_id=10, fill=mat1, name="HEU")
cell10.region = +surf107 & -surf108 & +surf36 & -surf35

# Cell: HEU
cell11 = openmc.Cell(cell_id=11, fill=mat1, name="HEU")
cell11.region = +surf107 & -surf108 & +surf40 & -surf39

# Cell: HEU
cell12 = openmc.Cell(cell_id=12, fill=mat1, name="HEU")
cell12.region = +surf107 & -surf108 & +surf44 & -surf43

# Cell: HEU
cell13 = openmc.Cell(cell_id=13, fill=mat1, name="HEU")
cell13.region = +surf105 & -surf108 & +surf49 & -surf48

# Cell: HEU
cell14 = openmc.Cell(cell_id=14, fill=mat1, name="HEU")
cell14.region = +surf105 & -surf108 & +surf53 & -surf52

# Cell: HEU
cell15 = openmc.Cell(cell_id=15, fill=mat1, name="HEU")
cell15.region = +surf105 & -surf108 & +surf57 & -surf56

# Cell: HEU
cell16 = openmc.Cell(cell_id=16, fill=mat1, name="HEU")
cell16.region = +surf105 & -surf108 & +surf61 & -surf60

# Cell: HEU
cell17 = openmc.Cell(cell_id=17, fill=mat1, name="HEU")
cell17.region = +surf105 & -surf108 & +surf65 & -surf64

# Cell: HEU
cell18 = openmc.Cell(cell_id=18, fill=mat1, name="HEU")
cell18.region = +surf105 & -surf108 & +surf69 & -surf68

# Cell: U10
cell19 = openmc.Cell(cell_id=19, fill=mat2, name="U10")
cell19.region = -surf105 & +surf70 & -surf48

# Cell: U10
cell20 = openmc.Cell(cell_id=20, fill=mat2, name="U10")
cell20.region = -surf107 & +surf48 & -surf18

# Cell: U10
cell21 = openmc.Cell(cell_id=21, fill=mat2, name="U10")
cell21.region = -surf106 & +surf18 & -surf17

# Cell: U10
cell22 = openmc.Cell(cell_id=22, fill=mat2, name="U10")
cell22.region = +surf103 & -surf106 & +surf17 & -surf16

# Cell: U10
cell23 = openmc.Cell(cell_id=23, fill=mat2, name="U10")
cell23.region = +surf103 & -surf107 & +surf16 & -surf14

# Cell: U10
cell24 = openmc.Cell(cell_id=24, fill=mat2, name="U10")
cell24.region = -surf100 & +surf17 & -surf1

# Cell: U10
cell25 = openmc.Cell(cell_id=25, fill=mat2, name="U10")
cell25.region = +surf102 & -surf103 & +surf17 & -surf3

# Cell: U10
cell26 = openmc.Cell(cell_id=26, fill=mat2, name="U10")
cell26.region = +surf102 & -surf104 & +surf3 & -surf2

# Cell: NatU
cell27 = openmc.Cell(cell_id=27, fill=mat3, name="NatU")
cell27.region = +surf103 & -surf108 & +surf4 & -surf3

# Cell: NatU
cell28 = openmc.Cell(cell_id=28, fill=mat3, name="NatU")
cell28.region = +surf103 & -surf108 & +surf6 & -surf5

# Cell: NatU
cell29 = openmc.Cell(cell_id=29, fill=mat3, name="NatU")
cell29.region = +surf103 & -surf108 & +surf8 & -surf7

# Cell: NatU
cell30 = openmc.Cell(cell_id=30, fill=mat3, name="NatU")
cell30.region = +surf103 & -surf108 & +surf10 & -surf9

# Cell: NatU
cell31 = openmc.Cell(cell_id=31, fill=mat3, name="NatU")
cell31.region = +surf103 & -surf108 & +surf12 & -surf11

# Cell: NatU
cell32 = openmc.Cell(cell_id=32, fill=mat3, name="NatU")
cell32.region = +surf107 & -surf108 & +surf15 & -surf13

# Cell: NatU
cell33 = openmc.Cell(cell_id=33, fill=mat3, name="NatU")
cell33.region = +surf107 & -surf108 & +surf21 & -surf20

# Cell: NatU
cell34 = openmc.Cell(cell_id=34, fill=mat3, name="NatU")
cell34.region = +surf107 & -surf108 & +surf23 & -surf22

# Cell: NatU
cell35 = openmc.Cell(cell_id=35, fill=mat3, name="NatU")
cell35.region = +surf107 & -surf108 & +surf25 & -surf24

# Cell: NatU
cell36 = openmc.Cell(cell_id=36, fill=mat3, name="NatU")
cell36.region = +surf107 & -surf108 & +surf27 & -surf26

# Cell: NatU
cell37 = openmc.Cell(cell_id=37, fill=mat3, name="NatU")
cell37.region = +surf107 & -surf108 & +surf29 & -surf28

# Cell: NatU
cell38 = openmc.Cell(cell_id=38, fill=mat3, name="NatU")
cell38.region = +surf107 & -surf108 & +surf31 & -surf30

# Cell: NatU
cell39 = openmc.Cell(cell_id=39, fill=mat3, name="NatU")
cell39.region = +surf107 & -surf108 & +surf33 & -surf32

# Cell: NatU
cell40 = openmc.Cell(cell_id=40, fill=mat3, name="NatU")
cell40.region = +surf107 & -surf108 & +surf35 & -surf34

# Cell: NatU
cell41 = openmc.Cell(cell_id=41, fill=mat3, name="NatU")
cell41.region = +surf107 & -surf108 & +surf37 & -surf36

# Cell: NatU
cell42 = openmc.Cell(cell_id=42, fill=mat3, name="NatU")
cell42.region = +surf107 & -surf108 & +surf39 & -surf38

# Cell: NatU
cell43 = openmc.Cell(cell_id=43, fill=mat3, name="NatU")
cell43.region = +surf107 & -surf108 & +surf41 & -surf40

# Cell: NatU
cell44 = openmc.Cell(cell_id=44, fill=mat3, name="NatU")
cell44.region = +surf107 & -surf108 & +surf43 & -surf42

# Cell: NatU
cell45 = openmc.Cell(cell_id=45, fill=mat3, name="NatU")
cell45.region = +surf107 & -surf108 & +surf45 & -surf44

# Cell: NatU
cell46 = openmc.Cell(cell_id=46, fill=mat3, name="NatU")
cell46.region = +surf107 & -surf108 & +surf47 & -surf46

# Cell: NatU
cell47 = openmc.Cell(cell_id=47, fill=mat3, name="NatU")
cell47.region = +surf105 & -surf108 & +surf50 & -surf49

# Cell: NatU
cell48 = openmc.Cell(cell_id=48, fill=mat3, name="NatU")
cell48.region = +surf105 & -surf108 & +surf52 & -surf51

# Cell: NatU
cell49 = openmc.Cell(cell_id=49, fill=mat3, name="NatU")
cell49.region = +surf105 & -surf108 & +surf54 & -surf53

# Cell: NatU
cell50 = openmc.Cell(cell_id=50, fill=mat3, name="NatU")
cell50.region = +surf105 & -surf108 & +surf56 & -surf55

# Cell: NatU
cell51 = openmc.Cell(cell_id=51, fill=mat3, name="NatU")
cell51.region = +surf105 & -surf108 & +surf58 & -surf57

# Cell: NatU
cell52 = openmc.Cell(cell_id=52, fill=mat3, name="NatU")
cell52.region = +surf105 & -surf108 & +surf60 & -surf59

# Cell: NatU
cell53 = openmc.Cell(cell_id=53, fill=mat3, name="NatU")
cell53.region = +surf105 & -surf108 & +surf62 & -surf61

# Cell: NatU
cell54 = openmc.Cell(cell_id=54, fill=mat3, name="NatU")
cell54.region = +surf105 & -surf108 & +surf64 & -surf63

# Cell: NatU
cell55 = openmc.Cell(cell_id=55, fill=mat3, name="NatU")
cell55.region = +surf105 & -surf108 & +surf66 & -surf65

# Cell: NatU
cell56 = openmc.Cell(cell_id=56, fill=mat3, name="NatU")
cell56.region = +surf105 & -surf108 & +surf68 & -surf67

# Cell: D38
cell57 = openmc.Cell(cell_id=57, fill=mat4, name="D38")
cell57.region = +surf106 & -surf108 & +surf18 & -surf16

# Cell: D38
cell58 = openmc.Cell(cell_id=58, fill=mat4, name="D38")
cell58.region = -surf105 & +surf72 & -surf70

# Cell: D38
cell59 = openmc.Cell(cell_id=59, fill=mat4, name="D38")
cell59.region = +surf105 & -surf108 & +surf72 & -surf71

# Cell: D38
cell60 = openmc.Cell(cell_id=60, fill=mat4, name="D38")
cell60.region = +surf108 & -surf109 & +surf72 & -surf1 & +surf112 & +surf115 & +surf118 & +surf121 & +surf124 & +surf127

# Cell: D38
cell61 = openmc.Cell(cell_id=61, fill=mat4, name="D38")
cell61.region = +surf105 & -surf108 & +surf3 & -surf1

# Cell: D38
cell62 = openmc.Cell(cell_id=62, fill=mat4, name="D38")
cell62.region = +surf102 & -surf105 & +surf2 & -surf1

# Cell: D38
cell63 = openmc.Cell(cell_id=63, fill=mat4, name="D38")
cell63.region = +surf104 & -surf105 & +surf3 & -surf2

# Cell: D38
cell64 = openmc.Cell(cell_id=64, fill=mat4, name="D38")
cell64.region = -surf110 & +surf72 & -surf1

# Cell: D38
cell65 = openmc.Cell(cell_id=65, fill=mat4, name="D38")
cell65.region = -surf113 & +surf72 & -surf1

# Cell: D38
cell66 = openmc.Cell(cell_id=66, fill=mat4, name="D38")
cell66.region = -surf116 & +surf72 & -surf1

# Cell: D38
cell67 = openmc.Cell(cell_id=67, fill=mat4, name="D38")
cell67.region = -surf119 & +surf72 & -surf1

# Cell: D38
cell68 = openmc.Cell(cell_id=68, fill=mat4, name="D38")
cell68.region = -surf122 & +surf72 & -surf1

# Cell: D38
cell69 = openmc.Cell(cell_id=69, fill=mat4, name="D38")
cell69.region = -surf125 & +surf72 & -surf1

# Cell: NatUp6
cell70 = openmc.Cell(cell_id=70, fill=mat5, name="NatUp6")
cell70.region = +surf105 & -surf108 & +surf71 & -surf69

# Cell: SS347
cell71 = openmc.Cell(cell_id=71, fill=mat6, name="SS347")
cell71.region = +surf111 & -surf112 & +surf72 & -surf1

# Cell: SS347
cell72 = openmc.Cell(cell_id=72, fill=mat6, name="SS347")
cell72.region = +surf114 & -surf115 & +surf72 & -surf1

# Cell: SS347
cell73 = openmc.Cell(cell_id=73, fill=mat6, name="SS347")
cell73.region = +surf117 & -surf118 & +surf72 & -surf1

# Cell: SS347
cell74 = openmc.Cell(cell_id=74, fill=mat6, name="SS347")
cell74.region = +surf120 & -surf121 & +surf72 & -surf1

# Cell: SS347
cell75 = openmc.Cell(cell_id=75, fill=mat6, name="SS347")
cell75.region = +surf123 & -surf124 & +surf72 & -surf1

# Cell: SS347
cell76 = openmc.Cell(cell_id=76, fill=mat6, name="SS347")
cell76.region = +surf126 & -surf127 & +surf72 & -surf1

# Cell: SS304
cell77 = openmc.Cell(cell_id=77, fill=mat7, name="SS304")
cell77.region = +surf101 & -surf102 & +surf17 & -surf1

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
outer_region = outer_region & ~cell68.region
outer_region = outer_region & ~cell69.region
outer_region = outer_region & ~cell70.region
outer_region = outer_region & ~cell71.region
outer_region = outer_region & ~cell72.region
outer_region = outer_region & ~cell73.region
outer_region = outer_region & ~cell74.region
outer_region = outer_region & ~cell75.region
outer_region = outer_region & ~cell76.region
outer_region = outer_region & ~cell77.region
outer_cell = openmc.Cell(cell_id=78, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24, cell25, cell26, cell27, cell28, cell29, cell30, cell31, cell32, cell33, cell34, cell35, cell36, cell37, cell38, cell39, cell40, cell41, cell42, cell43, cell44, cell45, cell46, cell47, cell48, cell49, cell50, cell51, cell52, cell53, cell54, cell55, cell56, cell57, cell58, cell59, cell60, cell61, cell62, cell63, cell64, cell65, cell66, cell67, cell68, cell69, cell70, cell71, cell72, cell73, cell74, cell75, cell76, cell77, outer_cell])
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
