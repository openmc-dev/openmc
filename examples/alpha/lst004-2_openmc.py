"""
LST004-2: Water-reflected 60-cm-diameter cylindrical tank with 10% enriched uranyl nitrate solution @ H/X=771 (Run 29)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 29
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 5.977800e-07)
mat1.add_nuclide("U235", 7.418100e-05)
mat1.add_nuclide("U236", 7.408800e-08)
mat1.add_nuclide("U238", 6.607400e-04)
mat1.add_nuclide("H1", 5.721600e-02)
mat1.add_element("N", 2.814100e-03)
mat1.add_nuclide("O16", 3.785000e-02)
mat1.add_s_alpha_beta("c_H_in_H2O")

# Stainless Steel
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("C", 4.373600e-05)
mat2.add_element("Si", 1.062700e-03)
mat2.add_element("Mn", 1.156100e-03)
mat2.add_element("P", 4.317000e-05)
mat2.add_element("S", 2.978200e-06)
mat2.add_element("Ni", 8.340300e-03)
mat2.add_element("Cr", 1.677500e-02)
mat2.add_element("Fe", 5.942100e-02)

# Water
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("H1", 6.665800e-02)
mat3.add_nuclide("O16", 3.332900e-02)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Air
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_element("N", 3.901600e-05)
mat4.add_nuclide("O16", 1.040900e-05)

materials = openmc.Materials([mat1, mat2, mat3, mat4])

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=29.5)
surf2 = openmc.ZCylinder(surface_id=2, r=29.8)
surf3 = openmc.ZCylinder(surface_id=3, r=59.8, boundary_type="vacuum")
# Hc
surf4 = openmc.ZPlane(surface_id=4, z0=46.70)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Soln
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1 & -surf4

# Air
cell2 = openmc.Cell(cell_id=2, fill=mat4)
cell2.region = -surf1 & +surf4

# SST
cell3 = openmc.Cell(cell_id=3, fill=mat2)
cell3.region = +surf1 & -surf2

# H2O
cell4 = openmc.Cell(cell_id=4, fill=mat3)
cell4.region = +surf2 & -surf3

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4])
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
source.space = openmc.stats.Point((0.0, 0.0, 23.35))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
