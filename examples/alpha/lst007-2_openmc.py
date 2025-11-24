"""
LST007-2: Unreflected 60-cm-diameter cylindrical tank with 10% enriched uranyl nitrate solution @ H/X=770 (Run 30)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 30
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 5.984000e-07)
mat1.add_nuclide("U235", 7.425700e-05)
mat1.add_nuclide("U236", 7.416500e-08)
mat1.add_nuclide("U238", 6.614200e-04)
mat1.add_nuclide("H1", 5.717600e-02)
mat1.add_element("N", 2.815600e-03)
mat1.add_nuclide("O16", 3.783600e-02)
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

# Air
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_element("N", 3.901600e-05)
mat3.add_nuclide("O16", 1.040900e-05)

materials = openmc.Materials([mat1, mat2, mat3])

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=29.5)
surf2 = openmc.ZCylinder(surface_id=2, r=29.8, boundary_type="vacuum")
# Hc
surf3 = openmc.ZPlane(surface_id=3, z0=54.20)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Soln
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1 & -surf3

# Air
cell2 = openmc.Cell(cell_id=2, fill=mat3)
cell2.region = -surf1 & +surf3

# SST
cell3 = openmc.Cell(cell_id=3, fill=mat2)
cell3.region = +surf1 & -surf2

root_universe = openmc.Universe(cells=[cell1, cell2, cell3])
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
source.space = openmc.stats.Point((0.0, 0.0, 27.1))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
