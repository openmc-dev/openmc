"""
LST021-1: Unreflected 80-cm diameter cylindrical tank with 10% enriched uranyl nitrate solution @ H/X=971 (Run 215)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 215
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 5.004200e-07)
mat1.add_nuclide("U235", 6.209800e-05)
mat1.add_nuclide("U236", 6.202100e-08)
mat1.add_nuclide("U238", 5.531200e-04)
mat1.add_nuclide("H1", 6.029700e-02)
mat1.add_element("N", 1.815700e-03)
mat1.add_nuclide("O16", 3.653500e-02)
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
mat3.add_element("N", 3.901400e-05)
mat3.add_nuclide("O16", 1.041000e-05)

materials = openmc.Materials([mat1, mat2, mat3])

# ==============================================================================
# Geometry
# ==============================================================================

# Tank/inner
surf1 = openmc.ZCylinder(surface_id=1, r=39.505)
# Tank/outer
surf2 = openmc.ZCylinder(surface_id=2, r=39.815)
# Base plate
surf3 = openmc.model.RectangularParallelepiped(-60.2, 39.8, -50.0, 50.0, -19.0, -16.0)
# Hole in base plate
surf4 = openmc.ZCylinder(surface_id=4, x0=tr, y0=24.8, r=7.76)
# boundary condition
surf5 = openmc.ZCylinder(surface_id=5, r=79.815, boundary_type="vacuum")
# Hc
surf6 = openmc.ZPlane(surface_id=6, z0=43.98)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Soln
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1 & -surf6

# SST
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = +surf1 & -surf2

# SST
cell3 = openmc.Cell(cell_id=3, fill=mat2)
cell3.region = -surf3 & +surf4

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
source.space = openmc.stats.Point((0.0, 0.0, 21.99))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
