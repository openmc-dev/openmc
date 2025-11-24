"""
LST016-3: Water reflected 28-cm-thick slab tank with 10% enriched uranyl nitrate solution @ H/X=608 (Run 125)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 125
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 7.655500e-07)
mat1.add_nuclide("U235", 9.499900e-05)
mat1.add_nuclide("U236", 9.488100e-08)
mat1.add_nuclide("U238", 8.461700e-04)
mat1.add_nuclide("H1", 5.780000e-02)
mat1.add_element("N", 2.365800e-03)
mat1.add_nuclide("O16", 3.764100e-02)
mat1.add_s_alpha_beta("c_H_in_H2O")

# Stainless Steel
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("C", 7.156700e-05)
mat2.add_element("Si", 7.141500e-04)
mat2.add_element("Mn", 9.909500e-04)
mat2.add_element("P", 5.087900e-05)
mat2.add_element("S", 1.042400e-05)
mat2.add_element("Ni", 8.560000e-03)
mat2.add_element("Cr", 1.672500e-02)
mat2.add_element("Fe", 5.956000e-02)

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

# Slab tank/inner
surf1 = openmc.model.RectangularParallelepiped(-14.04, 14.04, -34.515, 34.515, 0.0, 149.75)
# Slab tank/outer
surf2 = openmc.model.RectangularParallelepiped(-16.57, 16.57, -37.045, 37.045, -2.039999999999992, 152.63)
# Water reflector
surf3 = openmc.model.RectangularParallelepiped(-46.57, 46.57, -67.045, 67.045, -32.03999999999999, 172.63, boundary_type="vacuum")
# Hc
surf4 = openmc.ZPlane(surface_id=4, z0=51.37)

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

# Water
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
source.space = openmc.stats.Point((0.0, 0.0, 25.685))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
