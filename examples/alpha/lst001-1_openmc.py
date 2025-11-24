"""
LST001-1: SHEBA-II
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Fuel
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 6.785500e-07)
mat1.add_nuclide("U235", 1.237700e-04)
mat1.add_nuclide("U236", 1.208500e-06)
mat1.add_nuclide("U238", 2.350800e-03)
mat1.add_nuclide("H1", 5.617900e-02)
mat1.add_nuclide("O16", 3.296700e-02)
mat1.add_element("F", 5.103500e-03)
mat1.add_s_alpha_beta("c_H_in_H2O")

# Air
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("N", 3.521400e-05)
mat2.add_nuclide("O16", 1.509200e-05)

# SS304L
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_element("Cr", 1.634800e-02)
mat3.add_element("Mn", 1.719200e-03)
mat3.add_element("Fe", 6.003800e-02)
mat3.add_element("Ni", 7.241800e-03)

materials = openmc.Materials([mat1, mat2, mat3])

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=2.54)
surf2 = openmc.ZCylinder(surface_id=2, r=3.175)
surf3 = openmc.ZCylinder(surface_id=3, r=24.4475)
surf4 = openmc.ZCylinder(surface_id=4, r=25.4, boundary_type="vacuum")
surf5 = openmc.ZPlane(surface_id=5, z0=44.8)

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Soln
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = +surf2 & -surf3 & -surf5

# Air
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = +surf2 & -surf3 & +surf5

# Air
cell3 = openmc.Cell(cell_id=3, fill=mat2)
cell3.region = -surf1 & -surf4

# SS304L
cell4 = openmc.Cell(cell_id=4, fill=mat3)
cell4.region = +surf1 & +surf3 & -surf4

# SS304L
cell5 = openmc.Cell(cell_id=5, fill=mat3)
cell5.region = +surf1 & -surf2 & -surf3

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5])
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
source.space = openmc.stats.Box((-6.0, -6.0, 21.4), (6.0, 6.0, 23.4))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
