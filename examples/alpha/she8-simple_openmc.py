"""
SHE-8 (JAERI-1257) Homogenized model
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Core
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U235", 3.415000e-05)
mat1.add_nuclide("U238", 1.379000e-04)
mat1.add_element("C", 7.911000e-02)
mat1.add_nuclide("H1", 7.191000e-05)
mat1.add_nuclide("O16", 3.867000e-04)
mat1.add_s_alpha_beta("c_Graphite")
mat1.add_s_alpha_beta("c_H_in_H2O")

# Reflector
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("C", 7.732000e-02)
mat2.add_nuclide("H1", 7.137000e-05)
mat2.add_nuclide("O16", 3.673000e-05)
mat2.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([mat1, mat2])

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=28.73)
# surf2: Unsupported surface type "pri" with params ['6', '150', '0', '75', '129.9038', '-75', '129.9038', '-150', '0', '-75', '-129.9038', '75', '-129.9038', '-120', '120', 'tr', '0', '0', '0', '0', '0', '1', '0', '1', '0']

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Core
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1 & -surf2

# Refl
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = +surf1 & -surf2

root_universe = openmc.Universe(cells=[cell1, cell2])
geometry = openmc.Geometry(root_universe)

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 100
settings.batches = 520
settings.inactive = 21
settings.run_mode = "eigenvalue"

source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
