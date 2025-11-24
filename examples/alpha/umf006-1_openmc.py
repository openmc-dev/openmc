"""
U233-MET-FAST-006: 5.740 kg 233U(98.1) in FLATTOP
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# U-233 per Table 2
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_element("P", 1.842000e+01)
mat1.add_nuclide("U233", 9.813000e+01)
mat1.add_nuclide("U234", 1.240000e+00)
mat1.add_nuclide("U235", 3.000000e-02)
mat1.add_nuclide("U238", 6.000000e-01)

# Nat-U
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("U234", 2.403960e-06)
mat2.add_nuclide("U235", 3.461702e-04)
mat2.add_nuclide("U238", 4.773062e-02)

materials = openmc.Materials([mat1, mat2])

# ==============================================================================
# Geometry
# ==============================================================================

# Dimensions per
surf1 = openmc.Sphere(surface_id=1, r=4.2058)
# Section 3.2
surf2 = openmc.Sphere(surface_id=2, r=24.1194, boundary_type="vacuum")

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# U233
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1

# Tu
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = +surf1 & -surf2

root_universe = openmc.Universe(cells=[cell1, cell2])
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
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
