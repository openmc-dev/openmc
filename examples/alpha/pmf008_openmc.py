"""
PU-MET-FAST-008: THOR (9.587 kg Pu(5.10)-1.01Ga 9 15.29 g/cc in 24.57 cm Th @ 11.58 g/cc)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("Pu239", 3.604900e-02)
mat1.add_nuclide("Pu240", 1.956200e-03)
mat1.add_nuclide("Pu241", 1.145900e-04)
mat1.add_element("Ga", 1.333800e-03)

mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("Th232", 3.005400e-02)

materials = openmc.Materials([mat1, mat2])

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=5.310)
# THK = 24.57 cm
surf2 = openmc.Sphere(surface_id=2, r=29.880, boundary_type="vacuum")

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# core
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1

# refl
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
