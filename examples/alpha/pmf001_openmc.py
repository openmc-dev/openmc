"""
PU-MET-FAST-001: JEZEBEL (17.020 kg Pu(95.48)-1.02Ga @ 15.61 g/cc)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_element("Ga", 1.375200e-03)
mat1.add_nuclide("Pu239", 3.704700e-02)
mat1.add_nuclide("Pu240", 1.751200e-03)
mat1.add_nuclide("Pu241", 1.167400e-04)

materials = openmc.Materials([mat1])

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=6.3849, boundary_type="vacuum")

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# alloy
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1

root_universe = openmc.Universe(cells=[cell1])
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
