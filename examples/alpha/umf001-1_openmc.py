"""
U233-MET-FAST-001: 16.535 kg 233U(98.1) @ 18.424 g/cc
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# per Table 1
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U233", 4.671866e-02)
mat1.add_nuclide("U234", 5.903510e-04)
mat1.add_nuclide("U235", 1.428269e-05)
mat1.add_nuclide("U238", 2.856537e-04)

materials = openmc.Materials([mat1])

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=5.9838, boundary_type="vacuum")

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# U233
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
