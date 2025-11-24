"""
HEU-MET-FAST-001:  GODIVA (52.420 kg U(93.71) @ 18.74 g/cc)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 4.918400e-04)
mat1.add_nuclide("U235", 4.499400e-02)
mat1.add_nuclide("U238", 2.498400e-03)

materials = openmc.Materials([mat1])

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=8.7407, boundary_type="vacuum")

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# CORE
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
