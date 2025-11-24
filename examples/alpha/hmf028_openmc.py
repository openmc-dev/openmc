"""
HEU-MET-FAST-028-1:  17.840 kg Oy(93.24) sphere in 7.09" Nat-U (Flattop)
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Oy(93.24) @ 18.62 g/cc
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 4.886900e-04)
mat1.add_nuclide("U235", 4.448200e-02)
mat1.add_nuclide("U238", 2.703800e-03)

# Natural-U @ 19.00 g/cc
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("U234", 2.643800e-06)
mat2.add_nuclide("U235", 3.461000e-04)
mat2.add_nuclide("U238", 4.772100e-02)

materials = openmc.Materials([mat1, mat2])

# ==============================================================================
# Geometry
# ==============================================================================


# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# Oy
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
