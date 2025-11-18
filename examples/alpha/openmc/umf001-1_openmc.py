"""
U233-MET-FAST-001: 16.535 kg 233U(98.1) @ 18.424 g/cc
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

# per Table 1
mat1 = openmc.Material(material_id=1, name="per Table 1")
mat1.set_density("atom/b-cm", 1.194240e+02)
mat1.add_nuclide("F", 18.424)
mat1.add_nuclide("U233", 98.13)
mat1.add_nuclide("U234", 1.24)
mat1.add_nuclide("U235", 0.03)
mat1.add_nuclide("U238", 0.60)
mat1.add_nuclide("Table", 1)

materials = openmc.Materials([mat1])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=5.9838)


# Cell: U233
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="U233")
cell0.region = -surf1

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 15000
settings.batches = 4400
settings.inactive = 100
settings.run_mode = "eigenvalue"

# Source definition
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)
settings.source = source

settings.export_to_xml()

# ==============================================================================
# Tallies
# ==============================================================================

tallies = openmc.Tallies()
tallies.export_to_xml()

# ==============================================================================
# Run OpenMC
# ==============================================================================

openmc.run()
