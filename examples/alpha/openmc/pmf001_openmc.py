"""
PU-MET-FAST-001: JEZEBEL (17.020 kg Pu(95.48)-1.02Ga @ 15.61 g/cc)
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1, name="")
mat1.set_density("atom/b-cm", 4.029014e-02)
mat1.add_nuclide("Ga", 1.3752e-3)
mat1.add_nuclide("Pu239", 3.7047e-2)
mat1.add_nuclide("Pu240", 1.7512e-3)
mat1.add_nuclide("Pu241", 1.1674e-4)

materials = openmc.Materials([mat1])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=6.3849)


# Cell: alloy
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="alloy")
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
