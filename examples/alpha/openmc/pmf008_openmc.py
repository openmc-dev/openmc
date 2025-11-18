"""
PU-MET-FAST-008: THOR (9.587 kg Pu(5.10)-1.01Ga 9 15.29 g/cc in 24.57 cm Th @ 11.58 g/cc)
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1, name="")
mat1.set_density("atom/b-cm", 3.945359e-02)
mat1.add_nuclide("Pu239", 3.6049e-2)
mat1.add_nuclide("Pu240", 1.9562e-3)
mat1.add_nuclide("Pu241", 1.1459e-4)
mat1.add_nuclide("Ga", 1.3338e-3)

mat2 = openmc.Material(material_id=2, name="")
mat2.set_density("atom/b-cm", 3.005400e-02)
mat2.add_nuclide("Th232", 3.0054e-2)

materials = openmc.Materials([mat1, mat2])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=5.310)

# THK = 24.57 cm
surf2 = openmc.Sphere(surface_id=2, r=29.880)


# Cell: core
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="core")
cell0.region = -surf1

# Cell: refl
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="refl")
cell1.region = +surf1 & -surf2

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1])
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
