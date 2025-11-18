"""
SHE-8 (JAERI-1257) Homogenized model
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

# Core
mat1 = openmc.Material(material_id=1, name="Core")
mat1.set_density("atom/b-cm", 7.966875e-02)
mat1.add_nuclide("U235", 0.3415e-4)
mat1.add_nuclide("U238", 1.379e-4)
mat1.add_nuclide("C0", 0.7911e-1)
mat1.add_nuclide("O16", 0.3867e-3)
mat1.add_s_alpha_beta("c_H_in_H2O")

# Reflector
mat2 = openmc.Material(material_id=2, name="Reflector")
mat2.set_density("atom/b-cm", 7.735673e-02)
mat2.add_nuclide("C0", 0.7732e-1)
mat2.add_nuclide("O16", 0.3673e-4)
mat2.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([mat1, mat2])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=28.73)

# COG surface type "pri" with parameters: 6 150 0 75 129.9038 -75 129.9038 -150 0 -75 -129.9038 75 -129.9038 -120 120 tr 0 0 0 0 0 1 0 1 0
# This surface type requires manual translation to OpenMC
surf2 = openmc.Sphere(surface_id=2, r=1.0)  # PLACEHOLDER - REPLACE THIS


# Cell: Core
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Core")
cell0.region = -surf1 & -surf2

# Cell: Refl
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="Refl")
cell1.region = +surf1 & -surf2

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 100
settings.batches = 520
settings.inactive = 21
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
