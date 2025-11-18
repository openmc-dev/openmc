"""
IEU-MET-FAST-010: ZPR-6/9 Loading 11; U9 Benchmark Model
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1, name="")
mat1.set_density("atom/b-cm", 3.879827e-02)
mat1.add_nuclide("U235", 3.48543e-3)
mat1.add_nuclide("U238", 3.52634e-2)
mat1.add_nuclide("U234", 3.34211e-5)
mat1.add_nuclide("U236", 1.60227e-5)

mat2 = openmc.Material(material_id=2, name="")
mat2.set_density("atom/b-cm", 3.891287e-02)
mat2.add_nuclide("U235", 3.47666e-3)
mat2.add_nuclide("U238", 3.53869e-2)
mat2.add_nuclide("U234", 3.33328e-5)
mat2.add_nuclide("U236", 1.59801e-5)

mat3 = openmc.Material(material_id=3, name="")
mat3.set_density("atom/b-cm", 3.879827e-02)
mat3.add_nuclide("U235", 3.48543e-3)
mat3.add_nuclide("U238", 3.52634e-2)
mat3.add_nuclide("U234", 3.34211e-5)
mat3.add_nuclide("U236", 1.60227e-5)

mat4 = openmc.Material(material_id=4, name="")
mat4.set_density("atom/b-cm", 3.891287e-02)
mat4.add_nuclide("U235", 3.47666e-3)
mat4.add_nuclide("U238", 3.53869e-2)
mat4.add_nuclide("U234", 3.33328e-5)
mat4.add_nuclide("U236", 1.59801e-5)

mat5 = openmc.Material(material_id=5, name="")
mat5.set_density("atom/b-cm", 3.974255e-02)
mat5.add_nuclide("U235", 8.31425e-5)
mat5.add_nuclide("U238", 3.72293e-2)
mat5.add_nuclide("Cr", 1.72627e-3)
mat5.add_nuclide("Ni", 7.03836e-4)

mat6 = openmc.Material(material_id=6, name="")
mat6.set_density("atom/b-cm", 3.942384e-02)
mat6.add_nuclide("U235", 8.24312e-5)
mat6.add_nuclide("U238", 3.69113e-2)
mat6.add_nuclide("Cr", 1.72627e-3)
mat6.add_nuclide("Ni", 7.03836e-4)

mat7 = openmc.Material(material_id=7, name="")
mat7.set_density("atom/b-cm", 4.084171e-02)
mat7.add_nuclide("U235", 8.56868e-5)
mat7.add_nuclide("U238", 3.83747e-2)
mat7.add_nuclide("Cr", 1.69502e-3)
mat7.add_nuclide("Ni", 6.86302e-4)

mat8 = openmc.Material(material_id=8, name="")
mat8.set_density("atom/b-cm", 5.967147e-03)
mat8.add_nuclide("Cr", 1.18909e-3)
mat8.add_nuclide("Ni", 4.80175e-4)
mat8.add_nuclide("Fe", 4.27914e-3)
mat8.add_nuclide("C0", 1.87415e-5)

mat9 = openmc.Material(material_id=9, name="")
mat9.set_density("atom/b-cm", 5.919372e-03)
mat9.add_nuclide("Cr", 1.17957e-3)
mat9.add_nuclide("Ni", 4.76329e-4)
mat9.add_nuclide("Fe", 4.24488e-3)
mat9.add_nuclide("C0", 1.85933e-5)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================


# Cell: CORE1
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="CORE1")
cell0.region = +surf6 & -surf7 & -surf12

# Cell: CORE2
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="CORE2")
cell1.region = +surf7 & -surf8 & -surf12

# Cell: CORE3
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="CORE3")
cell2.region = +surf5 & -surf6 & -surf12

# Cell: CORE4
cell3 = openmc.Cell(cell_id=3, fill=mat4, name="CORE4")
cell3.region = +surf4 & -surf5 & -surf12

# Cell: DUUAX
cell4 = openmc.Cell(cell_id=4, fill=mat5, name="DUUAX")
cell4.region = +surf8 & -surf9 & -surf12

# Cell: LDLAX
cell5 = openmc.Cell(cell_id=5, fill=mat6, name="LDLAX")
cell5.region = +surf3 & -surf4 & -surf12

# Cell: DURAD
cell6 = openmc.Cell(cell_id=6, fill=mat7, name="DURAD")
cell6.region = +surf2 & -surf10 & +surf12 & -surf13

# Cell: MATRX
cell7 = openmc.Cell(cell_id=7, fill=mat8, name="MATRX")
cell7.region = +surf1 & -surf3 & -surf12

# Cell: MATRX
cell8 = openmc.Cell(cell_id=8, fill=mat8, name="MATRX")
cell8.region = +surf9 & -surf11 & -surf12

# Cell: MATRX
cell9 = openmc.Cell(cell_id=9, fill=mat8, name="MATRX")
cell9.region = +surf1 & -surf2 & +surf12 & -surf13

# Cell: MATRX
cell10 = openmc.Cell(cell_id=10, fill=mat8, name="MATRX")
cell10.region = +surf10 & -surf11 & +surf12 & -surf13

# Cell: MTRIX
cell11 = openmc.Cell(cell_id=11, fill=mat9, name="MTRIX")
cell11.region = +surf1 & -surf11 & +surf13 & -surf14

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary
# TODO: Adjust dimensions to encompass your entire geometry
boundary_box = openmc.model.RectangularParallelepiped(
    -200, 200, -200, 200, -200, 200,  # xmin, xmax, ymin, ymax, zmin, zmax
    boundary_type="vacuum")

# Create outer void cell (everything outside geometry but inside boundary)
# Particles are killed at the vacuum boundary
outer_region = -boundary_box
outer_region = outer_region & ~cell0.region
outer_region = outer_region & ~cell1.region
outer_region = outer_region & ~cell2.region
outer_region = outer_region & ~cell3.region
outer_region = outer_region & ~cell4.region
outer_region = outer_region & ~cell5.region
outer_region = outer_region & ~cell6.region
outer_region = outer_region & ~cell7.region
outer_region = outer_region & ~cell8.region
outer_region = outer_region & ~cell9.region
outer_region = outer_region & ~cell10.region
outer_region = outer_region & ~cell11.region
outer_cell = openmc.Cell(cell_id=12, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, outer_cell])
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

# Enable delayed neutron kinetics and alpha eigenvalue calculations
settings.calculate_prompt_k = True
settings.calculate_alpha = True

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
