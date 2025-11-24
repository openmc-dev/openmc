"""
PU-MET-INTER-002: ZPR-6/10 Loading 24; Benchmark Model
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("Pu238", 1.693800e-09)
mat1.add_nuclide("Pu239", 1.062990e-03)
mat1.add_nuclide("Pu240", 5.030370e-05)
mat1.add_nuclide("Pu241", 1.700340e-06)
mat1.add_nuclide("Pu242", 1.155920e-07)
mat1.add_nuclide("Am241", 2.986130e-06)
mat1.add_element("Cr", 9.841990e-03)
mat1.add_element("Al", 1.116590e-04)
mat1.add_element("Mn", 8.052540e-04)
mat1.add_element("Co", 3.719240e-06)
mat1.add_element("Ni", 4.085280e-03)
mat1.add_element("C", 2.582130e-02)
mat1.add_element("Cu", 1.012310e-04)
mat1.add_element("Fe", 3.526070e-02)
mat1.add_element("Mo", 7.427920e-05)
mat1.add_element("Si", 2.795770e-04)

mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("Cr", 1.572770e-02)
mat2.add_element("Mn", 1.560810e-03)
mat2.add_element("Co", 3.735820e-06)
mat2.add_element("Ni", 6.946570e-03)
mat2.add_element("C", 2.359480e-04)
mat2.add_element("Cu", 2.981150e-05)
mat2.add_element("Fe", 5.552210e-02)
mat2.add_element("Mo", 1.521200e-05)
mat2.add_element("Si", 9.454900e-04)

mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_element("Cr", 1.562860e-02)
mat3.add_element("Mn", 1.420570e-03)
mat3.add_element("Co", 3.779320e-06)
mat3.add_element("Ni", 6.784100e-03)
mat3.add_element("C", 2.371060e-04)
mat3.add_element("Cu", 5.040520e-05)
mat3.add_element("Fe", 5.519640e-02)
mat3.add_element("Mo", 4.640230e-05)
mat3.add_element("Si", 8.633600e-04)

mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_element("Cr", 1.827580e-03)
mat4.add_element("Mn", 5.503500e-04)
mat4.add_element("Co", 3.772220e-06)
mat4.add_element("Ni", 7.614350e-04)
mat4.add_element("C", 4.273150e-04)
mat4.add_element("Cu", 2.904570e-05)
mat4.add_element("Fe", 7.665320e-02)
mat4.add_element("Mo", 1.462350e-05)
mat4.add_element("Si", 1.451480e-04)

mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")
mat5.add_element("Cr", 1.187270e-03)
mat5.add_element("Mn", 1.056990e-04)
mat5.add_element("Ni", 4.795390e-04)
mat5.add_element("C", 1.870780e-05)
mat5.add_element("Cu", 1.718490e-05)
mat5.add_element("Fe", 4.272430e-03)
mat5.add_element("Mo", 8.236970e-06)
mat5.add_element("Si", 6.817470e-05)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5])

# ==============================================================================
# Geometry
# ==============================================================================


# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# CORE
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = +surf1 & -surf2 & -surf5

# SSAXR
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = +surf2 & -surf3 & -surf5

# SSRDR
cell3 = openmc.Cell(cell_id=3, fill=mat3)
cell3.region = +surf1 & -surf3 & +surf5 & -surf6

# FERDR
cell4 = openmc.Cell(cell_id=4, fill=mat4)
cell4.region = +surf1 & -surf3 & +surf6 & -surf7

# MATRX
cell5 = openmc.Cell(cell_id=5, fill=mat5)
cell5.region = +surf1 & -surf3 & +surf7 & -surf8

# MATRX
cell6 = openmc.Cell(cell_id=6, fill=mat5)
cell6.region = +surf3 & -surf4 & -surf8

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6])
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
source.space = openmc.stats.Point((0.0001, 0.0001, 0.0001))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
