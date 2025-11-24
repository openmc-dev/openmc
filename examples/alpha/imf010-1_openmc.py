"""
IEU-MET-FAST-010: ZPR-6/9 Loading 11; U9 Benchmark Model
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U235", 3.485430e-03)
mat1.add_nuclide("U238", 3.526340e-02)
mat1.add_nuclide("U234", 3.342110e-05)
mat1.add_nuclide("U236", 1.602270e-05)
mat1.add_element("Cr", 1.860670e-03)
mat1.add_element("Ni", 7.490530e-04)
mat1.add_element("Fe", 6.656880e-03)
mat1.add_element("C", 9.629900e-05)
mat1.add_element("Mo", 1.183270e-05)
mat1.add_element("Mn", 1.635090e-04)
mat1.add_element("Cu", 2.388330e-05)
mat1.add_element("H", 1.945800e-05)
mat1.add_element("Si", 7.160010e-05)
mat1.add_element("Cl", 3.335060e-05)
mat1.add_element("F", 9.912200e-05)

mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("U235", 3.476660e-03)
mat2.add_nuclide("U238", 3.538690e-02)
mat2.add_nuclide("U234", 3.333280e-05)
mat2.add_nuclide("U236", 1.598010e-05)
mat2.add_element("Cr", 1.874180e-03)
mat2.add_element("Ni", 7.546540e-04)
mat2.add_element("Fe", 6.707520e-03)
mat2.add_element("C", 9.641660e-05)
mat2.add_element("Mo", 1.197990e-05)
mat2.add_element("Mn", 1.648230e-04)
mat2.add_element("Cu", 2.423290e-05)
mat2.add_element("H", 1.944980e-05)
mat2.add_element("Si", 7.420380e-05)
mat2.add_element("Cl", 3.346860e-05)
mat2.add_element("F", 9.947260e-05)

mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("U235", 3.485430e-03)
mat3.add_nuclide("U238", 3.526340e-02)
mat3.add_nuclide("U234", 3.342110e-05)
mat3.add_nuclide("U236", 1.602270e-05)
mat3.add_element("Cr", 1.860670e-03)
mat3.add_element("Ni", 7.490530e-04)
mat3.add_element("Fe", 6.656880e-03)
mat3.add_element("C", 9.626060e-05)
mat3.add_element("Mo", 1.183270e-05)
mat3.add_element("Mn", 1.635090e-04)
mat3.add_element("Cu", 2.388330e-05)
mat3.add_element("H", 1.944710e-05)
mat3.add_element("Si", 7.160010e-05)
mat3.add_element("Cl", 3.333140e-05)
mat3.add_element("F", 9.906500e-05)

mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_nuclide("U235", 3.476660e-03)
mat4.add_nuclide("U238", 3.538690e-02)
mat4.add_nuclide("U234", 3.333280e-05)
mat4.add_nuclide("U236", 1.598010e-05)
mat4.add_element("Cr", 1.874180e-03)
mat4.add_element("Ni", 7.546540e-04)
mat4.add_element("Fe", 6.707520e-03)
mat4.add_element("C", 9.637810e-05)
mat4.add_element("Mo", 1.197990e-05)
mat4.add_element("Mn", 1.648230e-04)
mat4.add_element("Cu", 2.423290e-05)
mat4.add_element("H", 1.943880e-05)
mat4.add_element("Si", 7.420380e-05)
mat4.add_element("Cl", 3.344930e-05)
mat4.add_element("F", 9.941540e-05)

mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")
mat5.add_nuclide("U235", 8.314250e-05)
mat5.add_nuclide("U238", 3.722930e-02)
mat5.add_element("Cr", 1.726270e-03)
mat5.add_element("Ni", 7.038360e-04)
mat5.add_element("Fe", 6.266560e-03)
mat5.add_element("Al", 4.167470e-04)
mat5.add_element("C", 3.991490e-05)
mat5.add_element("Mo", 1.122060e-05)
mat5.add_element("Mn", 1.468590e-04)
mat5.add_element("Cu", 2.135370e-05)
mat5.add_element("H", 2.854270e-06)
mat5.add_element("Si", 8.312370e-05)
mat5.add_element("Cl", 4.933940e-06)
mat5.add_element("F", 1.461070e-05)

mat6 = openmc.Material(material_id=6)
mat6.set_density("sum")
mat6.add_nuclide("U235", 8.243120e-05)
mat6.add_nuclide("U238", 3.691130e-02)
mat6.add_element("Cr", 1.726270e-03)
mat6.add_element("Ni", 7.038360e-04)
mat6.add_element("Fe", 6.266070e-03)
mat6.add_element("Al", 4.167470e-04)
mat6.add_element("C", 3.976750e-05)
mat6.add_element("Mo", 1.122060e-05)
mat6.add_element("Mn", 1.468590e-04)
mat6.add_element("Cu", 2.135370e-05)
mat6.add_element("H", 2.817770e-06)
mat6.add_element("Si", 8.312370e-05)
mat6.add_element("Cl", 4.871180e-06)
mat6.add_element("F", 1.442490e-05)

mat7 = openmc.Material(material_id=7)
mat7.set_density("sum")
mat7.add_nuclide("U235", 8.568680e-05)
mat7.add_nuclide("U238", 3.837470e-02)
mat7.add_element("Cr", 1.695020e-03)
mat7.add_element("Ni", 6.863020e-04)
mat7.add_element("Fe", 6.132440e-03)
mat7.add_element("C", 3.779350e-05)
mat7.add_element("Mo", 1.076710e-05)
mat7.add_element("Mn", 1.468030e-04)
mat7.add_element("Cu", 2.168200e-05)
mat7.add_element("H", 2.538030e-06)
mat7.add_element("Si", 7.696760e-05)
mat7.add_element("Cl", 4.394540e-06)
mat7.add_element("F", 1.301320e-05)

mat8 = openmc.Material(material_id=8)
mat8.set_density("sum")
mat8.add_element("Cr", 1.189090e-03)
mat8.add_element("Ni", 4.801750e-04)
mat8.add_element("Fe", 4.279140e-03)
mat8.add_element("C", 1.874150e-05)
mat8.add_element("Mo", 8.256530e-06)
mat8.add_element("Mn", 1.059050e-04)
mat8.add_element("Cu", 1.723190e-05)
mat8.add_element("Si", 6.831050e-05)
mat8.add_element("F", 1.301320e-05)

mat9 = openmc.Material(material_id=9)
mat9.set_density("sum")
mat9.add_element("Cr", 1.179570e-03)
mat9.add_element("Ni", 4.763290e-04)
mat9.add_element("Fe", 4.244880e-03)
mat9.add_element("C", 1.859330e-05)
mat9.add_element("Mo", 8.190760e-06)
mat9.add_element("Mn", 1.050570e-04)
mat9.add_element("Cu", 1.709440e-05)
mat9.add_element("Si", 6.776300e-05)
mat9.add_element("F", 1.301320e-05)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9])

# ==============================================================================
# Geometry
# ==============================================================================


# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# CORE1
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = +surf6 & -surf7 & -surf12

# CORE2
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = +surf7 & -surf8 & -surf12

# CORE3
cell3 = openmc.Cell(cell_id=3, fill=mat3)
cell3.region = +surf5 & -surf6 & -surf12

# CORE4
cell4 = openmc.Cell(cell_id=4, fill=mat4)
cell4.region = +surf4 & -surf5 & -surf12

# DUUAX
cell5 = openmc.Cell(cell_id=5, fill=mat5)
cell5.region = +surf8 & -surf9 & -surf12

# LDLAX
cell6 = openmc.Cell(cell_id=6, fill=mat6)
cell6.region = +surf3 & -surf4 & -surf12

# MATRX
cell7 = openmc.Cell(cell_id=7, fill=mat8)
cell7.region = +surf1 & -surf3 & -surf12

# MATRX
cell8 = openmc.Cell(cell_id=8, fill=mat8)
cell8.region = +surf9 & -surf11 & -surf12

# DURAD
cell9 = openmc.Cell(cell_id=9, fill=mat7)
cell9.region = +surf2 & -surf10 & +surf12 & -surf13

# MATRX
cell10 = openmc.Cell(cell_id=10, fill=mat8)
cell10.region = +surf1 & -surf2 & +surf12 & -surf13

# MATRX
cell11 = openmc.Cell(cell_id=11, fill=mat8)
cell11.region = +surf10 & -surf11 & +surf12 & -surf13

# MTRIX
cell12 = openmc.Cell(cell_id=12, fill=mat9)
cell12.region = +surf1 & -surf11 & +surf13 & -surf14

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10, cell11, cell12])
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
