"""
HEU-MET-INTER-001: ZPR-9/34 Loading 303; Benchmark Model
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U235", 1.137620e-03)
mat1.add_nuclide("U238", 6.593900e-05)
mat1.add_nuclide("U234", 1.115760e-05)
mat1.add_nuclide("U236", 5.349150e-06)
mat1.add_element("Cr", 2.483970e-03)
mat1.add_element("Ni", 1.096670e-03)
mat1.add_element("Fe", 7.006900e-02)
mat1.add_element("Al", 1.109440e-03)
mat1.add_element("C", 3.275000e-04)
mat1.add_element("Mo", 8.762460e-06)
mat1.add_element("Mn", 5.361540e-04)
mat1.add_element("Cu", 1.812890e-05)
mat1.add_element("H", 1.025620e-06)
mat1.add_element("Si", 1.254780e-04)
mat1.add_element("Mg", 3.194670e-05)
mat1.add_element("Cl", 1.749460e-06)
mat1.add_element("F", 5.305090e-06)

mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_nuclide("U235", 1.132520e-03)
mat2.add_nuclide("U238", 6.564320e-05)
mat2.add_nuclide("U234", 1.110770e-05)
mat2.add_nuclide("U236", 5.325150e-06)
mat2.add_element("Cr", 2.504980e-03)
mat2.add_element("Ni", 1.110580e-03)
mat2.add_element("Fe", 7.026370e-02)
mat2.add_element("Al", 1.089720e-03)
mat2.add_element("C", 3.281090e-04)
mat2.add_element("Mo", 9.421530e-06)
mat2.add_element("Mn", 5.362480e-04)
mat2.add_element("Cu", 1.907340e-05)
mat2.add_element("H", 1.001680e-06)
mat2.add_element("Si", 1.273230e-04)
mat2.add_element("Mg", 3.137870e-05)
mat2.add_element("Cl", 1.752420e-06)
mat2.add_element("F", 5.314080e-06)

mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("U235", 1.028540e-03)
mat3.add_nuclide("U238", 5.961660e-05)
mat3.add_nuclide("U234", 1.008800e-05)
mat3.add_nuclide("U236", 4.836090e-06)
mat3.add_element("Cr", 2.322850e-03)
mat3.add_element("Ni", 1.084000e-03)
mat3.add_element("Fe", 6.933730e-02)
mat3.add_element("Al", 1.178030e-03)
mat3.add_element("C", 3.341740e-04)
mat3.add_element("Mo", 1.493620e-05)
mat3.add_element("Mn", 4.928720e-04)
mat3.add_element("Cu", 2.567280e-05)
mat3.add_element("H", 1.988350e-06)
mat3.add_element("Si", 1.273920e-04)
mat3.add_element("Mg", 3.392250e-05)
mat3.add_element("Cl", 3.478580e-06)
mat3.add_element("F", 1.054850e-05)

mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_nuclide("U235", 1.038160e-03)
mat4.add_nuclide("U238", 6.017430e-05)
mat4.add_nuclide("U234", 1.018210e-05)
mat4.add_nuclide("U236", 4.881420e-06)
mat4.add_element("Cr", 2.359180e-03)
mat4.add_element("Ni", 1.103320e-03)
mat4.add_element("Fe", 6.962390e-02)
mat4.add_element("Al", 1.148850e-03)
mat4.add_element("C", 3.328370e-04)
mat4.add_element("Mo", 1.496640e-05)
mat4.add_element("Mn", 4.963390e-04)
mat4.add_element("Cu", 2.571450e-05)
mat4.add_element("H", 2.044630e-06)
mat4.add_element("Si", 1.283580e-04)
mat4.add_element("Mg", 3.308200e-05)
mat4.add_element("Cl", 3.487620e-06)
mat4.add_element("F", 1.057600e-05)

mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")
mat5.add_element("Cr", 1.498360e-02)
mat5.add_element("Ni", 6.610340e-03)
mat5.add_element("Fe", 5.323250e-02)
mat5.add_element("Al", 1.153770e-03)
mat5.add_element("C", 7.237040e-05)
mat5.add_element("Mo", 8.510610e-06)
mat5.add_element("Mn", 1.444660e-03)
mat5.add_element("Cu", 1.752320e-05)
mat5.add_element("Si", 1.022720e-03)
mat5.add_element("Mg", 3.322220e-05)

mat6 = openmc.Material(material_id=6)
mat6.set_density("sum")
mat6.add_element("Cr", 1.493340e-02)
mat6.add_element("Ni", 6.456530e-03)
mat6.add_element("Fe", 5.278970e-02)
mat6.add_element("Al", 2.352450e-03)
mat6.add_element("C", 8.216370e-05)
mat6.add_element("Mo", 4.986400e-05)
mat6.add_element("Mn", 1.326600e-03)
mat6.add_element("Cu", 4.261410e-05)
mat6.add_element("Si", 9.320720e-04)
mat6.add_element("Mg", 1.113000e-04)

mat7 = openmc.Material(material_id=7)
mat7.set_density("sum")
mat7.add_nuclide("U235", 8.594760e-05)
mat7.add_nuclide("U238", 3.849880e-02)
mat7.add_element("Cr", 1.860220e-03)
mat7.add_element("Ni", 8.089170e-04)
mat7.add_element("Fe", 6.866420e-03)
mat7.add_element("C", 2.960460e-05)
mat7.add_element("Mo", 8.190330e-06)
mat7.add_element("Mn", 1.605880e-04)
mat7.add_element("Cu", 1.709350e-05)
mat7.add_element("H", 2.424230e-06)
mat7.add_element("Si", 1.474200e-04)
mat7.add_element("Cl", 4.204060e-06)
mat7.add_element("F", 1.244950e-05)

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7])

# ==============================================================================
# Geometry
# ==============================================================================


# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# CORE1
cell1 = openmc.Cell(cell_id=1, fill=mat1)
cell1.region = -surf1 & +surf6 & -surf7

# CORE2
cell2 = openmc.Cell(cell_id=2, fill=mat2)
cell2.region = -surf1 & +surf5 & -surf6

# CORE3
cell3 = openmc.Cell(cell_id=3, fill=mat3)
cell3.region = -surf1 & +surf4 & -surf5

# CORE4
cell4 = openmc.Cell(cell_id=4, fill=mat4)
cell4.region = -surf1 & +surf7 & -surf8

# AXSST
cell5 = openmc.Cell(cell_id=5, fill=mat5)
cell5.region = -surf1 & -surf3 & +surf8

# AXSST
cell6 = openmc.Cell(cell_id=6, fill=mat5)
cell6.region = -surf1 & -surf3 & -surf4

# RDSST
cell7 = openmc.Cell(cell_id=7, fill=mat6)
cell7.region = +surf1 & -surf2 & -surf3

# DURRS
cell8 = openmc.Cell(cell_id=8, fill=mat7)
cell8.region = +surf2 & -surf3

root_universe = openmc.Universe(cells=[cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8])
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
