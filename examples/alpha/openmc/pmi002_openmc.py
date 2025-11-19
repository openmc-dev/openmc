"""
PU-MET-INTER-002: ZPR-6/10 Loading 24; Benchmark Model
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1, name="")
mat1.set_density("atom/b-cm", 1.114996e-03)
mat1.add_nuclide("Pu238", 1.69380e-9)
mat1.add_nuclide("Pu239", 1.06299e-3)
mat1.add_nuclide("Pu240", 5.03037e-5)
mat1.add_nuclide("Pu241", 1.70034e-6)

mat2 = openmc.Material(material_id=2, name="")
mat2.set_density("atom/b-cm", 2.423882e-02)
mat2.add_element("Cr", 1.57277e-2)
mat2.add_element("Mn", 1.56081e-3)
mat2.add_element("Co", 3.73582e-6)
mat2.add_element("Ni", 6.94657e-3)

mat3 = openmc.Material(material_id=3, name="")
mat3.set_density("atom/b-cm", 2.383705e-02)
mat3.add_element("Cr", 1.56286e-2)
mat3.add_element("Mn", 1.42057e-3)
mat3.add_element("Co", 3.77932e-6)
mat3.add_element("Ni", 6.78410e-3)

mat4 = openmc.Material(material_id=4, name="")
mat4.set_density("atom/b-cm", 3.143137e-03)
mat4.add_element("Cr", 1.82758e-3)
mat4.add_element("Mn", 5.50350e-4)
mat4.add_element("Co", 3.77222e-6)
mat4.add_element("Ni", 7.61435e-4)

mat5 = openmc.Material(material_id=5, name="")
mat5.set_density("atom/b-cm", 1.791216e-03)
mat5.add_element("Cr", 1.18727e-3)
mat5.add_element("Mn", 1.05699e-4)
mat5.add_element("Ni", 4.79539e-4)
mat5.add_element("C", 1.87078e-5)
mat5.add_s_alpha_beta("c_Graphite")

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================


# Cell: CORE
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="CORE")
cell0.region = +surf1 & -surf2 & -surf5

# Cell: SSAXR
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="SSAXR")
cell1.region = +surf2 & -surf3 & -surf5

# Cell: SSRDR
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="SSRDR")
cell2.region = +surf1 & -surf3 & +surf5 & -surf6

# Cell: FERDR
cell3 = openmc.Cell(cell_id=3, fill=mat4, name="FERDR")
cell3.region = +surf1 & -surf3 & +surf6 & -surf7

# Cell: MATRX
cell4 = openmc.Cell(cell_id=4, fill=mat5, name="MATRX")
cell4.region = +surf1 & -surf3 & +surf7 & -surf8

# Cell: MATRX
cell5 = openmc.Cell(cell_id=5, fill=mat5, name="MATRX")
cell5.region = +surf3 & -surf4 & -surf8

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
outer_cell = openmc.Cell(cell_id=6, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, outer_cell])
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
source.space = openmc.stats.Point((0.0001, 0.0001, 0.0001))
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
