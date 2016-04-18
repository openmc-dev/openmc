.. _usersguide_summary:

===================
Summary File Format
===================

The current revision of the summary file format is 1.

**/filetype** (*char[]*)

    String indicating the type of file.

**/revision** (*int*)

    Revision of the summary file format. Any time a change is made in the
    format, this integer is incremented.

**/version_major** (*int*)

    Major version number for OpenMC

**/version_minor** (*int*)

    Minor version number for OpenMC

**/version_release** (*int*)

    Release version number for OpenMC

**/date_and_time** (*char[]*)

    Date and time the summary was written.

**/n_procs** (*int*)

    Number of MPI processes used.

**/n_particles** (*int8_t*)

    Number of particles used per generation.

**/n_batches** (*int*)

    Number of batches to simulate.

**/n_inactive** (*int*)

    Number of inactive batches. Only present if /run_mode is set to
    'k-eigenvalue'.

**/n_active** (*int*)

    Number of active batches. Only present if /run_mode is set to
    'k-eigenvalue'.

**/gen_per_batch** (*int*)

    Number of generations per batch. Only present if /run_mode is set to
    'k-eigenvalue'.

**/geometry/n_cells** (*int*)

    Number of cells in the problem.

**/geometry/n_surfaces** (*int*)

    Number of surfaces in the problem.

**/geometry/n_universes** (*int*)

    Number of unique universes in the problem.

**/geometry/n_lattices** (*int*)

    Number of lattices in the problem.

**/geometry/cells/cell <uid>/index** (*int*)

    Index in cells array used internally in OpenMC.

**/geometry/cells/cell <uid>/name** (*char[]*)

    Name of the cell.

**/geometry/cells/cell <uid>/universe** (*int*)

    Universe assigned to the cell. If none is specified, the default
    universe (0) is assigned.

**/geometry/cells/cell <uid>/fill_type** (*char[]*)

    Type of fill for the cell. Can be 'normal', 'universe', or 'lattice'.

**/geometry/cells/cell <uid>/material** (*int* or *int[]*)

    Unique ID of the material(s) assigned to the cell. This dataset is present
    only if fill_type is set to 'normal'.  The value '-1' signifies void
    material.  The data is an array if the cell uses distributed materials,
    otherwise it is a scalar.

**/geometry/cells/cell <uid>/offset** (*int[]*)

    Offsets used for distribcell tally filter. This dataset is present only if
    fill_type is set to 'universe'.

**/geometry/cells/cell <uid>/translation** (*double[3]*)

    Translation applied to the fill universe. This dataset is present only if
    fill_type is set to 'universe'.

**/geometry/cells/cell <uid>/rotation** (*double[3]*)

    Angles in degrees about the x-, y-, and z-axes for which the fill universe
    should be rotated. This dataset is present only if fill_type is set to
    'universe'.

**/geometry/cells/cell <uid>/lattice** (*int*)

    Unique ID of the lattice which fills the cell. Only present if fill_type is
    set to 'lattice'.

**/geometry/cells/cell <uid>/region** (*char[]*)

    Region specification for the cell.

**/geometry/cells/cell <uid>/distribcell_index** (*int*)

    Index of this cell in distribcell filter arrays.

**/geometry/surfaces/surface <uid>/index** (*int*)

    Index in surfaces array used internally in OpenMC.

**/geometry/surfaces/surface <uid>/name** (*char[]*)

    Name of the surface.

**/geometry/surfaces/surface <uid>/type** (*char[]*)

    Type of the surface. Can be 'x-plane', 'y-plane', 'z-plane', 'plane',
    'x-cylinder', 'y-cylinder', 'sphere', 'x-cone', 'y-cone', 'z-cone', or
    'quadric'.

**/geometry/surfaces/surface <uid>/coefficients** (*double[]*)

    Array of coefficients that define the surface. See :ref:`surface_element`
    for what coefficients are defined for each surface type.

**/geometry/surfaces/surface <uid>/boundary_condition** (*char[]*)

    Boundary condition applied to the surface. Can be 'transmission', 'vacuum',
    'reflective', or 'periodic'.

**/geometry/universes/universe <uid>/index** (*int*)

    Index in the universes array used internally in OpenMC.

**/geometry/universes/universe <uid>/cells** (*int[]*)

    Array of unique IDs of cells that appear in the universe.

**/geometry/lattices/lattice <uid>/index** (*int*)

    Index in the lattices array used internally in OpenMC.

**/geometry/lattices/lattice <uid>/name** (*char[]*)

    Name of the lattice.

**/geometry/lattices/lattice <uid>/type** (*char[]*)

    Type of the lattice, either 'rectangular' or 'hexagonal'.

**/geometry/lattices/lattice <uid>/pitch** (*double[]*)

    Pitch of the lattice.

**/geometry/lattices/lattice <uid>/outer** (*int*)

    Outer universe assigned to lattice cells outside the defined range.

**/geometry/lattices/lattice <uid>/offsets** (*int[]*)

    Offsets used for distribcell tally filter.

**/geometry/lattices/lattice <uid>/universes** (*int[]*)

    Three-dimensional array of universes assigned to each cell of the lattice.

**/geometry/lattices/lattice <uid>/dimension** (*int[]*)

    The number of lattice cells in each direction. This dataset is present only
    when the 'type' dataset is set to 'rectangular'.

**/geometry/lattices/lattice <uid>/lower_left** (*double[]*)

    The coordinates of the lower-left corner of the lattice. This dataset is
    present only when the 'type' dataset is set to 'rectangular'.

**/geometry/lattices/lattice <uid>/n_rings** (*int*)

    Number of radial ring positions in the xy-plane. This dataset is present
    only when the 'type' dataset is set to 'hexagonal'.

**/geometry/lattices/lattice <uid>/n_axial** (*int*)

    Number of lattice positions along the z-axis. This dataset is present only
    when the 'type' dataset is set to 'hexagonal'.

**/geometry/lattices/lattice <uid>/center** (*double[]*)

    Coordinates of the center of the lattice. This dataset is present only when
    the 'type' dataset is set to 'hexagonal'.

**/n_materials** (*int*)

    Number of materials in the problem.

**/materials/material <uid>/index** (*int*)

    Index in materials array used internally in OpenMC.

**/materials/material <uid>/name** (*char[]*)

    Name of the material.

**/materials/material <uid>/atom_density** (*double[]*)

    Total atom density of the material in atom/b-cm.

**/materials/material <uid>/nuclides** (*char[][]*)

    Array of nuclides present in the material, e.g., 'U-235.71c'.

**/materials/material <uid>/nuclide_densities** (*double[]*)

    Atom density of each nuclide.

**/materials/material <uid>/sab_names** (*char[][]*)

    Names of S(:math:`\alpha`,:math:`\beta`) tables assigned to the material.

**/tallies/n_tallies** (*int*)

    Number of tallies in the problem.

**/tallies/n_meshes** (*int*)

    Number of meshes in the problem.

**/tallies/mesh <uid>/index** (*int*)

    Index in the meshes array used internally in OpenMC.

**/tallies/mesh <uid>/type** (*char[]*)

    Type of the mesh. The only valid option is currently 'regular'.

**/tallies/mesh <uid>/dimension** (*int[]*)

    Number of mesh cells in each direction.

**/tallies/mesh <uid>/lower_left** (*double[]*)

    Coordinates of the lower-left corner of the mesh.

**/tallies/mesh <uid>/upper_right** (*double[]*)

    Coordinates of the upper-right corner of the mesh.

**/tallies/mesh <uid>/width** (*double[]*)

    Width of a single mesh cell in each direction.

**/tallies/tally <uid>/index** (*int*)

    Index in tallies array used internally in OpenMC.

**/tallies/tally <uid>/name** (*char[]*)

    Name of the tally.

**/tallies/tally <uid>/n_filters** (*int*)

    Number of filters applied to the tally.

**/tallies/tally <uid>/filter <j>/type** (*char[]*)

    Type of the j-th filter. Can be 'universe', 'material', 'cell', 'cellborn',
    'surface', 'mesh', 'energy', 'energyout', or 'distribcell'.

**/tallies/tally <uid>/filter <j>/offset** (*int*)

    Filter offset (used for distribcell filter).

**/tallies/tally <uid>/filter <j>/paths** (*char[][]*)

    The paths traversed through the CSG tree to reach each distribcell
    instance (for 'distribcell' filters only). This consists of the integer
    IDs for each universe, cell and lattice delimited by '->'. Each lattice
    cell is specified by its (x,y) or (x,y,z) indices.

**/tallies/tally <uid>/filter <j>/n_bins** (*int*)

    Number of bins for the j-th filter.

**/tallies/tally <uid>/filter <j>/bins** (*int[]* or *double[]*)

    Value for each filter bin of this type.

**/tallies/tally <uid>/nuclides** (*char[][]*)

    Array of nuclides to tally. Note that if no nuclide is specified in the user
    input, a single 'total' nuclide appears here.

**/tallies/tally <uid>/n_score_bins** (*int*)

    Number of scoring bins for a single nuclide. In general, this can be greater
    than the number of user-specified scores since each score might have
    multiple scoring bins, e.g., scatter-PN.

**/tallies/tally <uid>/moment_orders** (*char[][]*)

    Tallying moment orders for Legendre and spherical harmonic tally expansions
    (*e.g.*, 'P2', 'Y1,2', etc.).

**/tallies/tally <uid>/score_bins** (*char[][]*)

    Scoring bins for the tally.
