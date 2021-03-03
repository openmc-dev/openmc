.. _io_summary:

===================
Summary File Format
===================

The current version of the summary file format is 6.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the summary
               file format.
             - **openmc_version** (*int[3]*) -- Major, minor, and release
               version number for OpenMC.
             - **git_sha1** (*char[40]*) -- Git commit SHA-1 hash.
             - **date_and_time** (*char[]*) -- Date and time the summary was
               written.

**/geometry/**

:Attributes: - **n_cells** (*int*) -- Number of cells in the problem.
             - **n_surfaces** (*int*) -- Number of surfaces in the problem.
             - **n_universes** (*int*) -- Number of unique universes in the
               problem.
             - **n_lattices** (*int*) -- Number of lattices in the problem.
             - **dagmc** (*int*) -- Indicates that a DAGMC geometry was used
               if present.

**/geometry/cells/cell <uid>/**

:Datasets: - **name** (*char[]*) -- User-defined name of the cell.
           - **universe** (*int*) -- Universe assigned to the cell. If none is
             specified, the default universe (0) is assigned.
           - **fill_type** (*char[]*) -- Type of fill for the cell. Can be
             'material', 'universe', or 'lattice'.
           - **material** (*int* or *int[]*) -- Unique ID of the material(s)
             assigned to the cell. This dataset is present only if fill_type is
             set to 'normal'.  The value '-1' signifies void material.  The data
             is an array if the cell uses distributed materials, otherwise it is
             a scalar.
           - **temperature** (*double[]*) -- Temperature of the cell in Kelvin.
           - **translation** (*double[3]*) -- Translation applied to the fill
             universe. This dataset is present only if fill_type is set to
             'universe'.
           - **rotation** (*double[3]*) -- Angles in degrees about the x-, y-,
             and z-axes for which the fill universe should be rotated. This
             dataset is present only if fill_type is set to 'universe'.
           - **lattice** (*int*) -- Unique ID of the lattice which fills the
             cell. Only present if fill_type is set to 'lattice'.
           - **region** (*char[]*) -- Region specification for the cell.

**/geometry/surfaces/surface <uid>/**

:Datasets: - **name** (*char[]*) -- Name of the surface.
           - **type** (*char[]*) -- Type of the surface. Can be 'x-plane',
             'y-plane', 'z-plane', 'plane', 'x-cylinder', 'y-cylinder',
             'z-cylinder', 'sphere', 'x-cone', 'y-cone', 'z-cone', or 'quadric'.
           - **coefficients** (*double[]*) -- Array of coefficients that define
             the surface. See :ref:`surface_element` for what coefficients are
             defined for each surface type.
           - **boundary_condition** (*char[]*) -- Boundary condition applied to
             the surface. Can be 'transmission', 'vacuum', 'reflective', or
             'periodic'.

**/geometry/universes/universe <uid>/**

:Datasets:
           - **cells** (*int[]*) -- Array of unique IDs of cells that appear in
             the universe.

**/geometry/lattices/lattice <uid>/**

:Datasets: - **name** (*char[]*) -- Name of the lattice.
           - **type** (*char[]*) -- Type of the lattice, either 'rectangular' or
             'hexagonal'.
           - **pitch** (*double[]*) -- Pitch of the lattice in centimeters.
           - **outer** (*int*) -- Outer universe assigned to lattice cells
             outside the defined range.
           - **universes** (*int[][][]*) -- Three-dimensional array of universes
             assigned to each cell of the lattice.
           - **dimension** (*int[]*) -- The number of lattice cells in each
             direction. This dataset is present only when the 'type' dataset is
             set to 'rectangular'.
           - **lower_left** (*double[]*) -- The coordinates of the lower-left
             corner of the lattice. This dataset is present only when the 'type'
             dataset is set to 'rectangular'.
           - **n_rings** (*int*) -- Number of radial ring positions in the
             xy-plane. This dataset is present only when the 'type' dataset is
             set to 'hexagonal'.
           - **n_axial** (*int*) -- Number of lattice positions along the
             z-axis. This dataset is present only when the 'type' dataset is set
             to 'hexagonal'.
           - **center** (*double[]*) -- Coordinates of the center of the
             lattice. This dataset is present only when the 'type' dataset is
             set to 'hexagonal'.

**/materials/**

:Attributes: - **n_materials** (*int*) -- Number of materials in the problem.


**/materials/material <uid>/**

:Datasets: - **name** (*char[]*) -- Name of the material.
           - **atom_density** (*double[]*) -- Total atom density of the material
             in atom/b-cm.
           - **nuclides** (*char[][]*) -- Array of nuclides present in the
             material, e.g., 'U235'. This data set is only present if nuclides
             are used.
           - **nuclide_densities** (*double[]*) -- Atom density of each nuclide.
             This data set is only present if 'nuclides' data set is present.
           - **macroscopics** (*char[][]*) -- Array of macroscopic data sets
             present in the material. This dataset is only present if
             macroscopic data sets are used in multi-group mode.
           - **sab_names** (*char[][]*) -- Names of
             S(:math:`\alpha,\beta`) tables assigned to the material.

:Attributes: - **volume** (*double[]*) -- Volume of this material [cm^3]. Only
               present if ``volume`` supplied
             - **temperature** (*double[]*) -- Temperature of this material [K].
               Only present in ``temperature`` supplied
             - **depletable** (*int[]*) -- ``1`` if the material can be depleted,
               ``0`` otherwise. Always present

**/nuclides/**

:Attributes: - **n_nuclides** (*int*) -- Number of nuclides in the problem.

:Datasets: - **names** (*char[][]*) -- Names of nuclides.
           - **awrs** (*float[]*) -- Atomic weight ratio of each nuclide.

**/macroscopics/**

:Attributes:
             - **n_macroscopics** (*int*) -- Number of macroscopic data sets
               in the problem.

:Datasets: - **names** (*char[][]*) -- Names of the macroscopic data sets.

**/tallies/tally <uid>/**

:Datasets: - **name** (*char[]*) -- Name of the tally.
