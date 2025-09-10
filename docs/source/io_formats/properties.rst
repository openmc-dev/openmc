.. _io_properties:

======================
Properties File Format
======================

The current version of the properties file format is 1.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the
               statepoint file format.
             - **openmc_version** (*int[3]*) -- Major, minor, and release
               version number for OpenMC.
             - **git_sha1** (*char[40]*) -- Git commit SHA-1 hash.
             - **date_and_time** (*char[]*) -- Date and time the summary was
               written.
             - **path** (*char[]*) -- Path to directory containing input files.

**/geometry/**

:Attributes: - **n_cells** (*int*) -- Number of cells in the problem.

**/geometry/cells/cell <uid>/**

:Datasets: - **temperature** (*double[]*) -- Temperature of the cell in [K].
           - **density_mult** (*double[]*) -- Unitless density multipliers for
             the cell. The cell density is equal to the density multiplier
             times the density of the material filling the cell.

**/materials/**

:Attributes: - **n_materials** (*int*) -- Number of materials in the problem.

**/materials/material <uid>/**

:Attributes: - **atom_density** (*double*) -- Total density in [atom/b-cm].
             - **mass_density** (*double*) -- Total density in [g/cm^3].
