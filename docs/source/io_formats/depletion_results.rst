.. _io_depletion_results:

=============================
Depletion Results File Format
=============================

The current version of the depletion results file format is 1.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the
               statepoint file format.

:Datasets: - **eigenvalues** (*double[][]*) -- k-eigenvalues at each
             time/stage. This array has shape (number of timesteps, number of
             stages).
           - **number** (*double[][][][]*) -- Total number of atoms. This array
             has shape (number of timesteps, number of stages, number of
             materials, number of nuclides).
           - **reaction rates** (*double[][][][][]*) -- Reaction rates used to
             build depletion matrices. This array has shape (number of
             timesteps, number of stages, number of materials, number of
             nuclides, number of reactions).
           - **time** (*double[][2]*) -- Time in [s] at beginning/end of each
             step.

**/materials/<id>/**

:Attributes: - **index** (*int*) -- Index used in results for this material
             - **volume** (*double*) -- Volume of this material in [cm^3]

**/nuclides/<name>/**

:Attributes: - **atom number index** (*int*) -- Index in array of total atoms
               for this nuclide
             - **reaction rate index** (*int*) -- Index in array of reaction
               rates for this nuclide

**/reactions/<name>/**

:Attributes: - **index** (*int*) -- Index user in results for this reaction
