.. _io_depletion_results:

=============================
Depletion Results File Format
=============================

The current version of the depletion results file format is 1.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the
               statepoint file format.

:Datasets: - **eigenvalues** (*double[][][2]*) -- k-eigenvalues at each
             time/stage. This array has shape (number of timesteps, number of
             stages, value). The last axis contains the eigenvalue and the 
             associated uncertainty
           - **number** (*double[][][][]*) -- Total number of atoms. This array
             has shape (number of timesteps, number of stages, number of
             materials, number of nuclides).
           - **reaction rates** (*double[][][][][]*) -- Reaction rates used to
             build depletion matrices. This array has shape (number of
             timesteps, number of stages, number of materials, number of
             nuclides, number of reactions).
           - **time** (*double[][2]*) -- Time in [s] at beginning/end of each
             step.
           - **depletion time** (*double[]*) -- Average process time in [s] 
             spent depleting a material across all burnable materials and,
             if applicable, MPI processes.

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

.. note::

    The reaction rates for some isotopes not originally present may
    be non-zero, but should be negligible compared to other atoms.
    This can be controlled by changing the
    :class:`openmc.deplete.Operator` ``dilute_initial`` attribute.
