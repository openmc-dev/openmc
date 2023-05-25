.. _io_weight_windows:

====================
Weight Window Format
====================

The current revision of the weight window file format is 1.0.

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file.
             - **version** (*int[2]*) -- Major and minor version of the weight
               window file format.

**/weight_windows/**

:Attributes: - **n_weight_windows** (*int*) -- Number of weight window objects in the file.
             - **ids** (*int[]*) -- Unique IDs of weight window objects in the file.

**/weight_windows/weight_windows_<uid>/**

:Datasets: - **mesh** (*int*) -- ID of the mesh associated with the weight window object.
           - **particle_type** (*char[]*)  -- Particle type to which the weight windows apply.
           - **energy_bounds** (*double[]*) -- Energy bounds of the weight windows in [eV]
           - **lower_ww_bounds** (*double[]*) -- Weight window lower bounds.
           - **upper_ww_bounds** (*double[]*) -- Weight window upper bounds.
           - **survival_ratio** (*double*) -- Weight window survival ratio.
           - **max_lower_bound_ratio** (*double*) -- Maximum particle weight to lower weight window bound ratio.
           - **max_split** (*int*) -- Maximum number of splits per weight window check.
           - **weight_cutoff** (*double*) -- Particle weight cutoff.

**/meshes/**

:Attributes: - **n_meshes** (*int*) -- Number of meshes in the file.
             - **ids** (*int[]*) -- User-defined unique ID of each mesh.

**/meshes/mesh <uid>/**

Please see the section on **/tallies/meshes/** in the :doc:`statepoint`.
