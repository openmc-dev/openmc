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

**/weight_windows/weight_windows <uid>/**

:Datasets: - **mesh** (*int*) -- ID of the mesh associated with the weight window object.
           - **particle_type** (*char[]*)  -- Particle type to which the weight windows apply.
           - **energy_bounds** (*double[]*) -- Energy bounds of the weight windows
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

:Datasets: - **type** (*char[]*) -- Type of mesh.
           - **dimension** (*int*) -- Number of mesh cells in each dimension.
           - **Regular Mesh Only:**
              - **lower_left** (*double[]*) -- Coordinates of lower-left corner of
                mesh.
              - **upper_right** (*double[]*) -- Coordinates of upper-right corner
                of mesh.
              - **width** (*double[]*) -- Width of each mesh cell in each
                dimension.
           - **Rectilinear Mesh Only:**
              - **x_grid** (*double[]*) -- Mesh divisions along the x-axis.
              - **y_grid** (*double[]*) -- Mesh divisions along the y-axis.
              - **z_grid** (*double[]*) -- Mesh divisions along the z-axis.
           - **Cylindrical & Spherical Mesh Only:**
              - **r_grid** (*double[]*) -- The mesh divisions along the r-axis.
              - **phi_grid** (*double[]*) -- The mesh divisions along the phi-axis.
              - **origin** (*double[]*) -- The origin in cartesian coordinates.
           - **Spherical Mesh Only:**
              - **theta_grid** (*double[]*) -- The mesh divisions along the theta-axis.
           - **Unstructured Mesh Only:**
              - **filename** (*char[]*) -- Name of the mesh file.
              - **library** (*char[]*) -- Mesh library used to represent the
                                          mesh ("moab" or "libmesh").
              - **length_multiplier** (*double*) Scaling factor applied to the mesh.
              - **volumes** (*double[]*) -- Volume of each mesh cell.
              - **vertices** (*double[]*) -- x, y, z values of the mesh vertices.
              - **connectivity** (*int[]*) -- Connectivity array for the mesh
                cells.
              - **element_types** (*int[]*) -- Mesh element types.

