.. _capi:

=====
C API
=====

.. c:function:: void openmc_calculate_volumes()

   Run a stochastic volume calculation

.. c:function:: int openmc_cell_get_fill(int32_t index, int* type, int32_t** indices, int32_t* n)

   Get the fill for a cell

   :param int32_t index: Index in the cells array
   :param int* type: Type of the fill
   :param int32_t** indices: Array of material indices for cell
   :param int32_t n: Length of indices array
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_cell_get_id(int32_t index, int32_t* id)

   Get the ID of a cell

   :param int32_t index: Index in the cells array
   :param int32_t* id: ID of the cell
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_cell_set_fill(int32_t index, int type, int32_t n, int32_t* indices)

   Set the fill for a cell

   :param int32_t index: Index in the cells array
   :param int type: Type of the fill
   :param int32_t n: Length of indices array
   :param int32_t* indices: Array of material indices for cell
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_cell_set_id(int32_t index, int32_t id)

   Set the ID of a cell

   :param int32_t index: Index in the cells array
   :param int32_t id: ID of the cell
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_cell_set_temperature(index index, double T, int32_t* instance)

   Set the temperature of a cell.

   :param int32_t index: Index in the cells array
   :param double T: Temperature in Kelvin
   :param int32_t* instance: Which instance of the cell. To set the temperature
                             for all instances, pass a null pointer.
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_energy_filter_get_bins(int32_t index, double** energies, int32_t* n)

   Return the bounding energies for an energy filter

   :param int32_t index: Index in the filters array
   :param double** energies: Bounding energies of the bins for the energy filter
   :param int32_t* n: Number of energies specified
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_energy_filter_set_bins(int32_t index, int32_t n, double* energies)

   Set the bounding energies for an energy filter

   :param int32_t index: Index in the filters array
   :param int32_t n: Number of energies specified
   :param double* energies: Bounding energies of the bins for the energy filter
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_extend_cells(int32_t n, int32_t* index_start, int32_t* index_end)

   Extend the cells array by n elements

   :param int32_t n: Number of cells to create
   :param int32_t* index_start: Index of first new cell
   :param int32_t* index_end: Index of last new cell
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_extend_filters(int32_t n, int32_t* index_start, int32_t* index_end)

   Extend the filters array by n elements

   :param int32_t n: Number of filters to create
   :param int32_t* index_start: Index of first new filter
   :param int32_t* index_end: Index of last new filter
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_extend_materials(int32_t n, int32_t* index_start, int32_t* index_end)

   Extend the materials array by n elements

   :param int32_t n: Number of materials to create
   :param int32_t* index_start: Index of first new material
   :param int32_t* index_end: Index of last new material
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_extend_tallies(int32_t n, int32_t* index_start, int32_t* index_end)

   Extend the tallies array by n elements

   :param int32_t n: Number of tallies to create
   :param int32_t* index_start: Index of first new tally
   :param int32_t* index_end: Index of last new tally
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_filter_get_id(int32_t index, int32_t* id)

   Get the ID of a filter

   :param int32_t index: Index in the filters array
   :param int32_t* id: ID of the filter
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_filter_set_id(int32_t index, int32_t id)

   Set the ID of a filter

   :param int32_t index: Index in the filters array
   :param int32_t id: ID of the filter
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: void openmc_finalize()

   Finalize a simulation

.. c:function:: void openmc_find(double* xyz, int rtype, int32_t* id, int32_t* instance)

   Determine the ID of the cell/material containing a given point

   :param double[3] xyz: Cartesian coordinates
   :param int rtype: Which ID to return (1=cell, 2=material)
   :param int32_t* id: ID of the cell/material found. If a material is requested
                       and the point is in a void, the ID is 0. If an error
                       occurs, the ID is -1.
   :param int32_t* instance: If a cell is repetaed in the geometry, the instance
                             of the cell that was found and zero otherwise.

.. c:function:: int openmc_get_cell_index(int32_t id, int32_t* index)

   Get the index in the cells array for a cell with a given ID

   :param int32_t id: ID of the cell
   :param int32_t* index: Index in the cells array
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_filter_index(int32_t id, int32_t* index)

   Get the index in the filters array for a filter with a given ID

   :param int32_t id: ID of the filter
   :param int32_t* index: Index in the filters array
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: void openmc_get_filter_next_id(int32_t* id)

   Get an integer ID that has not been used by any filters.

   :param int32_t* id: Unused integer ID

.. c:function:: int openmc_get_keff(double k_combined[])

   :param double[2] k_combined: Combined estimate of k-effective
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_material_index(int32_t id, int32_t* index)

   Get the index in the materials array for a material with a given ID

   :param int32_t id: ID of the material
   :param int32_t* index: Index in the materials array
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_nuclide_index(char name[], int* index)

   Get the index in the nuclides array for a nuclide with a given name

   :param char[] name: Name of the nuclide
   :param int* index: Index in the nuclides array
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_tally_index(int32_t id, int32_t* index)

   Get the index in the tallies array for a tally with a given ID

   :param int32_t id: ID of the tally
   :param int32_t* index: Index in the tallies array
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: void openmc_hard_reset()

   Reset tallies, timers, and pseudo-random number generator state

.. c:function:: void openmc_init(int intracomm)

   Initialize OpenMC

   :param int intracomm: MPI intracommunicator

.. c:function:: int openmc_load_nuclide(char name[])

   Load data for a nuclide from the HDF5 data library.

   :param char[] name: Name of the nuclide.
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_add_nuclide(int32_t index, char name[], double density)

   Add a nuclide to an existing material. If the nuclide already exists, the
   density is overwritten.

   :param int32_t index: Index in the materials array
   :param char[] name: Name of the nuclide
   :param double density: Density in atom/b-cm
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_get_densities(int32_t index, int* nuclides[], double* densities[], int* n)

   Get density for each nuclide in a material.

   :param int32_t index: Index in the materials array
   :param int** nuclides: Pointer to array of nuclide indices
   :param double** densities: Pointer to the array of densities
   :param int* n: Length of the array
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_get_id(int32_t index, int32_t* id)

   Get the ID of a material

   :param int32_t index: Index in the materials array
   :param int32_t* id: ID of the material
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_material_set_density(int32_t index, double density)

   Set the density of a material.

   :param int32_t index: Index in the materials array
   :param double density: Density of the material in atom/b-cm
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_set_densities(int32_t index, int n, char* name[], double density[])

   :param int32_t index: Index in the materials array
   :param int n: Length of name/density
   :param char** name: Array of nuclide names
   :param double[] density: Array of densities
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_set_id(int32_t index, int32_t id)

   Set the ID of a material

   :param int32_t index: Index in the materials array
   :param int32_t id: ID of the material
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_material_filter_get_bins(int32_t index, int32_t** bins, int32_t* n)

   Get the bins for a material filter

   :param int32_t index: Index in the filters array
   :param int32_t** bins: Index in the materials array for each bin
   :param int32_t* n: Number of bins
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_material_filter_set_bins(int32_t index, int32_t n, int32_t* bins)

   Set the bins for a material filter

   :param int32_t index: Index in the filters array
   :param int32_t n: Number of bins
   :param int32_t* bins: Index in the materials array for each bin
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_mesh_filter_set_mesh(int32_t index, int32_t index_mesh)

   Set the mesh for a mesh filter

   :param int32_t index: Index in the filters array
   :param int32_t index_mesh: Index in the meshes array
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_next_batch()

   Simulate next batch of particles. Must be called after openmc_simulation_init().

   :return: Integer indicating whether simulation has finished (negative) or not
            finished (zero).
   :rtype: int

.. c:function:: int openmc_nuclide_name(int index, char* name[])

   Get name of a nuclide

   :param int index: Index in the nuclides array
   :param char** name: Name of the nuclide
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: void openmc_plot_geometry()

   Run plotting mode.

.. c:function:: void openmc_reset()

   Resets all tally scores

.. c:function:: void openmc_run()

   Run a simulation

.. c:function:: void openmc_simulation_finalize()

   Finalize a simulation.

.. c:function:: void openmc_simulation_init()

   Initialize a simulation. Must be called after openmc_init().

.. c:function:: int openmc_tally_get_id(int32_t index, int32_t* id)

   Get the ID of a tally

   :param int32_t index: Index in the tallies array
   :param int32_t* id: ID of the tally
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_get_filters(int32_t index, int32_t** indices, int* n)

   Get filters specified in a tally

   :param int32_t index: Index in the tallies array
   :param int32_t** indices: Array of filter indices
   :param int* n: Number of filters
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_get_n_realizations(int32_t index, int32_t* n)

   :param int32_t index: Index in the tallies array
   :param int32_t* n: Number of realizations
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_get_nuclides(int32_t index, int* nuclides[], int* n)

   Get nuclides specified in a tally

   :param int32_t index: Index in the tallies array
   :param int** nuclides: Array of nuclide indices
   :param int* n: Number of nuclides
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_get_scores(int32_t index, int** scores, int* n)

   Get scores specified for a tally

   :param int32_t index: Index in the tallies array
   :param int** scores: Array of scores
   :param int* n: Number of scores
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_results(int32_t index, double** ptr, int shape_[3])

   Get a pointer to tally results array.

   :param int32_t index: Index in the tallies array
   :param double** ptr: Pointer to the results array
   :param int[3] shape_: Shape of the results array
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_set_filters(int32_t index, int n, int32_t* indices)

   Set filters for a tally

   :param int32_t index: Index in the tallies array
   :param int n: Number of filters
   :param int32_t* indices: Array of filter indices
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_set_id(int32_t index, int32_t id)

   Set the ID of a tally

   :param int32_t index: Index in the tallies array
   :param int32_t id: ID of the tally
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_set_nuclides(int32_t index, int n, char* nuclides[])

   Set the nuclides for a tally

   :param int32_t index: Index in the tallies array
   :param int n: Number of nuclides
   :param char** nuclides: Array of nuclide names
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_set_scores(int32_t index, int n, int* scores)

   Set scores for a tally

   :param int32_t index: Index in the tallies array
   :param int n: Number of scores
   :param int* scores: Array of scores
   :return: Return status (negative if an error occurred)
   :rtype: int
