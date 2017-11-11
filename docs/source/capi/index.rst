.. _capi:

=====
C API
=====

.. c:function:: void openmc_calculate_volumes()

   Run a stochastic volume calculation

.. c:function:: int openmc_cell_get_id(int32_t index, int32_t* id)

   Get the ID of a cell

   :param index: Index in the cells array
   :type index: int32_t
   :param id: ID of the cell
   :type id: int32_t*
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_cell_set_temperature(index index, double T, int32_t* instance)

   Set the temperature of a cell.

   :param index: Index in the cells array
   :type index: int32_t
   :param T: Temperature in Kelvin
   :type T: double
   :param instance: Which instance of the cell. To set the temperature for all
                    instances, pass a null pointer.
   :type instance: int32_t*
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: void openmc_finalize()

   Finalize a simulation

.. c:function:: void openmc_find(double* xyz, int rtype, int32_t* id, int32_t* instance)

   Determine the ID of the cell/material containing a given point

   :param xyz: Cartesian coordinates
   :type xyz: double[3]
   :param rtype: Which ID to return (1=cell, 2=material)
   :type rtype: int
   :param id: ID of the cell/material found. If a material is requested and the
              point is in a void, the ID is 0. If an error occurs, the ID is -1.
   :type id: int32_t*
   :param instance: If a cell is repetaed in the geometry, the instance of the
                    cell that was found and zero otherwise.
   :type instance: int32_t*

.. c:function:: int openmc_get_cell_index(int32_t id, int32_t* index)

   Get the index in the cells array for a cell with a given ID

   :param id: ID of the cell
   :type id: int32_t
   :param index: Index in the cells array
   :type index: int32_t*
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_keff(double k_combined[])

   :param k_combined: Combined estimate of k-effective
   :type k_combined: double[2]
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_nuclide_index(char name[], int* index)

   Get the index in the nuclides array for a nuclide with a given name

   :param name: Name of the nuclide
   :type name: char[]
   :param index: Index in the nuclides array
   :type index: int*
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_tally_index(int32_t id, int32_t* index)

   Get the index in the tallies array for a tally with a given ID

   :param id: ID of the tally
   :type id: int32_t
   :param index: Index in the tallies array
   :type index: int32_t*
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_get_material_index(int32_t id, int32_t* index)

   Get the index in the materials array for a material with a given ID

   :param id: ID of the material
   :type id: int32_t
   :param index: Index in the materials array
   :type index: int32_t*
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: void openmc_hard_reset()

   Reset tallies, timers, and pseudo-random number generator state

.. c:function:: void openmc_init(int intracomm)

   Initialize OpenMC

   :param intracomm: MPI intracommunicator
   :type intracomm: int

.. c:function:: int openmc_load_nuclide(char name[])

   Load data for a nuclide from the HDF5 data library.

   :param name: Name of the nuclide.
   :type name: char[]
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_add_nuclide(int32_t index, char name[], double density)

   Add a nuclide to an existing material. If the nuclide already exists, the
   density is overwritten.

   :param index: Index in the materials array
   :type index: int32_t
   :param name: Name of the nuclide
   :type name: char[]
   :param density: Density in atom/b-cm
   :type density: double
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_get_densities(int32_t index, int* nuclides[], double* densities[])

   Get density for each nuclide in a material.

   :param index: Index in the materials array
   :type index: int32_t
   :param nuclides: Pointer to array of nuclide indices
   :type nuclides: int**
   :param densities: Pointer to the array of densities
   :type densities: double**
   :param n: Length of the array
   :type n: int
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_get_id(int32_t index, int32_t* id)

   Get the ID of a material

   :param index: Index in the materials array
   :type index: int32_t
   :param id: ID of the material
   :type id: int32_t*
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_material_set_density(int32_t index, double density)

   Set the density of a material.

   :param index: Index in the materials array
   :type index: int32_t
   :param density: Density of the material in atom/b-cm
   :type density: double
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_set_densities(int32_t, n, char* name[], double density[])

   :param index: Index in the materials array
   :type index: int32_t
   :param n: Length of name/density
   :type n: int
   :param name: Array of nuclide names
   :type name: char**
   :param density: Array of densities
   :type density: double[]
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_nuclide_name(int index, char* name[])

   Get name of a nuclide

   :param index: Index in the nuclides array
   :type index: int
   :param name: Name of the nuclide
   :type name: char**
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: void openmc_plot_geometry()

   Run plotting mode.

.. c:function:: void openmc_reset()

   Resets all tally scores

.. c:function:: void openmc_run()

   Run a simulation

.. c:function:: int openmc_tally_get_id(int32_t index, int32_t* id)

   Get the ID of a tally

   :param index: Index in the tallies array
   :type index: int32_t
   :param id: ID of the tally
   :type id: int32_t*
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_get_nuclides(int32_t index, int* nuclides[], int* n)

   Get nuclides specified in a tally

   :param index: Index in the tallies array
   :type index: int32_t
   :param nuclides: Array of nuclide indices
   :type nuclides: int**
   :param n: Number of nuclides
   :type n: int*
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_results(int32_t index, double** ptr, int shape_[3])

   Get a pointer to tally results array.

   :param index: Index in the tallies array
   :type index: int32_t
   :param ptr: Pointer to the results array
   :type ptr: double**
   :param shape_: Shape of the results array
   :type shape_: int[3]
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: int openmc_tally_set_nuclides(int32_t index, int n, char* nuclides[])

   Set the nuclides for a tally

   :param index: Index in the tallies array
   :type index: int32_t
   :param n: Number of nuclides
   :type n: int
   :param nuclides: Array of nuclide names
   :type nuclides: char**
   :return: Return status (negative if an error occurred)
   :rtype: int
