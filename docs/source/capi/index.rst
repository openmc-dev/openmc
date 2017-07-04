.. _capi:

=====
C API
=====

.. c:function:: void openmc_calculate_voumes()

   Run a stochastic volume calculation

.. c:function:: int openmc_cell_set_temperature(int id, double T, int* instance)

   Set the temperature of a cell.

   :param id: ID of the cell
   :type id: int
   :param T: Temperature in Kelvin
   :type T: double
   :param instance: Which instance of the cell. To set the temperature for all
                    instances, pass a null pointer.
   :type instance: int*
   :return: Return status (negative if an error occurred)
   :rtype: int

.. c:function:: void openmc_finalize()

   Finalize a simulation

.. c:function:: void openmc_find(double* xyz, int rtype, int* id, int* instance)

   Determine the ID of the cell/material containing a given point

   :param xyz: Cartesian coordinates
   :type xyz: double[3]
   :param rtype: Which ID to return (1=cell, 2=material)
   :type rtype: int
   :param id: ID of the cell/material found. If a material is requested and the
              point is in a void, the ID is 0. If an error occurs, the ID is -1.
   :type id: int
   :param instance: If a cell is repetaed in the geometry, the instance of the
                    cell that was found and zero otherwise.
   :type instance: int

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

.. c:function:: int openmc_material_add_nuclide(int id, char name[], double density)

   Add a nuclide to an existing material. If the nuclide already exists, the
   density is overwritten.

   :param id: ID of the material
   :type id: int
   :param name: Name of the nuclide
   :type name: char[]
   :param density: Density in atom/b-cm
   :type density: double
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_material_get_densities(int id, double* ptr)

   Get an array of nuclide densities for a material.

   :param id: ID of the material
   :type id: int
   :param ptr: Pointer to the array of densities
   :type ptr: double*
   :return: Length of the array
   :rtype: int

.. c:function:: int openmc_material_set_density(int id, double density)

   Set the density of a material.

   :param id: ID of the material
   :type id: int
   :param density: Density of the material in atom/b-cm
   :type density: double
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: void openmc_plot_geometry()

   Run plotting mode.

.. c:function:: void openmc_reset()

   Resets all tally scores

.. c:function:: void openmc_run()

   Run a simulation

.. c:function:: int openmc_set_density(double xyz[3], double density)

   Set the density of a material at a given point.

   :param xyz: Cartesian coordinates
   :type xyz: double[3]
   :param density: Density of the material to set in atom/b-cm
   :type density: double
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: int openmc_set_temperature(double xyz[3], double T)

   Set the density of a cell at a given point.

   :param xyz: Cartesian coordinates
   :type xyz: double[3]
   :param T: Temperature of the cell to set in Kelvin
   :type T: double
   :return: Return status (negative if an error occurs)
   :rtype: int

.. c:function:: void openmc_tally_results(int id, double** ptr, int shape_[3])

   Get a pointer to tally results array.

   :param id: ID of the tally
   :type id: int
   :param ptr: Pointer to the results array
   :type ptr: double**
   :param shape_: Shape of the results array
   :type shape_: int[3]
