.. _capi:

=====
C API
=====

.. c:function:: int openmc_find(double* xyz, int rtype)

   Return the ID of the cell/material containing a given point

   :param xyz: Cartesian coordinates
   :type xyz: double[3]
   :param rtype: Which ID to return (1=cell, 2=material)
   :type rtype: int
   :rtype: int

.. c:function:: void openmc_init(int intracomm)

   Initialize OpenMC

   :param intracomm: MPI intracommunicator
   :type intracomm: int

.. c:function:: void openmc_finalize()

   Finalize a simulation

.. c:function:: void openmc_reset()

   Resets all tally scores

.. c:function:: void openmc_run()

   Run a simulation
