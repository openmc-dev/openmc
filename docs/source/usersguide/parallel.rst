.. _usersguide_parallel:

===================
Running in Parallel
===================

If you are running a simulation on a computer with multiple cores, multiple
sockets, or multiple nodes (i.e., a cluster), you can benefit from the fact that
OpenMC is able to use all available hardware resources if configured
correctly. OpenMC is capable of using both distributed-memory (`MPI
<http://mpi-forum.org/>`_) and shared-memory (`OpenMP
<http://www.openmp.org/>`_) parallelism. If you are on a single-socket
workstation or a laptop, using shared-memory parallelism is likely
sufficient. On a multi-socket node, cluster, or supercomputer, chances are you
will need to use both distributed-memory (across nodes) and shared-memory
(within a single node) parallelism.

----------------------------------
Shared-Memory Parallelism (OpenMP)
----------------------------------

When using OpenMP, multiple threads will be launched and each is capable of
simulating a particle independently of all other threads. The primary benefit of
using OpenMP within a node is that it requires very little extra memory per
thread. OpenMP can be turned on or off at configure-time; by default it is
turned on. The only requirement is that the C++ compiler you use must support
the OpenMP 3.1 or higher standard. Most recent compilers do support the use of
OpenMP.

To specify the number of threads at run-time, you can use the ``threads``
argument to :func:`openmc.run`::

  openmc.run(threads=8)

If you're running :ref:`scripts_openmc` directly from the command line, you can
use the ``-s`` or ``--threads`` command-line argument. Alternatively, you can
use the :envvar:`OMP_NUM_THREADS` environment variable. If you do not specify
the number of threads, the OpenMP library will try to determine how many
hardware threads are available on your system and use that many threads.

In general, it is recommended to use as many OpenMP threads as you have hardware
threads on your system. Notably, on a system with Intel hyperthreading, the
hyperthreads should be used and can be expected to provide a 10--30% performance
improvement over not using hyperthreads.

------------------------------------
Distributed-Memory Parallelism (MPI)
------------------------------------

MPI defines a library specification for message-passing between processes. There
are two major implementations of MPI, `OpenMPI <https://www.open-mpi.org/>`_ and
`MPICH <http://www.mpich.org/>`_. Both implementations are known to work with
OpenMC; there is no obvious reason to prefer one over the other. Building OpenMC
with support for MPI requires that you have one of these implementations
installed on your system. For instructions on obtaining MPI, see
:ref:`prerequisites`. Once you have an MPI implementation installed, compile
OpenMC following :ref:`usersguide_compile_mpi`.

To run a simulation using MPI, :ref:`scripts_openmc` needs to be called using
the `mpiexec <https://www.mpich.org/static/docs/v3.1/www1/mpiexec.html>`_
wrapper. For example, to run OpenMC using 32 processes:

.. code-block:: sh

   mpiexec -n 32 openmc

The same thing can be achieved from the Python API by supplying the ``mpi_args``
argument to :func:`openmc.run`::

   openmc.run(mpi_args=['mpiexec', '-n', '32'])

----------------------
Maximizing Performance
----------------------

There are a number of things you can do to ensure that you obtain optimal
performance on a machine when running in parallel:

- **Use OpenMP within each NUMA node**. Some large server processors have so
  many cores that the last level cache is split to reduce memory latency. For
  example, the Intel Xeon Haswell-EP_ architecture uses a snoop mode called
  *cluster on die* where the L3 cache is split in half. Thus, in general, you
  should use one MPI process per socket (and OpenMP within each socket), but for
  these large processors, you will want to go one step further and use one
  process per NUMA node. The Xeon Phi Knights Landing architecture uses a
  similar concept called `sub NUMA clustering
  <https://colfaxresearch.com/knl-numa/>`_.
- **Use a sufficiently large number of particles per generation**. Between
  fission generations, a number of synchronization tasks take place. If the
  number of particles per generation is too low and you are using many
  processes/threads, the synchronization time may become non-negligible.
- **Use hardware threading if available**.
- **Use process binding**. When running with MPI, you should ensure that
  processes are bound_ to a specific hardware region. This can be set using the
  ``-bind-to`` (MPICH) or ``--bind-to`` (OpenMPI) option to ``mpiexec``.
- **Turn off generation of tallies.out**. For large simulations with millions of
  tally bins or more, generating this ASCII file might consume considerable
  time. You can turn off generation of ``tallies.out`` via the
  :attr:`Settings.output` attribute::

     settings = openmc.Settings()
     settings.output = {'tallies': False}

.. _Haswell-EP: http://www.anandtech.com/show/8423/intel-xeon-e5-version-3-up-to-18-haswell-ep-cores-/4
.. _bound: https://wiki.mpich.org/mpich/index.php/Using_the_Hydra_Process_Manager#Process-core_Binding
