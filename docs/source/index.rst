===========================
The OpenMC Monte Carlo Code
===========================

OpenMC is a community-developed Monte Carlo neutron and photon transport
simulation code. It is capable of performing fixed source, k-eigenvalue, and
subcritical multiplication calculations on models built using either a
constructive solid geometry or CAD representation. OpenMC supports both
continuous-energy and multigroup transport. The continuous-energy particle
interaction data is based on a native HDF5 format that can be generated from ACE
files produced by NJOY. Parallelism is enabled via a hybrid MPI and OpenMP
programming model.

OpenMC was originally developed by members of the `Computational Reactor Physics
Group <http://crpg.mit.edu>`_ at the `Massachusetts Institute of Technology
<http://web.mit.edu>`_ starting in 2011. Various universities, laboratories, and
other organizations now contribute to the development of OpenMC. For more
information on OpenMC, feel free to send a message to the User's Group `mailing
list <https://groups.google.com/forum/?fromgroups=#!forum/openmc-users>`_.

.. admonition:: Recommended publication for citing
   :class: tip

   Paul K. Romano, Nicholas E. Horelik, Bryan R. Herman, Adam G. Nelson, Benoit
   Forget, and Kord Smith, "`OpenMC: A State-of-the-Art Monte Carlo Code for
   Research and Development <https://doi.org/10.1016/j.anucene.2014.07.048>`_,"
   *Ann. Nucl. Energy*, **82**, 90--97 (2015).

.. only:: html

   --------
   Contents
   --------

.. toctree::
    :maxdepth: 1

    quickinstall
    examples/index
    releasenotes/index
    methods/index
    usersguide/index
    devguide/index
    pythonapi/index
    capi/index
    io_formats/index
    publications
    license
