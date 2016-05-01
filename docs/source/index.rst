===========================
The OpenMC Monte Carlo Code
===========================

OpenMC is a Monte Carlo particle transport simulation code focused on neutron
criticality calculations. It is capable of simulating 3D models based on
constructive solid geometry with second-order surfaces. OpenMC supports either
continuous-energy or multi-group transport.  The continuous-energy
particle interaction data is based on ACE format cross sections, also used
in the MCNP and Serpent Monte Carlo codes.

OpenMC was originally developed by members of the `Computational Reactor Physics
Group`_ at the `Massachusetts Institute of Technology`_ starting
in 2011. Various universities, laboratories, and other organizations now
contribute to the development of OpenMC. For more information on OpenMC, feel
free to send a message to the User's Group `mailing list`_. Documentation for
the latest developmental version of the develop branch can be found on
`Read the Docs`_.

.. _Computational Reactor Physics Group: http://crpg.mit.edu
.. _Massachusetts Institute of Technology: http://web.mit.edu
.. _mailing list: https://groups.google.com/forum/?fromgroups=#!forum/openmc-users
.. _Read the Docs: http://openmc.readthedocs.io/en/latest/

.. only:: html

    --------
    Contents
    --------

.. toctree::
    :maxdepth: 1

    quickinstall
    releasenotes
    methods/index
    usersguide/index
    devguide/index
    pythonapi/index
    publications
    license
    developers
