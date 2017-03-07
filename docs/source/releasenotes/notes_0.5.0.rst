.. _notes_0.5.0:

==============================
Release Notes for OpenMC 0.5.0
==============================

-------------------
System Requirements
-------------------

There are no special requirements for running the OpenMC code. As of this
release, OpenMC has been tested on a variety of Linux distributions, Mac OS X,
and Microsoft Windows 7. Memory requirements will vary depending on the size of
the problem at hand (mostly on the number of nuclides in the problem).

------------
New Features
------------

- All user input options that formerly accepted "off" or "on" should now be
  "false" or "true" (the proper XML schema datatype).
- The <criticality> element is deprecated and was replaced with <eigenvalue>.
- Added 'events' score that returns number of events that scored to a tally.
- Restructured tally filter implementation and user input.
- Source convergence acceleration via CMFD (implemented with PETSc).
- Ability to read source files in parallel when number of particles is greater
  than that number of source sites.
- Cone surface types.

---------
Bug Fixes
---------

- 737b90_: Coincident surfaces from separate universes / particle traveling
  tangent to a surface.
- a819b4_: Output of surface neighbors in summary.out file.
- b11696_: Reading long attribute lists in XML input.
- 2bd46a_: Search for tallying nuclides when no default_xs specified.
- 7a1f08_: Fix word wrapping when writing messages.
- c0e3ec_: Prevent underflow when compiling with MPI=yes and DEBUG=yes.
- 6f8d9d_: Set default tally labels.
- 6a3a5e_: Fix problem with corner-crossing in lattices.

.. _737b90: https://github.com/mit-crpg/openmc/commit/737b90
.. _a819b4: https://github.com/mit-crpg/openmc/commit/a819b4
.. _b11696: https://github.com/mit-crpg/openmc/commit/b11696
.. _2bd46a: https://github.com/mit-crpg/openmc/commit/2bd46a
.. _7a1f08: https://github.com/mit-crpg/openmc/commit/7a1f08
.. _c0e3ec: https://github.com/mit-crpg/openmc/commit/c0e3ec
.. _6f8d9d: https://github.com/mit-crpg/openmc/commit/6f8d9d
.. _6a3a5e: https://github.com/mit-crpg/openmc/commit/6a3a5e

