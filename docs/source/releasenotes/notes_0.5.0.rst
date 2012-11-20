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

- Added 'events' score that returns number of events that scored to a tally.
- Restructured tally filter implementation and user input.
- Source convergence acceleration via CMFD (implemented with PETSC).
- Ability to read source files in parallel when number of particles is greater
  than that number of source sites.
- Cone surface types.

---------
Bug Fixes
---------

- 6f8d9d_: Set default tally labels.
- 6a3a5e_: Fix problem with corner-crossing in lattices.

.. _6f8d9d: https://github.com/mit-crpg/openmc/commit/6f8d9d
.. _6a3a5e: https://github.com/mit-crpg/openmc/commit/6a3a5e

