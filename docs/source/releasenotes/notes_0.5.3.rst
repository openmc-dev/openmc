.. _notes_0.5.3:

==============================
Release Notes for OpenMC 0.5.3
==============================

.. note::
    These release notes are for an upcoming release of OpenMC and are still
    subject to change.

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

- Special run mode --tallies removed.
- Particle restarts and state point restarts are both identified with the -r
  command line flag.
- New regression test suite.
- All memory leaks fixed.
- Shared-memory parallelism with OpenMP.

---------
Bug Fixes
---------

- 2b1e8a_: Normalize direction vector after reflecting particle.
- 5853d2_: Set blank default for cross section listing alias.
- e178c7_: Fix infinite loop with words greater than 80 characters in write_message.
- c18a6e_: Chcek for valid secondary mode on S(a,b) tables.
- 82c456_: Fix bug where last process could have zero particles.

.. _2b1e8a: https://github.com/mit-crpg/openmc/commit/2b1e8a
.. _5853d2: https://github.com/mit-crpg/openmc/commit/5853d2
.. _e178c7: https://github.com/mit-crpg/openmc/commit/e178c7
.. _c18a6e: https://github.com/mit-crpg/openmc/commit/c18a6e
.. _82c456: https://github.com/mit-crpg/openmc/commit/82c456
