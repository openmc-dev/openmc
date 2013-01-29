.. _notes_0.5.1:

==============================
Release Notes for OpenMC 0.5.1
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

- Proper install capability with 'make install'.
- Lattices can now be 2 or 3 dimensions.
- New scatter-PN score type.
- New kappa-fission score type.
- Ability to tally any reaction by specifying MT.

---------
Bug Fixes
---------

- 63bfd2_: Fix tracklength tallies with cell filter and universes.
- 88daf7_: Fix analog tallies with survival biasing.

.. _63bfd2: https://github.com/mit-crpg/openmc/commit/63bfd2
.. _88daf7: https://github.com/mit-crpg/openmc/commit/88daf7
