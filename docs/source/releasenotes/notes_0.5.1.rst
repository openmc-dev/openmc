.. _notes_0.5.1:

==============================
Release Notes for OpenMC 0.5.1
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

- Absorption and combined estimators for k-effective.
- Natural elements can now be specified in materials using <element> rather than
  <nuclide>.
- Support for multiple S(a,b) tables in a single material (e.g. BeO).
- Test suite using Python nosetests.
- Proper install capability with 'make install'.
- Lattices can now be 2 or 3 dimensions.
- New scatter-PN score type.
- New kappa-fission score type.
- Ability to tally any reaction by specifying MT.

---------
Bug Fixes
---------

- 94103e_: Two checks for outgoing energy filters.
- e77059_: Fix reaction name for MT=849.
- b0fe88_: Fix distance to surface for cones.
- 63bfd2_: Fix tracklength tallies with cell filter and universes.
- 88daf7_: Fix analog tallies with survival biasing.

.. _94103e: https://github.com/mit-crpg/openmc/commit/94103e
.. _e77059: https://github.com/mit-crpg/openmc/commit/e77059
.. _b0fe88: https://github.com/mit-crpg/openmc/commit/b0fe88
.. _63bfd2: https://github.com/mit-crpg/openmc/commit/63bfd2
.. _88daf7: https://github.com/mit-crpg/openmc/commit/88daf7
