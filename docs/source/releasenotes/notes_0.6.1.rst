.. _notes_0.6.1:

==============================
Release Notes for OpenMC 0.6.1
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

- Coarse mesh finite difference acceleration no longer requires PETSc
- Statepoint file numbering is now zero-padded
- Python scripts now compatible with Python 2 or 3
- Ability to run particle restarts in fixed source calculations
- Capability to filter box source by fissionable materials
- Nuclide/element names are now case insensitive in input files
- Improved treatment of resonance scattering for heavy nuclides

---------
Bug Fixes
---------

- 03e890_: Check for energy-dependent multiplicities in ACE files
- 4439de_: Fix distance-to-surface calculation for general plane surface
- 5808ed_: Account for differences in URR band probabilities at different energies
- 2e60c0_: Allow zero atom/weight percents in materials
- 3e0870_: Don't use PWD environment variable when setting path to input files
- dc4776_: Handle probability table resampling correctly
- 01178b_: Fix metastables nuclides in NNDC cross_sections.xml file
- 62ec43_: Don't read tallies.xml when OpenMC is run in plotting mode
- 2a95ef_: Prevent segmentation fault on "current" score without mesh filter

.. _03e890: https://github.com/mit-crpg/openmc/commit/03e890
.. _4439de: https://github.com/mit-crpg/openmc/commit/4439de
.. _5808ed: https://github.com/mit-crpg/openmc/commit/5808ed
.. _2e60c0: https://github.com/mit-crpg/openmc/commit/2e60c0
.. _3e0870: https://github.com/mit-crpg/openmc/commit/3e0870
.. _dc4776: https://github.com/mit-crpg/openmc/commit/dc4776
.. _01178b: https://github.com/mit-crpg/openmc/commit/01178b
.. _62ec43: https://github.com/mit-crpg/openmc/commit/62ec43
.. _2a95ef: https://github.com/mit-crpg/openmc/commit/2a95ef

------------
Contributors
------------

This release contains new contributions from the following people:

- `Sterling Harper <smharper@mit.edu>`_
- `Bryan Herman <bherman@mit.edu>`_
- `Adam Nelson <nelsonag@umich.edu>`_
- `Paul Romano <paul.k.romano@gmail.com>`_
- `Jon Walsh <walshjon@mit.edu>`_
- `Will Boyd <wbinventor@gmail.com>`_
