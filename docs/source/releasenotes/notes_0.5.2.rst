.. _notes_0.5.2:

==============================
Release Notes for OpenMC 0.5.2
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

- Python script for mesh tally plotting
- Isotopic abundances based on IUPAC 2009 when using <element>
- Particle restart file
- Code will abort after certain number of lost particles (defaults to 10)
- Region outside lattice can be filled with material (void by default)
- 3D voxel plots
- Full HDF5/PHDF5 support (including support in statepoint.py)
- Cell overlap checking with -g command line flag (or when plotting)

---------
Bug Fixes
---------

- 7632f3_: Fixed bug in statepoint.py for multiple generations per batch.
- f85ac4_: Fix infinite loop bug in error module.
- 49c36b_: Don't convert surface ids if surface filter is for current tallies.
- 5ccc78_: Fix bug in reassignment of bins for mesh filter.
- b1f52f_: Fixed bug in plot color specification.
- eae7e5_: Fixed many memory leaks.
- 10c1cc_: Minor CMFD fixes.
- afdb50_: Add compatibility for XML comments without whitespace.
- a3c593_: Fixed bug in use of free gas scattering for H-1.
- 3a66e3_: Fixed bug in 2D mesh tally implementation.
- ab0793_: Corrected PETSC_NULL references to their correct types.

.. _7632f3: https://github.com/mit-crpg/openmc/commit/7632f3
.. _f85ac4: https://github.com/mit-crpg/openmc/commit/f85ac4
.. _49c36b: https://github.com/mit-crpg/openmc/commit/49c36b
.. _5ccc78: https://github.com/mit-crpg/openmc/commit/5ccc78
.. _b1f52f: https://github.com/mit-crpg/openmc/commit/b1f52f
.. _eae7e5: https://github.com/mit-crpg/openmc/commit/eae7e5
.. _10c1cc: https://github.com/mit-crpg/openmc/commit/10c1cc
.. _afdb50: https://github.com/mit-crpg/openmc/commit/afdb50
.. _a3c593: https://github.com/mit-crpg/openmc/commit/a3c593
.. _3a66e3: https://github.com/mit-crpg/openmc/commit/3a66e3
.. _ab0793: https://github.com/mit-crpg/openmc/commit/ab0793
