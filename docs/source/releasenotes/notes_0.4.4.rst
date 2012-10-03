.. _notes_0.4.4:

==============================
Release Notes for OpenMC 0.4.4
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

- Suppress creation of summary.out and cross_sections.out by default with option
  to turn them on with <output> tag in settings.xml file.
- Ability to create HDF5 state points.
- Binary source file is now part of state point file by default.
- Enhanced state point usage and added state point Python scripts.
- Turning confidence intervals on affects k-effective.
- Option to specify <upper_right> for tally meshes.

---------
Bug Fixes
---------

- 7ee416_: Fixed bug with multi-line input using type='word'.
- 792eb3_: Fixed degrees of freedom for confidence intervals.
- 7fd617_: Fixed bug with restart runs in parallel.
- dc4a8f_: Fixed bug with fixed source restart runs.

.. _7ee416: https://github.com/mit-crpg/openmc/commit/7ee416
.. _792eb3: https://github.com/mit-crpg/openmc/commit/792eb3
.. _7fd617: https://github.com/mit-crpg/openmc/commit/7fd617
.. _dc4a8f: https://github.com/mit-crpg/openmc/commit/dc4a8f
