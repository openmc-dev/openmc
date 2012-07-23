.. _notes_0.4.3:

==============================
Release Notes for OpenMC 0.4.3
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

- Ability to tally reaction rates for individual nuclides within a material.
- Reduced memory usage by removing redundant storage or some cross-sections.
- 3bd35b_: Log-log interpolation for URR probability tables.
- Support to specify labels on tallies (nelsonag_).

---------
Bug Fixes
---------

- 3212f5_: Fixed issue with blank line at beginning of XML files.

.. _nelsonag: https://github.com/nelsonag
.. _3bd35b: https://github.com/mit-crpg/openmc/commit/3bd35b
.. _3212f5: https://github.com/mit-crpg/openmc/commit/3212f5
