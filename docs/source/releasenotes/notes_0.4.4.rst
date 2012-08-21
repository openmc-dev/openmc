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

- Option to specify <upper_right> for tally meshes.

---------
Bug Fixes
---------

- 792eb3_: Fixed degrees of freedom for confidence intervals.

.. _792eb3: https://github.com/mit-crpg/openmc/commit/792eb3
