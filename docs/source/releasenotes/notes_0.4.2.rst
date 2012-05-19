.. _notes_0.4.2:

==============================
Release Notes for OpenMC 0.4.2
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

- Option to not reduce tallies across processors at end of each batch.
- Uniform fission site method for reducing variance on local tallies.
- Reading/writing binary source files.
- Added more messages for <trace> or high verbosity.
- Estimator for diffusion coefficient.
- Ability to specify 'point' source type.
- Ability to change random number seed.
- User's can now specify units='sum' on a <density> tag. This tells the code
  that the total material density is the sum of the atom fractions listed for
  each nuclide on the material.

---------
Bug Fixes
---------

- 671f30_: Fixed surface currents on mesh not encompassing geometry.
- b2c40e_: Fixed bug in incoming energy filter for track-length tallies.
- 5524fd_: Mesh filter now works with track-length tallies.
- d050c7_: Added Bessel's correction to make estimate of variance unbiased.
- 2a5b9c_: Fixed regression in plotting.

.. _671f30: https://github.com/mit-crpg/openmc/commit/671f30
.. _b2c40e: https://github.com/mit-crpg/openmc/commit/b2c40e
.. _5524fd: https://github.com/mit-crpg/openmc/commit/5524fd
.. _d050c7: https://github.com/mit-crpg/openmc/commit/d050c7
.. _2a5b9c: https://github.com/mit-crpg/openmc/commit/2a5b9c
