.. _notes_0.5.4:

==============================
Release Notes for OpenMC 0.5.4
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

- New XML parsing backend (FoX)
- Ability to write particle track files
- Handle lost particles more gracefully (via particle track files)
- Source sites outside geometry are resampled
- Multiple random number generator streams
- plot_mesh_tally.py utility converted to use Tkinter rather than PyQt
- Script added to download ACE data from NNDC
- Mixed ASCII/binary cross_sections.xml now allowed
- Expanded options for writing source bank
- Re-enabled ability to use source file as starting source

---------
Bug Fixes
---------

- 32c03c_: Check for valid data in cross_sections.xml
- c71ef5_: Fix bug in statepoint.py
- 8884fb_: Check for all ZAIDs for S(a,b) tables
- b38af0_: Fix XML reading on multiple levels of input
- d28750_: Fix bug in convert_xsdir.py

.. _32c03c: https://github.com/mit-crpg/openmc/commit/32c03c
.. _c71ef5: https://github.com/mit-crpg/openmc/commit/c71ef5
.. _8884fb: https://github.com/mit-crpg/openmc/commit/8884fb
.. _b38af0: https://github.com/mit-crpg/openmc/commit/b38af0
.. _d28750: https://github.com/mit-crpg/openmc/commit/d28750

------------
Contributors
------------

This release contains new contributions from the following people:

- `Sterling Harper <smharper@mit.edu>`_
- `Bryan Herman <bherman@mit.edu>`_
- `Nick Horelik <nhorelik@mit.edu>`_
- `Adam Nelson <nelsonag@umich.edu>`_
- `Paul Romano <paul.k.romano@gmail.com>`_
- `Tuomas Viitanen <tuomas.viitanen@vtt.fi>`_
- `Jon Walsh <walshjon@mit.edu>`_
