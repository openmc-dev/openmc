.. _notes_0.5.4:

==============================
Release Notes for OpenMC 0.5.4
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

- Source sites outside geometry are resampled
- XML-Fortran backend replaced by FoX XML
- Ability to write particle track files
- Handle lost particles more gracefully (via particle track files)
- Multiple random number generator streams
- Mesh tally plotting utility converted to use Tkinter rather than PyQt
- Script added to download ACE data from NNDC
- Mixed ASCII/binary cross_sections.xml now allowed
- Expanded options for writing source bank
- Re-enabled ability to use source file as starting source
- S(a,b) recalculation avoided when same nuclide and S(a,b) table are accessed

---------
Bug Fixes
---------

- 32c03c_: Check for valid data in cross_sections.xml
- c71ef5_: Fix bug in statepoint.py
- 8884fb_: Check for all ZAIDs for S(a,b) tables
- b38af0_: Fix XML reading on multiple levels of input
- d28750_: Fix bug in convert_xsdir.py
- cf567c_: ENDF/B-VI data checked for compatibility
- 6b9461_: Fix p_valid sampling inside of sample_energy

.. _32c03c: https://github.com/mit-crpg/openmc/commit/32c03c
.. _c71ef5: https://github.com/mit-crpg/openmc/commit/c71ef5
.. _8884fb: https://github.com/mit-crpg/openmc/commit/8884fb
.. _b38af0: https://github.com/mit-crpg/openmc/commit/b38af0
.. _d28750: https://github.com/mit-crpg/openmc/commit/d28750
.. _cf567c: https://github.com/mit-crpg/openmc/commit/cf567c
.. _6b9461: https://github.com/mit-crpg/openmc/commit/6b9461

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
