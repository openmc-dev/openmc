========================
cross_sections.xml Files
========================

As a reminder, in order to run a simulation with OpenMC, you will need cross
section data for each nuclide in your problem. OpenMC is not currently
distributed with cross section data, so you will have to obtain cross section
data by other means. The `user's guide`_ offers some helpful advice on how you
can obtain cross sections.

When OpenMC starts up, it needs a cross_sections.xml file that tells it where to
find ACE format cross sections. The files in this directory are configured to
work with a few common cross section sources.

- **cross_sections_ascii.xml** -- This file matches ENDF/B-VII.0 cross sections
    distributed with MCNP5 / MCNP6 beta.

- **cross_sections_nndc.xml** -- This file matches ENDF/B-VII.1 cross sections
    distributed from the `NNDC website`_.

- **cross_sections_serpent.xml** -- This file matches ENDF/B-VII.0 cross
    sections distributed with Serpent 1.1.7.

- **cross_sections.xml** - This file matches ENDF/B-VII.0 cross sections
    distributed with MCNP5 / MCNP6 beta *that have been converted to binary*.

To use any of these files, you need to follow two steps:

1. Change the path on the ``<directory>`` element in the cross_sections.xml file
to the directory containing the ACE files.

2. Enter the absolute path of the cross_sections.xml on the ``<cross_sections>``
element in your settings.xml, or set the CROSS_SECTIONS environment variable to
the full path of the cross_sections.xml file.

.. _user's guide: http://mit-crpg.github.io/openmc/usersguide/install.html#cross-section-configuration
.. _NNDC website: http://www.nndc.bnl.gov/endf/b7.1/acefiles.html
