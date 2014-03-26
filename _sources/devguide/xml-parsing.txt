.. _devguide_xml-parsing:

=================
XML Input Parsing
=================

OpenMC relies on the FoX_ Fortran XML library for reading and intrepreting the
XML input files for geometry, materials, settings, tallies, etc. The use of an
XML format makes writing input files considerably more flexible than would
otherwise be possible.

With the FoX library, extending the user input files to include new tags is
fairly straightforward. The steps for modifying/adding input are as follows:

1. Add appropriate calls to procedures from the `xml_interface module`_, such as
``check_for_node``, ``get_node_value``, and ``get_node_array``. All input
reading is performed in the `input_xml module`_.

2. Make sure that your input can be categorized as one of the datatypes from
`XML Schema Part 2`_ and that parsing of the data appropriately reflects
this. For example, for a boolean_ value, true can be represented either by "true"
or by "1".

3. Add code to check the variable for any possible errors.

A set of `RELAX NG`_ schemata exists that enables real-time validation of input
files when using the GNU Emacs text editor. You should also modify the RELAX NG
schema for the file you changed (e.g. src/relaxng/geometry.rnc) so that
those who use Emacs can confirm whether their input is valid before they
run. You will need to be familiar with RELAX NG `compact syntax`_.

.. _FoX: https://github.com/andreww/fox
.. _xml_interface module: https://github.com/mit-crpg/openmc/blob/develop/src/xml_interface.F90
.. _input_xml module: https://github.com/mit-crpg/openmc/blob/develop/src/input_xml.F90
.. _XML Schema Part 2: http://www.w3.org/TR/xmlschema-2/
.. _boolean: http://www.w3.org/TR/xmlschema-2/#boolean
.. _RELAX NG: http://relaxng.org/
.. _compact syntax: http://relaxng.org/compact-tutorial-20030326.html
