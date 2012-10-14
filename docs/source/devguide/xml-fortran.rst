.. _devguide_xml-fortran:

=========================
xml-fortran Input Parsing
=========================

OpenMC relies on the xml-fortran package for reading and intrepreting the XML
input files for geometry, materials, settings, tallies, etc. The use of an XML
format makes writing input files considerably more flexible than would otherwise
be possible.

With the xml-fortran package, extending the user input files to include new tags
is fairly straightforward. A "template" file exists for each diferent type of
input file that tells xml-fortran what to expect in a file. These template files
can be found in the src/templates directory. The steps for modifying/adding
input are as follows:

1. Add a ``<variable>``` tag to the desired template file,
e.g. src/templates/geometry_t.xml. See the `xml-fortran documentation`_ for a
description of the acceptable fields.

2. In the input_xml module, any input given in your new tag will be read
automatically through a call to, e.g. read_xml_file_geometry_t. Whatever
variable name you specified should have the data available.

3. Add code in the appropriate subroutine to check the variable for any possible
errors.

4. Add a variable in OpenMC to copy the temporary variable into if there are no
errors.

A set of `RELAX NG`_ schemata exists that enables real-time validation of input
files when using the GNU Emacs text editor. You should also modify the RELAX NG
schema for the template you changed (e.g. src/templates/geometry.rnc) so that
those who use Emacs can confirm whether their input is valid before they
run. You will need to be familiar with RELAX NG `compact syntax`_.

.. _xml-fortran documentation: http://xml-fortran.sourceforge.net/documentation.html
.. _RELAX NG: http://relaxng.org/
.. _compact syntax: http://relaxng.org/compact-tutorial-20030326.html
