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
   this. For example, for a boolean_ value, true can be represented either by
   "true" or by "1".

3. Add code to check the variable for any possible errors.

A set of `RELAX NG`_ schemata exists that enables real-time validation of input
files when using the GNU Emacs text editor. You should also modify the RELAX NG
schema for the file you changed (e.g. src/relaxng/geometry.rnc) so that
those who use Emacs can confirm whether their input is valid before they
run. You will need to be familiar with RELAX NG `compact syntax`_.

Working with the FoX Submodule
==============================

The FoX_ library is included as a submodule_ in OpenMC. This means that for a
given commit in OpenMC, there is an associated commit id that links to FoX.
The actual FoX source code is maintained at mit-crpg/fox, branch openmc. When
cloning the OpenMC repo for the first time, you will notice that the directory
*src/xml/fox* is empty. To fetch the submodule source code, you can manually
enter the following from the root directory of OpenMC:

.. code-block:: sh

    git submodule init
    git submodule update

It should be noted that if the submodule is not initialized and updated, *cmake*
will automatically perform these commands if it cannot file the FoX source code.

If you navigate into the FoX source code in OpenMC, src/xml/fox, and check git
information, you will notice that you are in a completely different repo. Actually,
you are in a clone of mit-crpg/fox. If you have write access to this repo, you can 
make changes to the FoX source code, commit and push just like any other repo.
Just because you make changes to the FoX source code in OpenMC or in a standalone
repo, this does not mean that OpenMC will automatically fetch these changes. The
way submodules work is that they are just stored as a commit id. To save FoX xml
source changes to your OpenMC branch, do the following:

1. Go into src/xml/fox and check out the appropriate source code state

2. Navigate back out of fox subdirectory and type:

.. code-block:: sh

    git status

3. Make sure you see that git recognized that the state of FoX changed:

::

    # On branch fox_submodule
    # Changes not staged for commit:
    #   (use "git add <file>..." to update what will be committed)
    #   (use "git checkout -- <file>..." to discard changes in working directory)
    #
    #   modified:   fox (new commits)

4. Commit and push this change

Editing FoX on Personal Fork
============================

If you don't have write access to mit-crpg/fox and thus can't make a branch off of the openmc
branch there, you will need to fork mit-crpg/fox to your personal account. You need to then
link your branch in your OpenMC repo, to the *openmc* branch on your own personal FoX fork.
To do this, edit the *.gitmodules* file in the root folder of the repo. It contains the
following information:

::

    [submodule "src/xml/fox"]
        path = src/xml/fox
        url = git@github.com:mit-crpg/fox

Change the url remote to your own fork. The commit id should stay constant until you start
making modification to FoX yourself. Once you have made changes to your FoX fork and linked
the new commit id to your OpenMC branch, you can pull request your changes in by peforming
the following steps:

1. Create a pull request from your fork of FoX to mit-crpg/fox and wait until it
   is merged into the openmc branch.

2. In your OpenMC repo, change your *.gitmodules* file back to point at mit-crpg/fox.

3. Submit a pull request to mit-crpg/openmc

.. warning:: If you make changes to your FoX submodule inside of an OpenMC repo and do not 
    commit, do **not** run *git submodule update*. This may throw away any changes that
    were not committed.

.. _FoX: https://github.com/mit-crpg/fox
.. _xml_interface module: https://github.com/mit-crpg/openmc/blob/develop/src/xml_interface.F90
.. _input_xml module: https://github.com/mit-crpg/openmc/blob/develop/src/input_xml.F90
.. _XML Schema Part 2: http://www.w3.org/TR/xmlschema-2/
.. _boolean: http://www.w3.org/TR/xmlschema-2/#boolean
.. _RELAX NG: http://relaxng.org/
.. _compact syntax: http://relaxng.org/compact-tutorial-20030326.html
.. _submodule: http://git-scm.com/book/en/Git-Tools-Submodules
