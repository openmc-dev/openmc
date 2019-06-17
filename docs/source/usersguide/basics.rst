.. _usersguide_basics:

======================
Basics of Using OpenMC
======================

----------------
Creating a Model
----------------

When you build and install OpenMC, you will have an :ref:`scripts_openmc`
executable on your system. When you run ``openmc``, the first thing it will do
is look for a set of XML_ files that describe the model you want to
simulation. Three of these files are required and another three are optional, as
described below.

.. admonition:: Required
   :class: error

   :ref:`io_materials`
     This file describes what materials are present in the problem and what they
     are composed of. Additionally, it indicates where OpenMC should look for a
     cross section library.

   :ref:`io_geometry`
     This file describes how the materials defined in ``materials.xml`` occupy
     regions of space. Physical volumes are defined using constructive solid
     geometry, described in detail in :ref:`usersguide_geometry`.

   :ref:`io_settings`
     This file indicates what mode OpenMC should be run in, how many particles
     to simulate, the source definition, and a whole host of miscellaneous
     options.

.. admonition:: Optional
   :class: note

   :ref:`io_tallies`
     This file describes what physical quantities should be tallied during the
     simulation (fluxes, reaction rates, currents, etc.).

   :ref:`io_plots`
     This file gives specifications for producing slice or voxel plots of the
     geometry.

eXtensible Markup Language (XML)
--------------------------------

Unlike many other Monte Carlo codes which use an arbitrary-format ASCII file
with "cards" to specify a particular geometry, materials, and associated run
settings, the input files for OpenMC are structured in a set of `XML
<http://www.w3.org/XML/>`_ files. XML, which stands for eXtensible Markup
Language, is a simple format that allows data to be exchanged efficiently
between different programs and interfaces.

Anyone who has ever seen webpages written in HTML will be familiar with the
structure of XML whereby "tags" enclosed in angle brackets denote that a
particular piece of data will follow. Let us examine the follow example:

.. code-block:: xml

    <person>
      <firstname>John</firstname>
      <lastname>Smith</lastname>
      <age>27</age>
      <occupation>Health Physicist</occupation>
    </person>

Here we see that the first tag indicates that the following data will describe a
person. The nested tags *firstname*, *lastname*, *age*, and *occupation*
indicate characteristics about the person being described.

In much the same way, OpenMC input uses XML tags to describe the geometry, the
materials, and settings for a Monte Carlo simulation. Note that because the XML
files have a well-defined structure, they can be validated using the
:ref:`scripts_validate` script or using :ref:`Emacs nXML mode
<usersguide_nxml>`.

Creating Input Files
--------------------

.. currentmodule:: openmc

The most rudimentary option for creating input files is to simply write them
from scratch using the :ref:`XML format specifications
<io_file_formats_input>`. This approach will feel familiar to users of other
Monte Carlo codes such as MCNP and Serpent, with the added bonus that the XML
formats feel much more "readable". Alternatively, input files can be generated
using OpenMC's :ref:`Python API <pythonapi>`, which is introduced in the
following section.

----------
Python API
----------

OpenMC's :ref:`Python API <pythonapi>` defines a set of functions and classes
that roughly correspond to elements in the XML files. For example, the
:class:`openmc.Cell` Python class directly corresponds to the
:ref:`cell_element` in XML. Each XML file itself also has a corresponding class:
:class:`openmc.Geometry` for ``geometry.xml``, :class:`openmc.Materials` for
``materials.xml``, :class:`openmc.Settings` for ``settings.xml``, and so on. To
create a model then, one creates instances of these classes and then uses the
``export_to_xml()`` method, e.g., :meth:`Geometry.export_to_xml`. Most scripts
that generate a full model will look something like the following:

.. code-block:: Python

   # Create materials
   materials = openmc.Materials()
   ...
   materials.export_to_xml()

   # Create geometry
   geometry = openmc.Geometry()
   ...
   geometry.export_to_xml()

   # Assign simulation settings
   settings = openmc.Settings()
   ...
   settings.export_to_xml()

Once a model has been created and exported to XML, a simulation can be run either
by calling :ref:`scripts_openmc` directly from a shell or by using the
:func:`openmc.run()` function from Python.

Identifying Objects
-------------------

In the XML user input files, each object (cell, surface, tally, etc.) has to be
uniquely identified by a positive integer (ID) in the same manner as MCNP and
Serpent. In the Python API, integer IDs can be assigned but it is not strictly
required. When IDs are not explicitly assigned to instances of the OpenMC Python
classes, they will be automatically assigned.

.. _result_files:

-----------------------------
Viewing and Analyzing Results
-----------------------------

After a simulation has been completed by running :ref:`scripts_openmc`, you will
have several output files that were created:

``tallies.out``
  An ASCII file showing the mean and standard deviation of the mean for any
  user-defined tallies.

``summary.h5``
  An HDF5 file with a complete description of the geometry and materials used in
  the simulation.

``statepoint.#.h5``
  An HDF5 file with the complete results of the simulation, including tallies as
  well as the final source distribution. This file can be used both to
  view/analyze results as well as restart a simulation if desired.

For a simple simulation with few tallies, looking at the ``tallies.out`` file
might be sufficient. For anything more complicated (plotting results, finding a
subset of results, etc.), you will likely find it easier to work with the
statepoint file directly using the :class:`openmc.StatePoint` class. For more
details on working with statepoints, see :ref:`usersguide_statepoint`.

--------------
Physical Units
--------------

Unless specified otherwise, all length quantities are assumed to be in units of
centimeters, all energy quantities are assumed to be in electronvolts, and all
time quantities are assumed to be in seconds.

======= ============ ======
Measure Default unit Symbol
======= ============ ======
length  centimeter   cm
energy  electronvolt eV
time    second       s
======= ============ ======

------------------------------------
ERSN-OpenMC Graphical User Interface
------------------------------------

A third-party Java-based user-friendly graphical user interface for creating XML
input files called ERSN-OpenMC_ is developed and maintained by members of the
Radiation and Nuclear Systems Group at the Faculty of Sciences Tetouan, Morocco.
The GUI also allows one to automatically download prerequisites for installing and
running OpenMC.

.. _ERSN-OpenMC: https://github.com/EL-Bakkali-Jaafar/ERSN-OpenMC
