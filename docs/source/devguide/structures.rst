.. _devguide_structures:

===============
Data Structures
===============

The purpose of this section is to give you an overview of the major data
structures in OpenMC and how they are logically related. A majority of variables
in OpenMC are `derived types`_ (similar to a struct in C). These derived types
are defined in the various header modules, e.g. src/geometry_header.F90. Most
important variables are found in the `global module`_. Have a look through that
module to get a feel for what variables you'll often come across when looking at
OpenMC code.

--------
Particle
--------

Perhaps the variable that you will see most often is simply called ``p`` and is
of type(Particle). This variable stores information about a particle's physical
characteristics (coordinates, direction, energy), what cell and material it's
currently in, how many collisions it has undergone, etc. In practice, only one
particle is followed at a time so there is no array of type(Particle). The
Particle type is defined in the `particle_header module`_.

You will notice that the direction and angle of the particle is stored in a
linked list of type(LocalCoord). In geometries with multiple :ref:`universes`,
the coordinates in each universe are stored in this linked list. If universes or
lattices are not used in a geometry, only one LocalCoord is present in the
linked list.

The LocalCoord type has a component called cell which gives the index in the
``cells`` array in the `global module`_. The ``cells`` array is of type(Cell)
and stored information about each region defined by the user.

----
Cell
----

The Cell type is defined in the `geometry_header module`_ along with other
geometry-related derived types. Each cell in the problem is described in terms
of its bounding surfaces, which are listed on the ``surfaces`` component. The
absolute value of each item in the ``surfaces`` component contains the index of
the corresponding surface in the ``surfaces`` array defined in the `global
module`_. The sign on each item in the ``surfaces`` component indicates whether
the cell exists on the positive or negative side of the surface (see
:ref:`methods_geometry`).

Each cell can either be filled with another universe/lattice or with a
material. If it is filled with a material, the ``material`` component gives the
index of the material in the ``materials`` array defined in the `global
module`_.

-------
Surface
-------

The Surface type is defined in the `geometry_header module`_. A surface is
defined by a type (sphere, cylinder, etc.) and a list of coefficients for that
surface type. The simplest example would be a plane perpendicular to the xy, yz,
or xz plane which needs only one parameter. The ``type`` component indicates the
type through integer parameters such as SURF_SPHERE or SURF_CYL_Y (these are
defined in the `constants module`_). The ``coeffs`` component gives the
necessary coefficients to parameterize the surface type (see
:ref:`surface_element`).

--------
Material
--------

The Material type is defined in the `material_header module`_. Each material
contains a number of nuclides at a given atom density. Each item in the
``nuclide`` component corresponds to the index in the global ``nuclides`` array
(as usual, found in the `global module`_). The ``atom_density`` component is the
same length as the ``nuclides`` component and lists the corresponding atom
density in atom/barn-cm for each nuclide in the ``nuclides`` component.

If the material contains nuclides for which binding effects are important in
low-energy scattering, a :math:`S(\alpha,\beta)` can be associated with that
material through the ``sab_table`` component. Again, this component contains the
index in the ``sab_tables`` array from the `global module`_.

-------
Nuclide
-------

The Nuclide derived type stores cross section and interaction data for a nucleus
and is defined in the `ace_header module`_. The ``energy`` component is an array
that gives the discrete energies at which microscopic cross sections are
tabulated. The actual microscopic cross sections are stored in a separate
derived type, Reaction. An arrays of Reactions is present in the ``reactions``
component. There are a few summary microscopic cross sections stored in other
components, such as ``total``, ``elastic``, ``fission``, and ``nu_fission``.

If a Nuclide is fissionable, the prompt and delayed neutron yield and energy
distributions are also stored on the Nuclide type. Many nuclides also have
unresolved resonance probability table data. If present, this data is stored in
the component ``urr_data`` of derived type UrrData. A complete description of
the probability table method is given in :ref:`probability_tables`.

The list of nuclides present in a problem is stored in the ``nuclides`` array
defined in the `global module`_.

----------
SAlphaBeta
----------

The SAlphaBeta derived type stores :math:`S(\alpha,\beta)` data to account for
molecular binding effects when treating thermal scattering. Each SAlphaBeta
table is associated with a specific nuclide as identified in the ``zaid``
component. A complete description of the :math:`S(\alpha,\beta)` treatment can
be found in :ref:`sab_tables`.

---------
XsListing
---------

The XsListing derived type stores information on the location of an ACE cross
section table based on the data in cross_sections.xml and is defined in the
`ace_header module`_. For each ``<ace_table>`` you see in cross_sections.xml,
there is a XsListing with its information. When the user input is read, the
array ``xs_listings`` in the `global module`_ that is of derived type XsListing
is used to locate the ACE data to parse.

--------------
NuclideMicroXS
--------------

The NuclideMicroXS derived type, defined in the `ace_header module`_, acts as a
'cache' for microscopic cross sections. As a particle is traveling through
different materials, cross sections can be reused if the energy of the particle
hasn't changed. The components ``total``, ``elastic``, ``absorption``,
``fission``, and ``nu_fission`` represent those microscopic cross sections at
the current energy of the particle for a given nuclide. An array ``micro_xs`` in
the `global module`_ that is the same length as the ``nuclides`` array stores
these cached cross sections for each nuclide in the problem.

---------------
MaterialMacroXS
---------------

In addition to the NuclideMicroXS type, there is also a MaterialMacroXS derived
type, defined in the `ace_header module`_ that stored cached *macroscopic* cross
sections for the current material. These macroscopic cross sections are used for
both physics and tallying purposes. The variable ``material_xs`` in the `global
module`_ is of type MaterialMacroXS.


.. _derived types: http://nf.nci.org.au/training/FortranAdvanced/slides/slides.025.html
.. _global module: https://github.com/mit-crpg/openmc/blob/master/src/global.F90
.. _particle_header module: https://github.com/mit-crpg/openmc/blob/master/src/particle_header.F90
.. _geometry_header module: https://github.com/mit-crpg/openmc/blob/master/src/geometry_header.F90
.. _constants module: https://github.com/mit-crpg/openmc/blob/master/src/constants.F90
.. _material_header module: https://github.com/mit-crpg/openmc/blob/master/src/material_header.F90
.. _ace_header module: https://github.com/mit-crpg/openmc/blob/master/src/ace_header.F90
