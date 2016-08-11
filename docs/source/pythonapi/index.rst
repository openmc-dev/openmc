.. _pythonapi:

==========
Python API
==========

OpenMC includes a rich Python API that enables programmatic pre- and
post-processing. The easiest way to begin using the API is to take a look at the
example Jupyter_ notebooks provided. However, this assumes that you are already
familiar with Python and common third-party packages such as NumPy_. If you have
never programmed in Python before, there are many good tutorials available
online. We recommend going through the modules from Codecademy_ and/or the
`Scipy lectures`_. The full API documentation serves to provide more information
on a given module or class.

-------------------------
Example Jupyter Notebooks
-------------------------

.. toctree::
    :maxdepth: 1

    examples/post-processing
    examples/pandas-dataframes
    examples/tally-arithmetic
    examples/mgxs-part-i
    examples/mgxs-part-ii
    examples/mgxs-part-iii
    examples/mgxs-part-iv
    examples/nuclear-data

------------------------------------
:mod:`openmc` -- Basic Functionality
------------------------------------

Handling nuclear data
---------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.XSdata
   openmc.MGXSLibrary


Simulation Settings
-------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Source
   openmc.ResonanceScattering
   openmc.Settings

Material Specification
----------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Nuclide
   openmc.Element
   openmc.Macroscopic
   openmc.Material
   openmc.Materials

Building geometry
-----------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Plane
   openmc.XPlane
   openmc.YPlane
   openmc.ZPlane
   openmc.XCylinder
   openmc.YCylinder
   openmc.ZCylinder
   openmc.Sphere
   openmc.Cone
   openmc.XCone
   openmc.YCone
   openmc.ZCone
   openmc.Quadric
   openmc.Halfspace
   openmc.Intersection
   openmc.Union
   openmc.Complement
   openmc.Cell
   openmc.Universe
   openmc.RectLattice
   openmc.HexLattice
   openmc.Geometry

Many of the above classes are derived from several abstract classes:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Surface
   openmc.Region
   openmc.Lattice

One function is also available to create a hexagonal region defined by the
intersection of six surface half-spaces.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.make_hexagon_region

Constructing Tallies
--------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Filter
   openmc.Mesh
   openmc.Trigger
   openmc.Tally
   openmc.Tallies

Coarse Mesh Finite Difference Acceleration
------------------------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.CMFDMesh
   openmc.CMFD

Plotting
--------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Plot
   openmc.Plots

Running OpenMC
--------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.run
   openmc.plot_geometry

Post-processing
---------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Particle
   openmc.StatePoint
   openmc.Summary

Various classes may be created when performing tally slicing and/or arithmetic:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.arithmetic.CrossScore
   openmc.arithmetic.CrossNuclide
   openmc.arithmetic.CrossFilter
   openmc.arithmetic.AggregateScore
   openmc.arithmetic.AggregateNuclide
   openmc.arithmetic.AggregateFilter

---------------------------------
:mod:`openmc.stats` -- Statistics
---------------------------------

Univariate Probability Distributions
------------------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.stats.Univariate
   openmc.stats.Discrete
   openmc.stats.Uniform
   openmc.stats.Maxwell
   openmc.stats.Watt
   openmc.stats.Tabular
   openmc.stats.Legendre
   openmc.stats.Mixture

Angular Distributions
---------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.stats.UnitSphere
   openmc.stats.PolarAzimuthal
   openmc.stats.Isotropic
   openmc.stats.Monodirectional

Spatial Distributions
---------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.stats.Spatial
   openmc.stats.CartesianIndependent
   openmc.stats.Box
   openmc.stats.Point

----------------------------------------------------------
:mod:`openmc.mgxs` -- Multi-Group Cross Section Generation
----------------------------------------------------------

Energy Groups
-------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.mgxs.EnergyGroups

Multi-group Cross Sections
--------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclassinherit.rst

    openmc.mgxs.MGXS
    openmc.mgxs.AbsorptionXS
    openmc.mgxs.CaptureXS
    openmc.mgxs.Chi
    openmc.mgxs.ChiPrompt
    openmc.mgxs.FissionXS
    openmc.mgxs.InverseVelocity
    openmc.mgxs.KappaFissionXS
    openmc.mgxs.MultiplicityMatrixXS
    openmc.mgxs.NuFissionXS
    openmc.mgxs.NuFissionMatrixXS
    openmc.mgxs.NuScatterXS
    openmc.mgxs.NuScatterMatrixXS
    openmc.mgxs.PromptNuFissionXS
    openmc.mgxs.ScatterXS
    openmc.mgxs.ScatterMatrixXS
    openmc.mgxs.TotalXS
    openmc.mgxs.TransportXS

Multi-group Cross Section Libraries
-----------------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.mgxs.Library

-------------------------------------
:mod:`openmc.model` -- Model Building
-------------------------------------

TRISO Fuel Modeling
-------------------

Classes
+++++++

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.model.TRISO

Functions
+++++++++

.. autosummary::
   :toctree: generated
   :nosignatures:

   openmc.model.create_triso_lattice

--------------------------------------------
:mod:`openmc.data` -- Nuclear Data Interface
--------------------------------------------

Physical Data
-------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myfunction.rst

    openmc.data.atomic_mass

Core Classes
------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.data.IncidentNeutron
    openmc.data.Reaction
    openmc.data.Product
    openmc.data.Tabulated1D
    openmc.data.ThermalScattering
    openmc.data.CoherentElastic

Angle-Energy Distributions
--------------------------

.. autosummary::
    :toctree: generated
    :nosignatures:
    :template: myclass.rst

    openmc.data.AngleEnergy
    openmc.data.KalbachMann
    openmc.data.CorrelatedAngleEnergy
    openmc.data.UncorrelatedAngleEnergy
    openmc.data.NBodyPhaseSpace
    openmc.data.AngleDistribution
    openmc.data.EnergyDistribution
    openmc.data.ArbitraryTabulated
    openmc.data.GeneralEvaporation
    openmc.data.MaxwellEnergy
    openmc.data.Evaporation
    openmc.data.WattEnergy
    openmc.data.MadlandNix
    openmc.data.DiscretePhoton
    openmc.data.LevelInelastic
    openmc.data.ContinuousTabular

ACE Format
----------

Classes
+++++++

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.data.ace.Library
   openmc.data.ace.Table

Functions
+++++++++

.. autosummary::
   :toctree: generated
   :nosignatures:

   openmc.data.ace.ascii_to_binary

.. _Jupyter: https://jupyter.org/
.. _NumPy: http://www.numpy.org/
.. _Codecademy: https://www.codecademy.com/tracks/python
.. _Scipy lectures: https://scipy-lectures.github.io/
