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
   openmc.VolumeCalculation
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

Cross sections for nuclides, elements, and materials can be plotted using the
following function:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.plot_xs

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

.. _pythonapi_tallies:

Constructing Tallies
--------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Filter
   openmc.UniverseFilter
   openmc.MaterialFilter
   openmc.CellFilter
   openmc.CellFromFilter
   openmc.CellbornFilter
   openmc.SurfaceFilter
   openmc.MeshFilter
   openmc.MeshSurfaceFilter
   openmc.EnergyFilter
   openmc.EnergyoutFilter
   openmc.MuFilter
   openmc.PolarFilter
   openmc.AzimuthalFilter
   openmc.DistribcellFilter
   openmc.DelayedGroupFilter
   openmc.EnergyFunctionFilter
   openmc.LegendreFilter
   openmc.SpatialLegendreFilter
   openmc.SphericalHarmonicsFilter
   openmc.ZernikeFilter
   openmc.Mesh
   openmc.Trigger
   openmc.TallyDerivative
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

Geometry Plotting
-----------------

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
   openmc.calculate_volumes
   openmc.plot_geometry
   openmc.plot_inline
   openmc.search_for_keff

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
