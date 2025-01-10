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

   openmc.SourceBase
   openmc.IndependentSource
   openmc.FileSource
   openmc.CompiledSource
   openmc.MeshSource
   openmc.SourceParticle
   openmc.VolumeCalculation
   openmc.Settings

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.read_source_file
   openmc.write_source_file
   openmc.wwinp_to_wws

Material Specification
----------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

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
   openmc.XTorus
   openmc.YTorus
   openmc.ZTorus
   openmc.Halfspace
   openmc.Intersection
   openmc.Union
   openmc.Complement
   openmc.BoundingBox
   openmc.Cell
   openmc.Universe
   openmc.DAGMCUniverse
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
   openmc.MaterialFromFilter
   openmc.CellFilter
   openmc.CellFromFilter
   openmc.CellBornFilter
   openmc.CellInstanceFilter
   openmc.CollisionFilter
   openmc.SurfaceFilter
   openmc.MeshFilter
   openmc.MeshBornFilter
   openmc.ImportanceFilter
   openmc.MeshSurfaceFilter
   openmc.EnergyFilter
   openmc.EnergyoutFilter
   openmc.MuFilter
   openmc.MuSurfaceFilter
   openmc.PolarFilter
   openmc.AzimuthalFilter
   openmc.DistribcellFilter
   openmc.DelayedGroupFilter
   openmc.EnergyFunctionFilter
   openmc.LegendreFilter
   openmc.SpatialLegendreFilter
   openmc.SphericalHarmonicsFilter
   openmc.TimeFilter
   openmc.ZernikeFilter
   openmc.ZernikeRadialFilter
   openmc.ParticleFilter
   openmc.RegularMesh
   openmc.RectilinearMesh
   openmc.CylindricalMesh
   openmc.SphericalMesh
   openmc.UnstructuredMesh
   openmc.Trigger
   openmc.TallyDerivative
   openmc.SensitivityTally
   openmc.Sensitivity
   openmc.Tally
   openmc.Tallies

Geometry Plotting
-----------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.Plot
   openmc.ProjectionPlot
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
   openmc.ParticleList
   openmc.ParticleTrack
   openmc.StatePoint
   openmc.Summary
   openmc.Track
   openmc.Tracks

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.voxel_to_vtk

The following classes and functions are used for functional expansion reconstruction.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.ZernikeRadial

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.legendre_from_expcoef


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

Variance Reduction
------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass

   openmc.WeightWindows
   openmc.WeightWindowGenerator
   openmc.hdf5_to_wws


Coarse Mesh Finite Difference Acceleration
------------------------------------------

CMFD is implemented in OpenMC and allows users to accelerate fission source
convergence during inactive neutron batches. To use CMFD, the
:class:`openmc.cmfd.CMFDRun` class executes OpenMC through the C API, solving
the CMFD system between fission generations and modifying the source weights.
Note that the :mod:`openmc.cmfd` module is not imported by default with the
:mod:`openmc` namespace and needs to be imported explicitly.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.cmfd.CMFDMesh
   openmc.cmfd.CMFDRun

At the minimum, a CMFD mesh needs to be specified in order to run CMFD. Once the
mesh and other optional properties are set, a simulation can be run with CMFD
turned on using :meth:`openmc.cmfd.CMFDRun.run`.
