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
   openmc.ZernikeRadialFilter
   openmc.ParticleFilter
   openmc.Mesh
   openmc.Trigger
   openmc.TallyDerivative
   openmc.Tally
   openmc.Tallies

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

Coarse Mesh Finite Difference Acceleration
------------------------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   openmc.CMFDMesh
   openmc.CMFDRun

CMFD is implemented in OpenMC and allows users to accelerate fission source
convergence during inactive neutron batches. To run CMFD, the CMFDRun class should
be used. The following properties can be set through the CMFDRun class:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.CMFDRun.cmfd_begin
   openmc.CMFDRun.dhat_reset
   openmc.CMFDRun.cmfd_display
   openmc.CMFDRun.cmfd_downscatter
   openmc.CMFDRun.cmfd_feedback
   openmc.CMFDRun.cmfd_ktol
   openmc.CMFDRun.cmfd_mesh
   openmc.CMFDRun.norm
   openmc.CMFDRun.cmfd_adjoint_type
   openmc.CMFDRun.cmfd_power_monitor
   openmc.CMFDRun.cmfd_run_adjoint
   openmc.CMFDRun.cmfd_shift
   openmc.CMFDRun.cmfd_stol
   openmc.CMFDRun.cmfd_spectral
   openmc.CMFDRun.cmfd_reset
   openmc.CMFDRun.cmfd_write_matrices
   openmc.CMFDRun.gauss_seidel_tolerance

At the minimum, a CMFD mesh needs to be specified in order to run CMFD. Once
these properties are set, an OpenMC simulation can be run with CMFD turned on
with the function:

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   openmc.CMFDRun.run

