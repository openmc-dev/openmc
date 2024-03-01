------------------------------------------------------
:mod:`openmc.lib` -- Python bindings to the C/C++ API
------------------------------------------------------

.. automodule:: openmc.lib

Functions
---------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   broaden_wmp_polynomials
   calc_pn
   calc_rn
   calc_zn_rad
   calculate_volumes
   cell
   cells
   current_batch
   evaluate_legendre
   export_properties
   export_weight_windows
   filter
   filters
   finalize
   find_cell
   find_material
   global_bounding_box
   global_tallies
   hard_reset
   id_map
   import_properties
   import_weight_windows
   init
   is_initialized
   is_statepoint_batch
   iter_batches
   keff
   load_nuclide
   master
   material
   materials
   maxwell_spectrum
   mesh
   meshes
   next_batch
   normal_variate
   nuclide
   nuclides
   num_realizations
   plot
   plot_geometry
   property_map
   reset
   reset_timers
   rotate_angle
   run
   run_in_memory
   sample_external_source
   settings
   simulation_finalize
   simulation_init
   source_bank
   statepoint_write
   tallies
   tally
   watt_spectrum
   weight_windows

Classes
-------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   AzimuthalFilter
   CDLL
   Cell
   CellFilter
   CellInstanceFilter
   CellbornFilter
   CellfromFilter
   CollisionFilter
   CylindricalMesh
   DelayedGroupFilter
   DistribcellFilter
   EnergyFilter
   EnergyFunctionFilter
   EnergyoutFilter
   Filter
   LegendreFilter
   Material
   MaterialFilter
   MaterialFromFilter
   MeshFilter
   MeshSurfaceFilter
   MuFilter
   Nuclide
   POINTER
   ParticleFilter
   PathLike
   PolarFilter
   RectilinearMesh
   RegularMesh
   SpatialLegendreFilter
   SphericalHarmonicsFilter
   SphericalMesh
   Structure
   SurfaceFilter
   Tally
   UniverseFilter
   UnstructuredMesh
   WeightWindows
   ZernikeFilter
   ZernikeRadialFilter