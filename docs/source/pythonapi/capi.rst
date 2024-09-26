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

   calculate_volumes
   current_batch
   current_surface_file
   export_properties
   export_weight_windows
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
   is_statepoint_batch
   iter_batches
   keff
   load_nuclide
   master
   next_batch
   num_realizations
   plot_geometry
   property_map
   reset
   reset_timers
   run
   run_in_memory
   sample_external_source
   simulation_finalize
   simulation_init
   source_bank
   statepoint_load
   statepoint_write

Classes
-------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   AzimuthalFilter
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
   Mesh
   MeshFilter
   MeshBornFilter
   MeshSurfaceFilter
   MuFilter
   Nuclide
   ParticleFilter
   PolarFilter
   RectilinearMesh
   RegularMesh
   SpatialLegendreFilter
   SphericalHarmonicsFilter
   SphericalMesh
   SurfaceFilter
   Tally
   UniverseFilter
   UnstructuredMesh
   WeightWindows
   ZernikeFilter
   ZernikeRadialFilter

Data
----

.. data:: cells

   Mapping of cell ID to :class:`openmc.lib.Cell` instances.

   :type: dict

.. data:: filters

   Mapping of filter ID to :class:`openmc.lib.Filter` instances.

   :type: dict

.. data:: materials

   Mapping of material ID to :class:`openmc.lib.Material` instances.

   :type: dict

.. data:: meshes

   Mapping of mesh ID to :class:`openmc.lib.Mesh` instances.

   :type: dict

.. data:: nuclides

   Mapping of nuclide name to :class:`openmc.lib.Nuclide` instances.

   :type: dict

.. data:: tallies

   Mapping of tally ID to :class:`openmc.lib.Tally` instances.

   :type: dict

.. data:: weight_windows

   Mapping of weight window ID to :class:`openmc.lib.WeightWindows` instances.

   :type: dict
