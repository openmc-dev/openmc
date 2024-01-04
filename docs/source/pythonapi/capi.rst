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
   export_properties
   export_weight_windows
   import_weight_windows
   finalize
   find_cell
   find_material
   global_bounding_box
   global_tallies
   hard_reset
   id_map
   import_properties
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
   simulation_init
   simulation_finalize
   source_bank
   statepoint_write

Classes
-------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   Cell
   CylindricalMesh
   EnergyFilter
   MaterialFilter
   Material
   MeshFilter
   MeshSurfaceFilter
   Nuclide
   RectilinearMesh
   RegularMesh
   SphericalMesh
   Tally
