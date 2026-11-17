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
   feature_enabled
   finalize
   find_cell
   find_material
   dagmc_universe_cell_ids
   global_bounding_box
   global_tallies
   get_emission_spectrum
   get_emission_distribution
   get_angular_distribution
   sample_angular_distribution
   get_reaction_products
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
   run_random_ray
   sample_external_source
   slice_data
   slice_data_overlap_count
   slice_data_overlap_info
   set_emission_spectrum
   simulation_finalize
   simulation_init
   source_bank
   statepoint_load
   statepoint_write

Runtime secondary-distribution inspection
-----------------------------------------

:func:`get_reaction_products` discovers the native reaction products of a loaded
nuclide. :func:`get_emission_distribution` copies a selected native angle-energy
law into a tagged dictionary. These functions inspect the C++ objects used by
transport; they do not reopen nuclear-data files or reconstruct distributions
through the offline Python data API::

    # After Model.init_lib(), with Li6 present in the model:
    products = openmc.lib.get_reaction_products('Li6')
    for product in products:
        for distribution in product['distributions']:
            print(product['mt'], product['product_index'], product['particle'],
                  distribution['distribution_index'],
                  distribution['angle_energy_type'], distribution['energy_type'])

    law = openmc.lib.get_emission_distribution('Li6', 51)
    print(law['energy'])  # Native level-inelastic parameters, if present

Discovery returns one dictionary per product, including ``nuclide``,
``nuclide_index``, ``mt``, ``product_index``, ``particle``, ``emission_mode``,
``reference_frame`` (the reaction's stored frame), ``parent_nuclide_index``,
``yield``, ``applicability``, and ``distributions``. Product and distribution
indices refer to the loaded C++ vectors, including products substituted by D1S.
An empty applicability list means that no selection functions were stored;
the reader does not manufacture a constant function. Yield and applicability
are inspected independently of distribution shape and cannot be changed here.

Each distribution descriptor contains ``distribution_index``,
``angle_energy_type``, ``energy_type`` (None when not applicable),
``read_supported``, ``mutation_supported``, and ``mutation_reason``.
Whole-distribution mutation is not implemented. The separate
``energy_mutation_supported`` flag identifies existing ContinuousTabular
energy-only editing using ``set_emission_spectrum``. Capabilities describe API
support, not permission to edit in the current simulation state. Independent
angular editing is reported separately by ``angular_mutation_supported`` and
``angular_mutation_reason``; it preserves native grids and interpolation.
Use ``include_data=True`` to include each complete law under ``data``;
the default avoids copying large joint tables.

The schema (``schema_version=1``) preserves the native hierarchy:

* ``uncorrelated`` contains ``angle`` and ``energy``. An absent angle is tagged
  ``isotropic``; an absent energy object is None, as used by elastic kinematics.
* ``tabular-angle`` contains an incident grid and native cosine tables.
* ``continuous-tabular`` contains incident energy, breakpoints, interpolation,
  and outgoing tables with energy, PDF, CDF, interpolation, and discrete count.
* ``discrete-photon`` contains primary flag, energy, and atomic weight ratio.
* ``level-inelastic`` contains threshold and mass ratio.
* ``maxwell`` and ``evaporation`` contain tabulated ``theta`` and ``u``.
* ``watt`` contains separate tabulated ``a`` and ``b`` functions and ``u``.
* ``correlated`` retains conditional angle tables inside each outgoing-energy
  table. These are not independent angular marginals.
* ``kalbach-mann`` retains outgoing PDF/CDF and corresponding ``r`` and ``a``.
* ``nbody`` contains particle count, total mass, atomic weight ratio, Q, and
  an explicit isotropic ``angle``.
* ``decay-photon`` contains a native generic spectrum, isotropic ``angle``,
  and ownership provenance.
  Its storage is borrowed from a chain nuclide and can be shared by several
  products and consumers. This interface does not change that ownership.

Tabulated functions contain ``x``, ``y``, ``breakpoints``, and ``interpolation``;
polynomials contain coefficients in ascending power order. Public breakpoints
use one-based ENDF/HDF5 endpoints. Interpolation values are the native enums:
the existing loader reduces codes such as 12 and 22 to 2, discarding their
two-dimensional qualifier. ``sampling_interpolation`` separately describes
the behavior of the current sampler. The reader cannot recover metadata
discarded during loading.

Numerical arrays are independent NumPy copies. Native PDFs and CDFs are neither
normalized nor rebuilt by inspection. In particular, angular ``Tabular``
objects may already have been normalized by the native constructor. Generic
D1S discrete/mixture spectra expose the native alias tables and importance
weights instead of inventing original probabilities that are no longer stored.
All energies are in eV; angular abscissae are scattering cosines.

Reads require loaded nuclear data but do not require ``simulation_init``.
Serialize library calls and perform inspection outside batches. Snapshots
remain valid after finalization and editing them has no effect on transport.
No new setter is introduced. Thermal-scattering laws and external source
distributions are outside this ReactionProduct interface.

The C ABI uses ``openmc_secondary_snapshot_create``,
``openmc_secondary_snapshot_node``, and ``openmc_secondary_snapshot_free``.
Snapshots own a tagged tree of scalars, arrays, objects, and lists. Node views
are const pointers into the snapshot, never into live transport objects, and
remain valid only until that snapshot is freed. The Python wrapper copies all
values and frees the snapshot even if decoding fails.

Run ``examples/discover_secondary_distributions.py`` to print representative
MT/product/distribution addresses for the configured library's light,
structural, and uranium nuclides. The script loads and reads data without
running transport or modifying cross sections.

Runtime emission spectra
------------------------

:func:`get_emission_spectrum` and :func:`set_emission_spectrum` provide in-memory
access to a loaded nuclide's secondary energy spectrum. Currently, only
uncorrelated continuous-tabular spectra are supported. Select a spectrum using
the nuclide index, reaction MT number, and zero-based product and distribution
indices. For example, the NNDC HDF5 library represents U235 MT=18's prompt neutron
product as a continuous-tabular spectrum::

    import numpy as np
    import openmc.lib

    # After loading a model with U235 through Model.init_lib():
    openmc.lib.simulation_init()
    index = openmc.lib.nuclides['U235']._index
    spectrum = openmc.lib.get_emission_spectrum(index, 18)
    for table in spectrum['tables']:
        table['pdf'] *= np.exp(-table['energy'] / 3.e5)
    openmc.lib.set_emission_spectrum(index, 18, spectrum)
    changed = openmc.lib.get_emission_spectrum(index, 18)
    openmc.lib.next_batch()

The getter returns independent arrays. Edits take effect only when passed to
the setter, which validates and replaces the entire spectrum atomically. The
incident-energy grid and table count must match native data exactly. Outgoing
grids may change. Probabilities are normalized and CDFs are rebuilt; any CDF in
the supplied dictionary is ignored. For mixed tables, the first ``n_discrete``
entries are discrete line probabilities, followed by continuous densities in
1/eV. Histogram and linear-linear interpolation are supported. Discrete-only
tables cannot be adjacent to continuous tables when incident-energy sampling
uses unit-base interpolation.

Mutation requires :func:`simulation_init` and is only allowed between batches,
before :func:`simulation_finalize`. Serialize library calls. Under MPI, make
the same edit on every rank. Mutations affect subsequent sampling, not particles
already stored in the source bank, and are neither written to nuclear-data files
nor saved as edits in statepoint files.

The repository example ``examples/runtime_emission_spectrum.py`` compares two
runs with the same seed and an identical first batch. It saves original and
modified native PDFs, plots subsequent fission source-bank energies, and writes
per-batch mean energies and raw source energies. Run it with the NNDC data
configured through ``OPENMC_CROSS_SECTIONS``::

    OMP_NUM_THREADS=2 python examples/runtime_emission_spectrum.py --output demo

Add ``--event`` to exercise event-based transport. This deliberate spectrum
perturbation demonstrates runtime sensitivity; it does not provide a new
physical nuclear-data evaluation.

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
   ParentNuclideFilter
   ParticleFilter
   ParticleProductionFilter
   PolarFilter
   ReactionFilter
   RectilinearMesh
   RegularMesh
   SpatialLegendreFilter
   SphericalHarmonicsFilter
   SphericalMesh
   SolidRayTracePlot
   SurfaceFilter
   Tally
   TemperatureField
   TemporarySession
   UniverseFilter
   UnstructuredMesh
   WeightFilter
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

.. data:: plots

   Mapping of plot ID to :class:`openmc.lib.SolidRayTracePlot` instances.

   :type: dict

.. data:: nuclides

   Mapping of nuclide name to :class:`openmc.lib.Nuclide` instances.

   :type: dict

.. data:: tallies

   Mapping of tally ID to :class:`openmc.lib.Tally` instances.

   :type: dict

.. data:: temperature_field

   Temperatures in [K] of each element of the temperature field mesh, if one
   was specified. Empty when no temperature field is in use.

   :type: :class:`openmc.lib.TemperatureField`

.. data:: weight_windows

   Mapping of weight window ID to :class:`openmc.lib.WeightWindows` instances.

   :type: dict

Runtime angular inspection and validation
-----------------------------------------

``get_emission_distribution`` remains the complete native-law reader.
``get_angular_distribution`` returns its independent angular component for
uncorrelated laws and an explicit ``isotropic`` tag for N-body/decay photons.
It raises ValueError for correlated and Kalbach-Mann laws, directing callers to
the complete reader. No angular marginal is constructed and no setter is added.

The independent angular schema is::

    {"type": "isotropic"}
    {"type": "tabular-angle", "incident_energy": array,
     "sampling_interpolation": "linear-stochastic",
     "tables": [{"type": "tabular", "x": mu_array, "pdf": array,
                 "cdf": array, "interpolation": int, "integral": float}, ...]}

``x`` is mu = cos(theta) in the product's native reference frame. ``integral``
is the terminal CDF before normalization by the C++ Tabular constructor; the
returned PDF and CDF have already been divided by it. Inspection makes no
further changes. Both PDF and CDF control sampling: the CDF chooses an interval
and the PDF determines the inverse within that interval. Independently verify
their consistency; evaluated data may contain small rounding differences.

AngleDistribution stores no incident breakpoints/interpolation regions. Its
sampler linearly mixes the two adjacent incident tables by stochastic selection,
clamping to the first/last table outside the grid. Cosine interpolation uses the
native table enum: 1 is histogram and 2 is linear-linear in nuclear-data files.
Python Legendre/Uniform laws are converted to tables during HDF5 export; their
original coefficients are absent from the runtime object and cannot be recovered.
An isotropic *table* remains tabular in the snapshot; only an absent angle law
is tagged isotropic for an UncorrelatedAngleEnergy object.

Correlated snapshots retain ``tables[i]["angle"][j]`` alongside outgoing-energy
``energy[j]``, PDF, CDF, interpolation, and discrete-energy count. Each angular
entry uses the same Tabular schema above. The sampler selects the lower angular
table for histogram outgoing-energy interpolation; otherwise it chooses the
nearest outgoing CDF endpoint. It does not interpolate angular PDFs in outgoing
energy. Kalbach-Mann snapshots retain corresponding ``r[j]`` and ``a[j]``;
these are histogram-selected or linearly interpolated at the sampled outgoing
energy before unit-base mapping. Their conditional angular density is
``a/(2*sinh(a)) * (cosh(a*mu) + r*sinh(a*mu))``. No outgoing-energy integration
occurs in either reader.

``sample_angular_distribution`` samples independent laws through the same native
AngleDistribution sampler used by elastic transport, using caller-owned RNG
state. It also samples explicit isotropic laws. This diagnostic does not run a
collision, choose product applicability, apply thermal target motion or convert
the product frame to laboratory directions. Calls require loaded data, must be
serialized, and are rejected during batches. Correlated laws are rejected.

Run ``examples/runtime_angular_distribution.py --inventory`` for representative
native laws, checks of every angular table in the requested nuclides, and Li6
elastic sampling validation. Use ``--mt 24`` for the non-elastic Li6 example.
For joint laws use, for example, ``--nuclide Be9 --mt 16 --incident-table 12
--outgoing-point 10`` or ``--nuclide C12 --mt 5``. These plot explicit native
conditional points, not marginals. Plots label the product frame, integrate
histograms and linear PDFs using their native laws, and transform to theta with
``p_theta(theta) = p_mu(cos(theta))*sin(theta)``. The theta axis is shown in
degrees while its density is explicitly per radian.


Runtime sensitivity perturbations
--------------------------------

.. autofunction:: openmc.lib.get_reaction_info
.. autofunction:: openmc.lib.set_angular_distribution
.. autofunction:: openmc.lib.sample_emission_energy
.. autofunction:: openmc.lib.perturb_xs
.. autofunction:: openmc.lib.perturb_tabular
.. autofunction:: openmc.lib.nuclear_data_perturbation

The combined context validates candidates before application and captures
owned native rollback copies. Original PDF/CDF values are restored without
reconstruction, avoiding cumulative drift. Native setter/rebuild failures roll
back earlier changes under serialized library access. This is rollback
atomicity, not concurrency isolation. Library finalization invalidates backup
handles and must not occur inside a perturbation context. Mutations require
simulation initialization and are rejected during active batches.

XS candidates are full arrays on the existing nuclide grid, with the original
threshold and zeros below it preserved. Only the selected reaction changes;
derived total, absorption and photon-production arrays are rebuilt explicitly.
Redundant XS and nuclides with URR, multipole, resonance-scattering or fission
representations are excluded from the combined context because smooth-array
changes alone do not update their dependent data consistently. The older XS
setter remains available with its existing explicit-rebuild contract.

Energy shape changes support native uncorrelated ContinuousTabular laws.
Angular shape changes support existing independent tabular-angle laws. The
context preserves grids, interpolation, discrete-energy counts/probabilities,
reaction/product identity, yield, applicability and reference frames. Angular
Tabular objects have no discrete angular masses. Histogram and linear-linear
PDFs are integrated by their respective native rules. Exponential shape tilts
are applied at native grid points and the resulting native interpolant is
normalized; the exponential is not evaluated continuously between points.
An existing repeated outgoing-energy coordinate is retained only when the
native grid, interpolation and discrete count are unchanged.

Isotropic angular laws have no parameter to change and are not implicitly
converted to tables. LevelInelastic and DiscretePhoton energies remain read-only
because their parameters encode kinematics/level information. MaxwellEnergy,
Evaporation and WattEnergy remain read-only pending parameter-function and
energy-restriction validation adapters. NBodyPhaseSpace, CorrelatedAngleEnergy
and KalbachMann remain read-only; independent PDF changes would not preserve
their native joint structure. In particular, conditional angular tables and
Kalbach r/a parameters require separate joint-law adapters.

``examples/fusion_li6_runtime_uq.py`` demonstrates all three changes in one
loaded library session and generates plots, a controlled comparison CSV,
a small Latin-hypercube training CSV, metadata JSON and a repeated-run report.
All angular/energy data come through runtime APIs. No evaluated nuclear-data
files are rewritten or reprocessed. The coefficients describe controlled
sensitivity experiments, not an evaluated covariance uncertainty model.

For matched repeatability, finish one run, call ``hard_reset()``, set the desired
seed and enter the next context. ``run()`` finalizes simulation state but retains
loaded nuclear data. Extract tally copies before initializing another run.
The example uses one transport thread for bitwise tally comparisons; parallel
reductions can differ in the last floating-point bits.
