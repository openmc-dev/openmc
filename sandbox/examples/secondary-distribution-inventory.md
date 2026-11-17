# Native secondary-distribution discovery: ENDF/B-VIII.0

Generated using `sandbox/examples/discover_secondary_distributions.py` against
`/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml`.

Each row groups distributions with the same particle, emission mode, stored
reaction frame, and native law types. The example address uses runtime product
indices (which need not match numerical HDF5 product suffixes). Every reported
example was retrieved through `get_emission_distribution` from a loaded C++
object. Occurrences count distributions, not distinct reactions. This is a
read-only discovery run; no transport or cross-section modification is performed.

The frame is the reaction's stored flag. Transport treats elastic kinematics,
photons, and fission separately; it does not apply the inelastic-neutron frame
transformation to every product. A dash in the energy-law column means no
separate EnergyDistribution object, including elastic and joint-law cases.

Whole-law mutation is unavailable in Step 1. “Energy only” denotes the existing
ContinuousTabular setter. Watt is supported by the reader and tested using the
native HDF5 loader, but does not occur in the configured library scan. D1S is
verified separately with `tests/chain_ni.xml` and remains read-only.

| Nuclide | Example MT/product/distribution | Particle | Mode | Frame | AngleEnergy | Energy law | Occurrences | Read | Write |
|---|---|---|---|---|---|---|---:|---|---|
| H1 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| H1 | 102/0/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 1 | yes | energy only (legacy API) |
| H2 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| H2 | 16/0/0 | neutron | prompt | center-of-mass | nbody | — | 1 | yes | no |
| H2 | 102/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 1 | yes | no |
| He3 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| He4 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Li6 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Li6 | 24/0/0 | neutron | prompt | laboratory | uncorrelated | continuous-tabular | 1 | yes | energy only (legacy API) |
| Li6 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 31 | yes | no |
| Li6 | 57/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 1 | yes | no |
| Li6 | 102/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 3 | yes | no |
| Li7 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Li7 | 16/0/0 | neutron | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Li7 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 32 | yes | no |
| Li7 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 1 | yes | no |
| Li7 | 102/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 3 | yes | no |
| Be9 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Be9 | 16/0/0 | neutron | prompt | laboratory | correlated | — | 1 | yes | no |
| Be9 | 102/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 11 | yes | no |
| Be9 | 701/0/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 1 | yes | energy only (legacy API) |
| C12 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| C12 | 5/0/0 | neutron | prompt | center-of-mass | kalbach-mann | — | 1 | yes | no |
| C12 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 1 | yes | energy only (legacy API) |
| C12 | 28/0/0 | neutron | prompt | laboratory | uncorrelated | evaporation | 2 | yes | no |
| C12 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 12 | yes | no |
| C12 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 1 | yes | no |
| C12 | 102/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 3 | yes | no |
| O16 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| O16 | 4/0/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 23 | yes | no |
| O16 | 5/0/0 | neutron | prompt | center-of-mass | kalbach-mann | — | 10 | yes | no |
| O16 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 8 | yes | energy only (legacy API) |
| O16 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 7 | yes | no |
| O16 | 102/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 97 | yes | no |
| O16 | 108/0/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Al27 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Al27 | 5/0/0 | neutron | prompt | center-of-mass | kalbach-mann | — | 8 | yes | no |
| Al27 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 47 | yes | energy only (legacy API) |
| Al27 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 39 | yes | no |
| Al27 | 102/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 290 | yes | no |
| Al27 | 102/213/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 76 | yes | energy only (legacy API) |
| Fe54 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Fe54 | 5/0/0 | neutron | prompt | center-of-mass | correlated | — | 3 | yes | no |
| Fe54 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe54 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 10 | yes | no |
| Fe54 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 31 | yes | no |
| Fe54 | 102/0/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe54 | 601/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 341 | yes | no |
| Fe56 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Fe56 | 5/0/0 | neutron | prompt | center-of-mass | correlated | — | 3 | yes | no |
| Fe56 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe56 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 39 | yes | no |
| Fe56 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 367 | yes | no |
| Fe56 | 102/0/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe56 | 601/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 159 | yes | no |
| Fe57 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Fe57 | 5/0/0 | neutron | prompt | center-of-mass | correlated | — | 3 | yes | no |
| Fe57 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe57 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 14 | yes | no |
| Fe57 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 153 | yes | no |
| Fe57 | 102/0/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe57 | 601/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 189 | yes | no |
| Fe58 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Fe58 | 5/0/0 | neutron | prompt | center-of-mass | correlated | — | 3 | yes | no |
| Fe58 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe58 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 39 | yes | no |
| Fe58 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 269 | yes | no |
| Fe58 | 102/0/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| Fe58 | 601/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 93 | yes | no |
| U235 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| U235 | 5/0/0 | neutron | prompt | center-of-mass | correlated | — | 4 | yes | no |
| U235 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 4 | yes | energy only (legacy API) |
| U235 | 18/0/0 | neutron | prompt | laboratory | uncorrelated | continuous-tabular | 1 | yes | energy only (legacy API) |
| U235 | 18/1/0 | neutron | delayed | laboratory | uncorrelated | continuous-tabular | 6 | yes | energy only (legacy API) |
| U235 | 18/7/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| U235 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 39 | yes | no |
| U235 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 445 | yes | no |
| U235 | 801/0/0 | photon | prompt | laboratory | uncorrelated | discrete-photon | 131 | yes | no |
| U238 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| U238 | 5/0/0 | neutron | prompt | center-of-mass | correlated | — | 4 | yes | no |
| U238 | 5/1/0 | photon | prompt | center-of-mass | uncorrelated | continuous-tabular | 4 | yes | energy only (legacy API) |
| U238 | 18/0/0 | neutron | prompt | laboratory | uncorrelated | continuous-tabular | 1 | yes | energy only (legacy API) |
| U238 | 18/1/0 | neutron | delayed | laboratory | uncorrelated | continuous-tabular | 6 | yes | energy only (legacy API) |
| U238 | 18/7/0 | photon | prompt | laboratory | uncorrelated | continuous-tabular | 3 | yes | energy only (legacy API) |
| U238 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 39 | yes | no |
| U238 | 51/1/0 | photon | prompt | center-of-mass | uncorrelated | discrete-photon | 351 | yes | no |
| Ra223 | 2/0/0 | neutron | prompt | center-of-mass | uncorrelated | — | 1 | yes | no |
| Ra223 | 16/0/0 | neutron | prompt | laboratory | uncorrelated | evaporation | 10 | yes | no |
| Ra223 | 18/0/0 | neutron | total | laboratory | uncorrelated | maxwell | 1 | yes | no |
| Ra223 | 51/0/0 | neutron | prompt | center-of-mass | uncorrelated | level-inelastic | 14 | yes | no |

## Step 1 implementation and verification

New public Python functions: `get_reaction_products` and
`get_emission_distribution`. No new mutation API was added. Existing
`get_emission_spectrum` / `set_emission_spectrum` remain available.

Files added or updated specifically for Step 1 (the checkout already contained
uncommitted ContinuousTabular implementation, tests, and documentation):

- `src/secondary_data.cpp`, `include/openmc/secondary_data.h`: common resolver,
  native inspection adapters, and snapshot implementation.
- `include/openmc/capi.h`: tagged, owned snapshot C ABI.
- `src/nuclide.cpp`: existing spectrum API now uses the common resolver.
- `include/openmc/{chain,distribution,distribution_angle,distribution_energy,endf,secondary_correlated,secondary_kalbach,secondary_nbody,secondary_uncorrelated}.h`:
  internal adapter friendship; no changes to sampling or ownership semantics.
- `openmc/lib/secondary.py`, `openmc/lib/__init__.py`: Python decoding and exports.
- `CMakeLists.txt`, `tests/cpp_unit_tests/CMakeLists.txt`: source/test registration.
- `tests/cpp_unit_tests/test_secondary_data.cpp`: HDF5 loader-based coverage for
  all nine non-D1S law combinations, absent energy, generic D1S spectra, snapshot
  independence, unchanged native samples, C ABI errors, and in-batch rejection.
- `tests/unit_tests/test_lib_secondary.py`: real-data discovery and complete
  payload comparisons, copied arrays, invalid selections, and visibility of
  existing native ContinuousTabular mutation through the new getter.
- `tests/unit_tests/test_lib_secondary_chain.py`: D1S discovery, capabilities,
  lifecycle, and copied data surviving native finalization.
- `tests/unit_tests/test_lib_secondary_watt.py`: native Watt loading and Python
  readback from an explicitly synthetic HDF5 fixture.
- `docs/source/pythonapi/capi.rst`: schema, lifecycle, capabilities, native
  interpolation behavior, C ABI lifetime, and D1S ownership documentation.
- `sandbox/examples/discover_secondary_distributions.py` and this report.

Verification:

- `cmake --build build -j 4`: passed; RelWithDebInfo, strict floating point,
  MPI disabled.
- CTest `test_secondary_data`, `test_emission_spectrum`, `test_distribution`:
  all three executables passed. New C++ tests: 96 assertions in four test cases.
- Python new read tests plus existing `test_lib_emission_spectrum.py`:
  79 passed with `OMP_NUM_THREADS=2` and installed ENDF/B-VIII.0 data.
- Discovery script: passed for all 17 selected nuclides, including retrieval
  of the full native payload for each of the 90 reported representative rows.
- Existing `runtime_emission_spectrum.py`, 1,000 particles and four batches:
  completed successfully and confirmed an identical initial batch. This is
  legacy mutation compatibility coverage, not new mutation support.
- `git diff --check` and Python compilation: passed.

The installed data differs from the NNDC regression reference library. No
claim is made that the full NNDC regression suite or MPI configuration was
validated. Step 2 has not been implemented.
