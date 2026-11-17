# Read-only runtime angular inspection

The existing `get_reaction_products` / `get_emission_distribution` interface and
`SecondaryDataReader` snapshots are retained. No angular mutation is implemented.
`get_angular_distribution` provides an independent-law convenience reader;
`sample_angular_distribution` provides native angular sampling with caller-owned
RNG state. Neither changes live nuclear data or transport RNG streams.

## Native class inventory

| ReactionProduct law | Actual angular representation | Independent in native frame? |
|---|---|---|
| `UncorrelatedAngleEnergy` | `AngleDistribution` containing `vector<unique_ptr<Tabular>>`, or empty-angle isotropic fallback | Yes |
| `CorrelatedAngleEnergy` | `CorrTable::angle[j]`, native `Tabular` conditional on outgoing energy selection | No |
| `KalbachMann` | `KMTable::r[j]`, `a[j]`, analytic conditional density | No |
| `NBodyPhaseSpace` | Uniform cosine sampled independently of its phase-space energy | Yes, isotropic |
| `DecayPhotonAngleEnergy` (D1S chain products) | `Uniform(-1,1)` independent of borrowed photon energy spectrum | Yes, isotropic |

These are all laws constructed by `ReactionProduct` in
[src/reaction_product.cpp](../src/reaction_product.cpp). There is no native
reaction-angular `Isotropic`, `Legendre`, `Discrete`, or `Mixture` subclass.
The general C++ Distribution hierarchy also includes Discrete, Uniform,
PowerLaw, Maxwell, Watt, Normal, Equiprobable, Mixture, and DecaySpectrum, but
ReactionProduct angular readers do not instantiate them. Generic distributions
in decay-photon snapshots describe **energy**, not angle.

Python `openmc.data` can hold Uniform, Legendre, Discrete, and Mixture objects.
`AngleDistribution.to_hdf5` converts non-Tabular angular laws to tables.
`CorrelatedAngleEnergy.to_hdf5` also converts analytic angular laws to tables,
retaining a special zero interpolation code for discrete thermal angles.
Original Legendre coefficients do not survive into the runtime law.
No coefficients or analytic types are invented from those tables.

Thermal AngleEnergy classes (CoherentElasticAE, IncoherentElasticAE,
IncoherentElasticAEDiscrete, IncoherentInelasticAEDiscrete,
IncoherentInelasticAE, MixedElasticAE) are attached to thermal-scattering data,
not ReactionProducts; they are outside this API. External-source polar/azimuthal
distributions are also outside it.

## Exact sampler provenance

| HDF5 field | C++ storage | Snapshot | Sampling use |
|---|---|---|---|
| `angle/energy` | `AngleDistribution::energy_` | `angle.incident_energy` | `get_energy_index` computes linear fraction r; select upper table iff r > PRN |
| `angle/mu` row 0, attributes `offsets` | `Tabular::x_` | `angle.tables[i].x` | Inverse CDF starts at x[i]; final angular sample is clamped to [-1,1] |
| `angle/mu` row 1 | `Tabular::p_` after division by terminal CDF | `angle.tables[i].pdf` | Histogram: x[i]+(u-c[i])/p[i]; linear: quadratic inverse using endpoint PDF slope |
| `angle/mu` row 2 | `Tabular::c_` after division by terminal CDF | `angle.tables[i].cdf` | `lower_bound(c.begin()+1,c.end(),u)` chooses interval |
| `angle/mu` attribute `interpolation[i]` | `Tabular::interp_` | table `interpolation` | Selects histogram (1) or linear-linear (2) inverse |
| `angle/mu` last CDF | `Tabular::integral_` | table `integral` | Loader normalization divisor; not a second normalization to apply to returned PDF |
| absent `angle` group | empty `AngleDistribution` | `angle.type=isotropic` | `uniform_distribution(-1,1,seed)` |
| correlated `energy`, interpolation attributes | `energy_`, `breakpoints_`, `interpolation_` | law `incident_energy`, `breakpoints`, `interpolation` | Actual incident selection is linear stochastic, regardless of stored region metadata |
| correlated `energy_out` rows 0–2 | `CorrTable::e_out,p,c` | `tables[i].energy,pdf,cdf` | Select/invert outgoing CDF, then unit-base map continuous outgoing energy |
| correlated `energy_out` rows 3–4 and `mu` rows 0–2 | conditional `CorrTable::angle[j]` Tabular objects | `tables[i].angle[j]` using x/pdf/cdf/interpolation/integral | Histogram outgoing interpolation selects angle[k]; otherwise choose k or k+1 by distance of sampled u to outgoing CDF endpoints, then call Tabular::sample |
| correlated `energy_out` attributes `n_discrete_lines`, `interpolation`, `offsets` | outgoing table `n_discrete`, `interpolation` and array boundaries | outgoing table metadata | Separates discrete-energy selection from continuous inverse/unit-base mapping |
| Kalbach `distribution` rows 0–2 / attributes | `KMTable::e_out,p,c,n_discrete,interpolation` | outgoing tables | Same energy selection/inversion and unit-base mapping pattern |
| Kalbach `distribution` rows 3–4 | `KMTable::r,a` | outgoing table `r,a` | Histogram-selected or linearly interpolated at sampled outgoing energy before unit-base mapping |
| N-body attributes `n_particles,total_mass,atomic_weight_ratio,q_value` | `n_bodies_,mass_ratio_,A_,Q_` | same attribute names, plus isotropic angle tag | Parameters determine outgoing energy; angle is separately uniform |
| chain photon spectrum | borrowed `photon_energy_` | `spectrum`, `ownership`, isotropic angle tag | `DecayPhotonAngleEnergy::sample` samples spectrum, then `Uniform(-1,1)` |

The independent sampler path is
[src/secondary_uncorrelated.cpp](../src/secondary_uncorrelated.cpp) →
[src/distribution_angle.cpp](../src/distribution_angle.cpp) →
`Distribution::sample` → `Tabular::sample_unbiased` in
[src/distribution.cpp](../src/distribution.cpp).
`Tabular::init` copies HDF5 CDF values instead of reconstructing them, then
divides PDF and CDF by their terminal CDF. Inspection copies these live values.
CDF/PDF rounding inconsistencies therefore remain observable.

For correlated laws the path is
[src/secondary_correlated.cpp](../src/secondary_correlated.cpp)
`sample` → `sample_dist` → selected angular `Tabular::sample`.
The closest endpoint is measured in outgoing **CDF** distance, not energy
distance. The exposed table at an outgoing point is a native conditional table,
not an interpolated conditional density for arbitrary laboratory E_out.

For Kalbach the path is
[src/secondary_kalbach.cpp](../src/secondary_kalbach.cpp)
`sample` → `sample_params`. Given sampled r,a, one branch samples the symmetric
cosh density via `asinh(U(-1,1)*sinh(a))/a`; the other samples the forward
exponential using `log(u*exp(a)+(1-u)*exp(-a))/a`. Its density is
`a/(2*sinh(a))*(cosh(a*mu)+r*sinh(a*mu))`.
The stable plotting expression is algebraically equivalent; a=0 has a uniform
analytic limit. The native sampler itself divides by a and does not special-case
zero, so inspecting a zero parameter is not evidence of safe native sampling.

N-body angle sampling is explicit in
[src/secondary_nbody.cpp](../src/secondary_nbody.cpp); decay-photon angle
sampling is explicit in [src/chain.cpp](../src/chain.cpp).

Transport integration is in [src/physics.cpp](../src/physics.cpp):
`elastic_scatter` directly samples `products_[0].distribution_[0]->angle()`;
`inelastic_scatter` calls `ReactionProduct::sample` → `sample_dist` (including
applicability selection) → the selected AngleEnergy law. It then converts CM
energy/cosine to LAB where required and calls `rotate_angle`. Fission and
secondary-photon production also call product sampling and rotate directions.
Thermal target motion and material P0 scattering overrides can change final
transport directions; the reported native law precedes those operations.
The diagnostic samples a specified distribution and does not sample its
applicability or model an entire collision.

## Python schemas

```python
# Full uncorrelated law
{"type": "uncorrelated", "angle": angle, "energy": energy_or_None}
# Independent angle
{"type": "isotropic"}
{"type": "tabular-angle", "incident_energy": np_array,
 "sampling_interpolation": "linear-stochastic", "tables": [angular_table, ...]}
# Native angular table; x means mu=cos(theta)
{"type": "tabular", "x": np_array, "pdf": np_array, "cdf": np_array,
 "interpolation": int, "integral": float}
# Correlated law
{"type": "correlated", "incident_energy": np_array,
 "breakpoints": int_array, "interpolation": int_array,
 "sampling_interpolation": "linear-stochastic-unit-base",
 "tables": [{"energy": np_array, "pdf": np_array, "cdf": np_array,
             "interpolation": int, "n_discrete": int,
             "angle": [angular_table, ...]}, ...]}
# Kalbach law: same incident/outgoing structure, with r and a instead of angle
{"type": "kalbach-mann", "incident_energy": np_array,
 "breakpoints": int_array, "interpolation": int_array,
 "sampling_interpolation": "linear-stochastic-unit-base",
 "tables": [{"energy": np_array, "pdf": np_array, "cdf": np_array,
             "interpolation": int, "n_discrete": int,
             "r": np_array, "a": np_array}, ...]}
# N-body: existing native parameters, now also angle={"type":"isotropic"}
# Decay-photon: existing borrowed spectrum/ownership, now also isotropic angle
```

Outgoing tables also retain any pre-existing energy snapshot fields. Incident
breakpoints are one-based. The native loader reduces interpolation codes modulo
10, so discarded two-dimensional qualifiers cannot be recovered. Uncorrelated
AngleDistribution stores **no** incident breakpoints or interpolation regions;
none are fabricated. Arrays are independent copies, valid after finalization.
All energies are eV. Frames come from product discovery. Correlated/Kalbach
convenience reads raise a clear error directing callers to the complete reader;
no independent angular marginal exists in this API.

## ENDF/B-VIII.0 examples and complete Li6 neutron inventory

Configured library: `/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml`.
The following rows were read from loaded C++ snapshots. Each requested nuclide
has representative rows; all Li6 neutron laws are included. Incident-grid ranges
refer to angular grids, not cross-section grids or reaction thresholds.

| Nuclide | MT | Product | Distribution | Particle | Frame | Law | Angular law | Incident grid (eV) |
|---|---|---|---|---|---|---|---|---|
| H1 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 153 points, 1e-05–2e+07 |
| H1 | 102 | 0 | 0 | photon | laboratory | uncorrelated | tabular-angle | 2 points, 1e-05–2e+07 |
| H2 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 348 points, 1e-05–1.5e+08 |
| H2 | 16 | 0 | 0 | neutron | center-of-mass | nbody | isotropic | none |
| H2 | 102 | 0 | 0 | photon | laboratory | uncorrelated | tabular-angle | 2 points, 1e-05–1.5e+08 |
| Li6 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 216 points, 1e-05–2e+07 |
| Li6 | 24 | 0 | 0 | neutron | laboratory | uncorrelated | tabular-angle | 20 points, 4.318e+06–2e+07 |
| Li6 | 51 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.7515e+06–2e+07 |
| Li6 | 52 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 2.3354e+06–2e+07 |
| Li6 | 53 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 17 points, 2.5514e+06–2e+07 |
| Li6 | 54 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 2.9192e+06–2e+07 |
| Li6 | 55 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 3.5031e+06–2e+07 |
| Li6 | 56 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 4.0869e+06–2e+07 |
| Li6 | 57 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 4.159e+06–2e+07 |
| Li6 | 57 | 1 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 4.159e+06–2e+07 |
| Li6 | 58 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 4.6708e+06–2e+07 |
| Li6 | 59 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 5.2546e+06–2e+07 |
| Li6 | 60 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 5.8384e+06–2e+07 |
| Li6 | 61 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 6.4223e+06–2e+07 |
| Li6 | 62 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 7.0061e+06–2e+07 |
| Li6 | 63 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 7.59e+06–2e+07 |
| Li6 | 64 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 8.1738e+06–2e+07 |
| Li6 | 65 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 8.7577e+06–2e+07 |
| Li6 | 66 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 9.3415e+06–2e+07 |
| Li6 | 67 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 9.9253e+06–2e+07 |
| Li6 | 68 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.0509e+07–2e+07 |
| Li6 | 69 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.1093e+07–2e+07 |
| Li6 | 70 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.1677e+07–2e+07 |
| Li6 | 71 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.2261e+07–2e+07 |
| Li6 | 72 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.2845e+07–2e+07 |
| Li6 | 73 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.3428e+07–2e+07 |
| Li6 | 74 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.4012e+07–2e+07 |
| Li6 | 75 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.4596e+07–2e+07 |
| Li6 | 76 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.518e+07–2e+07 |
| Li6 | 77 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.5764e+07–2e+07 |
| Li6 | 78 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.6348e+07–2e+07 |
| Li6 | 79 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.6931e+07–2e+07 |
| Li6 | 80 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.7515e+07–2e+07 |
| Li6 | 81 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 2 points, 1.8099e+07–2e+07 |
| Li7 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 75 points, 1e-05–2e+07 |
| Li7 | 51 | 1 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 546277–2e+07 |
| Be9 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 94 points, 1e-05–2e+07 |
| Be9 | 16 | 0 | 0 | neutron | laboratory | correlated | conditional-tabular | 24 points, 1.74883e+06–2e+07 |
| Be9 | 102 | 0 | 0 | photon | laboratory | uncorrelated | tabular-angle | 2 points, 1e-05–2e+07 |
| C12 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 719 points, 1e-05–1.5e+08 |
| C12 | 5 | 0 | 0 | neutron | center-of-mass | kalbach-mann | conditional-r,a | 19 points, 1e-05–1.5e+08 |
| C12 | 5 | 1 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 1e-05–1.5e+08 |
| O16 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 1636 points, 1e-05–1.5e+08 |
| O16 | 4 | 0 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 9.43138e+06–3e+07 |
| O16 | 5 | 0 | 0 | neutron | center-of-mass | kalbach-mann | conditional-r,a | 16 points, 1e-05–1.5e+08 |
| Al27 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 60 points, 1e-05–1.5e+08 |
| Al27 | 5 | 0 | 0 | neutron | center-of-mass | kalbach-mann | conditional-r,a | 27 points, 1e-05–1.5e+08 |
| Al27 | 5 | 1 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 1e-05–1.5e+08 |
| Fe56 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 3976 points, 1e-05–1.5e+08 |
| Fe56 | 5 | 0 | 0 | neutron | center-of-mass | correlated | conditional-tabular | 46 points, 1e-05–1.5e+08 |
| Fe56 | 5 | 1 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 1e-05–1.5e+08 |
| U235 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 126 points, 1e-05–3e+07 |
| U235 | 5 | 0 | 0 | neutron | center-of-mass | correlated | conditional-tabular | 56 points, 1e-05–3e+07 |
| U235 | 5 | 1 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 1e-05–3e+07 |
| U235 | 18 | 1 | 0 | neutron | laboratory | uncorrelated | isotropic | none |
| U238 | 2 | 0 | 0 | neutron | center-of-mass | uncorrelated | tabular-angle | 137 points, 1e-05–3e+07 |
| U238 | 5 | 0 | 0 | neutron | center-of-mass | correlated | conditional-tabular | 68 points, 1e-05–3e+07 |
| U238 | 5 | 1 | 0 | photon | center-of-mass | uncorrelated | tabular-angle | 2 points, 1e-05–3e+07 |
| U238 | 18 | 1 | 0 | neutron | laboratory | uncorrelated | isotropic | none |

Li6 has neutron angular laws for MT=2, MT=24, and MT=51–81; every such law is
product=0, distribution=0 and UncorrelatedAngleEnergy. MT=2 has no energy law,
because transport supplies elastic kinematics. MT=24 has a ContinuousTabular
energy law and an independent **laboratory-frame** angular law. MT=51–81 use
LevelInelastic energy laws and center-of-mass angles. All angular table
interpolation codes are 2. MT=2,24,53 have nonuniform angular PDFs; the other
MT=51–81 laws store constant 0.5 PDFs in native tables (they remain tagged
tabular-angle). No reaction meaning was inferred from the MT labels.

At 14.1 MeV, MT=2 selects between 14.0 and 14.2 MeV with upper-table weight 0.5;
both tables have 92 mu points. MT=24 selects between 14 and 15 MeV with weight
0.1; both tables have 21 points. MT=53 brackets 13 and 20 MeV with weight
0.1571428571 (10 and 22 points). All remaining Li6 angular grids have two points
and are isotropic at all incident energies.

## Normalization and native sampler validation

Native histogram integrals use left-bin PDF times bin width; linear integrals
use the exact linear interval integral. No angular distribution in these eleven
nuclides has discrete angular mass. Discrete outgoing-energy lines are retained
separately; they are not angular probability masses. Conditional tabular checks
cover **every** outgoing point at **every** incident point. Kalbach checks select
first/middle/last incident and outgoing indices in each native law.

| Check | Count | Min integral | Max integral | Max absolute error |
|---|---:|---:|---:|---:|
| Independent angular tables | 22,214 | 0.999999876745 | 1.00000014698 | 1.4698e-7 |
| Correlated conditional angular tables | 44,242 | 0.99999976 | 1.00000026 | 2.6e-7 |
| Selected Kalbach conditional parameters | 165 | 1 | 1 | 4.44e-16 |

Tolerance for real-data tabular checks is 1e-6 to allow stored-data precision.
Explicit isotropic distributions integrate analytically to 1.
Selected conditional plots: Be9 MT=16 at E_in=12 MeV, E_out=0.8 MeV integrates
to 0.9999999825; C12 MT=5 at E_in=70 MeV, E_out=168.942 eV integrates to 1.

| Native independent sampler | E_in | Samples | PDF integral | chi-square / df | p-value |
|---|---:|---:|---:|---:|---:|
| Li6 MT=2, product=0, distribution=0 | 14.1 MeV | 200,000 | 0.999999953158 | 46.4043 / 39 | 0.193529 |
| Li6 MT=24, product=0, distribution=0 | 14.1 MeV | 200,000 | 0.9999999345 | 49.0333 / 39 | 0.130265 |

Both use seed=917, 40 fixed mu bins, exact bin integrals of the exposed PDF,
and the actual AngleDistribution::sample method used by transport. Expected
counts are not renormalized. Both native histograms are statistically consistent
with the PDFs. These results verify raw native angular sampling, not the final
laboratory directions of full collision tracking. Synthetic HDF5 tests additionally
exercise histogram and linear interpolation at both grid points and their
midpoint using 100,000 native samples per case and six-sigma bin bounds.

Plots: [Li6 MT=2](angular-validation/Li6-mt2-angular.png),
[Li6 MT=24](angular-validation/Li6-mt24-angular.png),
[Be9 conditional](angular-validation/Be9-mt16-angular.png),
[C12 Kalbach conditional](angular-validation/C12-mt5-angular.png).
The theta plots apply the sine Jacobian and explicitly label density per radian.

## Files changed for this extension

- `src/secondary_data.cpp`: explicit N-body/decay isotropic snapshots and independent native sampling diagnostic, with existing resolver/batch guard.
- `include/openmc/capi.h`: sampling diagnostic declaration.
- `openmc/lib/secondary.py`: independent-angle convenience reader and native diagnostic wrapper; existing snapshot APIs retained.
- `docs/source/pythonapi/capi.rst`: schema, interpolation, normalization, conditional semantics, diagnostic documentation.
- `sandbox/examples/runtime_angular_distribution.py`: real-data inventory, normalization, conditional plotting, native histogram validation.
- `tests/cpp_unit_tests/test_secondary_data.cpp`: histogram/linear HDF5 fixtures, snapshot checks, native histograms and batch rejection.
- `tests/unit_tests/test_lib_angular.py`: real-data normalization, copies, independent/conditional behavior, native sampler comparisons and invalid inputs.
- `tests/unit_tests/test_lib_angular_lifecycle.py`: reads/sampling before initialization, after loading and after finalization; snapshot lifetime.
- This report and four example PNGs.

The workspace already contained the secondary API and energy-mutation work.
Those existing changes are retained. No transport sampler implementation,
secondary energy setter, or reference file was changed by this extension.

## Verification and limitations

- Build: `cmake --build build -j 4` passed (RelWithDebInfo, strict FP enabled, MPI disabled).
- Formatting: clang-format 18.1.8; `git diff --check` passed.
- ENDF/B-VIII.0: 108 Python tests passed across angular, lifecycle, existing secondary, chain, Watt, and emission-spectrum test modules (15.76 s).
- C++: all 13 available CTest targets passed, including distribution, secondary-data and emission-spectrum tests. Full CTest runs 14 targets; `test_mcpl_stat_sum` is reported failed because its test skips when the optional MCPL library cannot be loaded. No MCPL code was changed.
- Existing `runtime_emission_spectrum` regression passed with installed NNDC data (1.29 s). Running it with ENDF/B-VIII.0 first produced a reference mismatch; NNDC is the prescribed reference library. The user's configured cross-section environment and reference files were preserved.
- Real-data inventory, Li6 MT=2/24 histogram comparisons, Be9 conditional and C12 Kalbach conditional plotting scripts completed successfully.
- No affected transport references require regeneration: this extension does not change transport sampling or nuclear data. No MPI build was available for validation.

Representations that cannot be recovered or fully validated here:

1. **Original Legendre/analytic coefficients:** removed during Python HDF5 tabulation. Runtime snapshots expose the actual native table; original coefficients cannot safely be inferred.
2. **Original discrete thermal-angle metadata:** the correlated C++ loader changes interpolation code 0 to histogram and constructs a Tabular. The runtime object no longer identifies original discrete angular probabilities. No such angles were found in the requested real ReactionProducts, and thermal scattering is outside scope. A general discrete-angle probability-sum guarantee is therefore not claimed. Recovering this would require a separate loader/schema change, not manufacturing discrete masses from these snapshots.
3. **General nonstandard Tabular interpolation:** the underlying class has log-linear (4) and log-log (5) inverse sampling branches, exposed by the same Tabular schema. Nuclear-data angular exporters emit only histogram (1) or linear-linear (2); all inspected real angular tables use those supported nuclear-data conventions. Native `Tabular::evaluate` treats non-histogram interpolation as linear even when its sampler uses logarithmic inversion. The example integrates codes 4/5 by their sampling interpolation, but they lack representative real-data and native-histogram validation in this report. Code 3 has no native inverse sampling branch and must not be treated as a safely sampled law.
4. **Arbitrary correlated conditional interpolation and angular marginals:** snapshots retain actual outgoing tables and conditional angular data; they do not construct a marginal or infer a simple analytic p(mu|E_in,E_out) between native points. Joint native histogram validation is not added; the required practical native tests cover independent Li6 MT=2 and MT=24.
5. **Final collision directions:** native angular validation does not validate target thermal motion, laboratory conversion, azimuthal rotation, or material P0 overrides; those are separate transport operations traced above.
6. The semantic RAG tools named by AGENTS.md were not exposed in this session's tool registry. Investigation used source searches across readers, serializers, reaction dispatch, physics consumers, and native samplers.

Read-only work stops here. No angular setter or runtime angular perturbation was implemented.
