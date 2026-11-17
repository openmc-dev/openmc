# DEMO-inspired Li₂O CSG implementation and verification report

Date: 6 November 2026. Repository: `/sandbox/openmc`, existing `xs-change` branch.

## Outcome

A full-360° parametric native-CSG tokamak, three blanket fidelities, volumetric DT source, verified reaction-specific production scores, heating/transport/shielding scores, automated tests and native-query VTK export have been implemented and exercised. This is a usable **engineering-inspired neutronics research model**, not a validated reconstruction of a frozen DEMO/HCPB reference CAD design or an optimised breeder blanket. Li₂O is a deliberate departure from the referenced HCPB ceramic blend.

See [README](README.md) for commands and configuration, [references](references.md) for published values versus assumptions, and [gallery](gallery.md) for images. Complete executable XML inputs and standard OpenMC statepoints are in each successful example directory.

## Phase A — Audit and research

The requested `tokamak_csg.py` prototype could not be located in this checkout or the searched `/sandbox` and `/tmp` trees; no prototype geometry was presumed correct. Existing uncommitted native-runtime and sensitivity work was preserved. No native C++, C API, nuclear-data file or older analysis directory was changed for this model.

Installed OpenMC is `0.1.dev16025+gfd99f66a2.d20260924`, using `/sandbox/openmc/openmc/lib/libopenmc.so`. The configured library is `/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml`. Individual evaluated files, compiled-library hash, Git revision/dirty state, source-code hash and model context are recorded per calculation in `nuclear_data_provenance.json`. Historical saved hashes represent the code actually used for those runs; subsequent CLI-only validation repair does not retroactively change provenance.

Local inspection verified the actual `openmc.model.Polygon(..., basis='rz')`, rotated-universe, native `TokamakSource`, rejection constraints, CellInstanceFilter, native tally SUM and source sampling APIs. No arbitrary elliptical torus was assumed. Primary DEMO design-space, magnet, HCPB, shielding, containment, divertor and materials references are cited in `references.md`; nonpublished detailed builds are explicitly assumptions.

## Phase B — Global geometry

Implemented plasma and chamber/SOL/divertor-connected vacuum, central supports and shield, six individually cased/insulated CS modules, 16 separate D-shaped TF winding/case/insulation regions, six cased PF coils, double-wall vessel and filler, thermal shield, and two inclined layered lower divertor targets with cooled cassette. The source plasma is a faceted Miller contour; it is not a Grad–Shafranov solution or an X-point magnetic equilibrium.

Sector-universe reuse represents all 360 degrees. Exact signed conical RZ halfspaces are used for inclined target and clipping boundaries. Early independent near-coincident polygon decomposition caused native navigation/debug errors; sharing adjacent interfaces and using exact signed boundaries resolved those errors. Surface merging is scoped to a universe so that rotated local-coordinate surfaces are not incorrectly conflated. Final executed configurations below completed native geometry-debug transport.

## Phase C — Blanket engineering and inventories

All levels retain first-wall armour/steel/coolant/rear wall, module gaps and sidewalls, backplate and manifolds. Three active cycles alternate Li₂O breeder, Be multiplier and cooling layers. Level 1 homogenises the active interior, level 2 resolves these radial layers, and level 3 adds grouped helium slots, ribs and breeder-unit subdivisions. There are 768 nominal blanket module instances, with independent scores and IDs. Poloidal segmentation is by separate inboard/outboard vertical stacks, truncated below to form the divertor opening. It is not a published commercial-module count.

Li₂O solid density is 2.013 g/cm³, packing fraction 0.64, with purge helium accounted for in the remaining volume. Default Li-6 enrichment is 60 atom% in each breeder cycle. Be is similarly packed at 0.64. Steel, winding packs, coolant plates and manifolds use explicit documented volume-fraction mixtures. Theoretical solid density and bed density are not interchanged.

A common-point **100,000-point native-query inventory calculation** compares all three levels. Versus level 2, estimated Li inventory changes are +0.188% (level 1) and −0.115% (level 3), with paired one-standard-error uncertainties 1.954% and 0.487%. Be changes are +0.414% ±2.133% and −0.179% ±0.516%. Thus the estimates are consistent with similar inventories, but do **not** prove exact equality or a sub-percent tolerance. Global homogenisation and lower truncation can produce real differences. Tables and methodology: [inventory](results/inventory/fidelity_inventory_comparison.csv).

## Phase D — Native transport and responses

All tabulated uncertainties are **one batch standard error**. Timings include initialisation, verification and result export, with two OpenMP threads and geometry debugging; they are not pure collision-kernel benchmarks. Peak RSS is process high-water memory and may include library/geometry initialisation. Construction of the three main models took approximately 1.0–1.1 seconds.

| Case | Histories | Unique cells / physical leaf instances | Surfaces | Materials | Lithium tritons/source | Verified Li+Be tritons/source | Elapsed s | Histories/s | Peak MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| global | 10,000 | 220 / 7,072 | 1537 | 13 | 1.225941 ± 0.017582 | 1.232607 ± 0.017546 | 13.27 | 753.3 | 945 |
| segmented | 20,000 | 348 / 13,216 | 2049 | 17 | 1.166260 ± 0.008255 | 1.177362 ± 0.008234 | 29.14 | 686.3 | 953 |
| detailed | 5,000 | 524 / 21,664 | 2529 | 14 | 1.135656 ± 0.018225 | 1.146997 ± 0.018181 | 11.34 | 441.1 | 874 |
| coupled_smoke | 400 | 220 / 7,072 | 1537 | 13 | 1.203883 ± 0.042690 | 1.210869 ± 0.042399 | 4.94 | 80.9 | 956 |

The global/segmented/detailed runs use 10 batches, and coupled smoke uses 4. All use seed 20261106. Segmented is the higher-statistics demonstration (20,000 histories), not a production-quality shielding campaign. Level-1 homogenisation changes transport and the results should not be treated as statistically identical across fidelities. Level-2/3 agreement is suggestive at these statistics, not a convergence study.

The 2 GW assumed fusion power corresponds to 7.092624×10²⁰ neutrons/s at 17.6 MeV/fusion, and 1602.27 MW source-neutron power. Raw scores and per-second conversions are both retained. TBR denotes verified tritons per DT source neutron, independent of source-power scaling.

Production uses current Li-6 MT105, Li-7 MT52–82 and Be-9 MT700/701 partial reaction rates, with verified unit triton yields for these pinned evaluations. The summed baseline production cross sections are compared against loaded MT205 before running. MT205 is a verification diagnostic, **not** the primary score. The Li-only response is separately retained. Production in oxygen and other structural constituents has not been exhaustively mapped; reported `TBR` is explicitly scoped to verified Li+Be production. It is not a complete all-isotope yield inventory.

All channel/module totals are summed **within each batch** before estimating uncertainty. The per-module batch sum is checked against the global production batch sum, preserving covariance. Inboard/outboard and poloidal-stack breeding reports, energy-binned flux, volume-integrated and mesh flux, boundary leakage, shield-interface signed current, (n,2n) event rate, component heating, and VV/TF/PF/CS exposure are saved. Neutron-only local heating and coupled neutron/photon heating use distinct documented score/estimator combinations. Coupled transport was demonstrated for 400 histories only; it is not a precise photon-heating assessment. Magnet scores that are zero at these histories are unresolved, not evidence of adequate shielding.

## Phase E — Visualisation

`results/detailed/visualization/reactor.vtm` is a component multiblock with VTU children, alongside a full VTI and extracted component VTP surfaces. Metadata retains cell, material, component, module, sector and inboard/outboard identity. **All classification comes from native `openmc.lib.find_cell` on the CSG.** No second independently defined reactor mesh is used for transport.

Both full-reactor and focused first-wall files passed VTK readback. The coarse full export has 737,280 voxels of approximately 39.58×39.58×31.25 cm. It cannot faithfully resolve millimetre armour. A focused export has 576,000 voxels at 0.1×0.75×0.5 cm, with 9,676 sampled tungsten voxels. These are sampled visualisations, not exact analytic surfaces. Sub-voxel structures can disappear; small holes in coarse rendered magnets are voxel artefacts, not intentional geometry penetrations.

XY, XZ, YZ, two additional toroidal cuts and inboard/outboard/divertor close-ups were generated and visually inspected. ParaView-readable VTK files were verified with VTK readers; **the ParaView GUI was not available and is not claimed to have been opened**. Headless OpenGL rendering failed in this container. Supplied cutaway/full views instead use Matplotlib Agg rendering of the same VTK-extracted surfaces. An optional ParaView script and colour table are provided. An exploded-view implementation was not added.

![Native CSG poloidal section](results/detailed/visualization/xz_full.png)
![Explicit breeder/multiplier/cooling subdivision](results/detailed/visualization/blanket_outboard_detail.png)
![Sampled full-reactor cutaway](results/detailed/visualization/cutaway.png)
![Poloidal breeding, batch standard errors](results/segmented/poloidal_breeding.png)

## Phase F — Tests and acceptance evidence

Automated tests cover solid/bed composition and density, enrichment, optional multiplier, signed RZ halfspaces, plasma volume, random recursive overlap/coverage, XML/native-surface types, sector periodicity, invalid config, module count and fidelity envelopes, native source sampling, asymmetric source-bank sampling, native debug transport and module conservation, VTK roundtrip, and CLI validation dispatch. The final result is **20 passed in 17.71 s**; see [results/tests.log](results/tests.log).

Geometry debug completed the three listed fidelity runs and the coupled smoke run without a reported lost-particle/overlap failure. Native cell-query checks also completed without failures. Finite sampling does **not** prove the absence of all overlaps: the 20k segmented debug report identifies 36 cells with fewer than ten overlap checks. Additional targeted tracks/points and higher-statistics qualification remain appropriate. Some thin or shielded cells receive little transport coverage.

The installed Zr96 library emits a negative probability-table warning at 600 K. It is retained in logs and the nuclear data were not modified to suppress it. Nearest-temperature selection is a documented model approximation, especially for cryogenic winding regions. Native plotter/OpenGL backend failures during development were avoided with native-query plots and CPU rendering; they are not represented as successful native image-rendering tests.

| Requested capability | Status/evidence |
|---|---|
| Full 360° native CSG, distinct major subsystems | Implemented; XML, native point queries and debug transport |
| Meaningful Li₂O HCPB engineering segmentation | Implemented at three approximation levels; no manufacturing qualification |
| Configurable dimensions/material/source/run/export settings | Implemented in typed config and JSON examples |
| Fixed-source transport, production and shielding tallies | Executed; deep-shielding statistics unresolved |
| 3D ParaView-compatible dataset | VTK readers verified; ParaView GUI check outstanding |
| Automated tests | Executed; see current test log |
| Published versus assumed dimensions | Explicitly separated in references |
| Exhaustive geometry/engineering validation | Not established by finite transport/point tests |
| Blanket thermal/mechanical or magnetic self-consistency | Not implemented or claimed |

## Remaining physics and engineering limitations

No ports, penetrations, remote-handling details, non-axisymmetric divertor cassette gaps, realistic coiled coolant routing, mechanical loads or magnetic equilibrium calculation are included. Cooling slots are grouped approximations; pebbles are homogenised. PF/CS placements and winding packs are engineering assumptions, not electromagnetic design results. EUROFER/SS316L/insulation and winding mixtures are representative, not certified batch specifications. The hot-water model omits a matched thermal scattering treatment. Magnet transport temperatures use available evaluated data, not a true 4.5 K nuclear-data library.

Li₂O changes lithium inventory, oxygen content and thermophysical behaviour relative to reference HCPB ceramics; neutronics calculations do not qualify its irradiation stability, swelling, purge chemistry or heat extraction. A TBR above unity here does not establish plant self-sufficiency after penetrations, outages, processing losses and inventory requirements. More histories, library sensitivity checks and an all-isotope production audit are needed before physical design conclusions. No UQ/sensitivity campaign was run as part of this model-development task.

Recommended next steps are to replace assumed subsystem builds with a selected engineering baseline; add realistic ports/cooling routes; target poorly sampled thin/interface regions with additional geometry tests; then run converged coupled transport and validated shielding/heating studies before optimisation or UQ. Reproduction commands are in the README, and every run saves its effective config and input XML.
