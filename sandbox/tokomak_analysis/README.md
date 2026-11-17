# DEMO-inspired Li₂O tokamak — native OpenMC CSG

This directory contains a configurable, full-360° neutronics model, its executed examples and validation evidence. It uses **native CSG**, with OpenMC transporting particles through the same cells that are queried for visualisation. There is no DAGMC geometry or second transport engine.

This is an engineering-inspired research starting point, **not a qualified DEMO engineering design**. Published parameters and explicit assumptions are separated in [references.md](references.md). Li₂O deliberately replaces the reference HCPB silicate/titanate breeder. See [implementation_report.md](implementation_report.md) for measured results, acceptance status and limitations, and [gallery.md](gallery.md) for images.

## Run from the repository root

Use the Python environment containing the compiled modified OpenMC fork:

```bash
cd /sandbox/openmc
export OPENMC_CROSS_SECTIONS=/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=2
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.main \
  --config sandbox/tokomak_analysis/examples/level2.json \
  --output sandbox/tokomak_analysis/results/my_segmented_run \
  --particles 2000 --batches 10 --validate --run --debug --plot --vtk
```

Without `--run`, the command builds and exports XML without transport. `--validate` performs recursive Python CSG overlap/coverage point checks. `--debug` enables native overlap checking during an actual `--run`; it is not an exhaustive geometric proof. `--plot` samples native cell queries for eight diagnostic slices. `--vtk` samples the native geometry for ParaView. Separate processes initialise/finalise each example. Do not use an output directory simultaneously from two processes.

Reproduce the executed transport cases:

```bash
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.main --fidelity 1 --particles 1000 --batches 10 --run --debug --output sandbox/tokomak_analysis/results/global
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.main --fidelity 2 --particles 2000 --batches 10 --run --debug --output sandbox/tokomak_analysis/results/segmented
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.main --fidelity 3 --particles 500 --batches 10 --run --debug --output sandbox/tokomak_analysis/results/detailed
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.main --config sandbox/tokomak_analysis/examples/coupled_smoke.json --run --debug --output sandbox/tokomak_analysis/results/coupled_smoke
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.main --fidelity 3 --plot --vtk --output sandbox/tokomak_analysis/results/detailed
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.main --config sandbox/tokomak_analysis/examples/first_wall_zoom.json --vtk --output sandbox/tokomak_analysis/results/first_wall_zoom
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.neutronics.plots sandbox/tokomak_analysis/results/segmented
/opt/openmc-venv/bin/python -m sandbox.tokomak_analysis.inventory
OMP_NUM_THREADS=2 /opt/openmc-venv/bin/python -m pytest -q tokomak_analysis/tests
```

Reusing an output directory updates its XML/results. Use a new directory to preserve an old calculation. Geometry sampling/export can load nuclear data but does not run histories. More histories are needed for individual module and especially deep-shielding observables; increase `particles` and/or `batches` in your JSON. No large simulation campaign or UQ scan is launched automatically.

## Files and responsibilities

| File/module | Purpose |
|---|---|
| `config.py` | Frozen typed dataclasses, dimensions, fractions, temperatures, source, execution and export settings; validation and JSON I/O |
| `materials.py` | Isotopic compositions, solid/bed densities and volume-consistent homogenisation |
| `plasma.py` | Faceted Miller D contour and consistent normal-offset contours |
| `geometry/primitives.py` | Revolved RZ polygons, exact signed conical RZ boundaries, sectors, boxes and component metadata |
| `geometry/assembly.py` | Sector universe hierarchy, component assembly and OpenMC model |
| `geometry/blanket.py`, `breeder_units.py` | Toroidal/poloidal modules, first wall, breeder/multiplier/cooling cycles, sidewalls, ribs and manifolds |
| `geometry/vessel.py` | Double-wall vessel, water/steel filler and thermal shield |
| `geometry/central_column.py`, `magnets.py` | Support, shielding, segmented CS, 16 discrete D-shaped TF coils, six PF coils |
| `geometry/divertor.py` | Two inclined W/CuCrZr/water/steel targets and cassette |
| `neutronics/source.py` | Plasma-volume DT sources and fusion-power normalisation |
| `neutronics/library.py` | Library hashing and pinned production-channel verification |
| `neutronics/tallies.py` | Production, module, flux, spectrum, current, heating and magnet tallies |
| `neutronics/postprocess.py` | Batch-correlated sums, units, structured results and mesh-flux VTK |
| `neutronics/plots.py` | Replot saved breeding contributions and blanket spectrum |
| `visualization/export_vtk.py` | Native cell-query slices and VTI/VTU/VTP/VTM export |
| `visualization/paraview.py` | Optional ParaView camera/colour/cutaway screenshot helper |
| `visualization/render.py` | CPU/Agg rendering of the exported VTK surfaces |
| `validation.py`, `inventory.py` | Coverage/overlap sampling and paired fidelity-inventory comparisons |
| `main.py` | Command-line orchestration; XML export and library transport |
| `data/` | Small copied evidence records for the verified evaluation mappings; no runtime dependency on earlier UQ directories |
| `tests/` | Geometry, materials, source, native transport and VTK tests |

## Geometry and fidelity

Coordinates are Cartesian X/Y/Z, the toroidal axis is Z, and **all input lengths are centimetres**. Plasma reference values are major radius 910 cm, minor radius 290 cm, elongation 1.65 and triangularity 0.33. The 32-point contour is revolved using native planes/cylinders/cones, not an imaginary elliptical-torus primitive. Scrape-off space connects to the lower divertor vacuum; no magnetic equilibrium or field-line connection is calculated.

The blanket build is 75.5 cm inboard and 98.2 cm outboard. A 0.2 cm tungsten armour is followed by separate 0.8 cm structural, 1 cm cooling and 0.2 cm rear-wall layers. Three active breeder/multiplier/cooling cycles precede a 3 cm backplate and 12 cm manifold approximation. There are 16 repeated sectors, three toroidal modules per sector and 16 vertical poloidal stacks (eight inboard and eight outboard): **768 nominal module instances**. The lower blanket truncation leaves a divertor opening. Sidewalls and inter-module gaps remain explicit.

- **Level 1:** Active blanket interiors homogenised with volume-derived composition; first wall, gaps, segmentation and back support retained.
- **Level 2:** Explicit radial breeder, multiplier and cooling layers in each module; cooling/structure volume fractions homogenised within their layers.
- **Level 3:** Additional grouped helium slots, structural ribs and breeder-unit subdivisions. These are CSG engineering approximations, not individual pebble contacts, curved piping or manufacturing drawings.

The six-module CS has casing/insulation, the 16 D-shaped TF coils have winding/insulation/case regions, and PF coils have separate case/winding regions. Vessel walls, shield/filler, central supports and a recognisable lower divertor are separate materials/cells. Detailed dimensions of these systems are assumptions where references do not establish them. Shared surfaces are deduplicated **within universes**, not globally across rotated coordinates.

## Change the model/settings

Edit a copy of `examples/level2.json`, then pass `--config your_file.json`. The dataclass in `config.py` lists every setting and default; unknown JSON keys fail rather than being ignored. CLI `--fidelity`, `--particles` and `--batches` override the JSON. Example programmatic use:

```python
from dataclasses import replace
from pathlib import Path
from sandbox.tokomak_analysis.config import Config
from sandbox.tokomak_analysis.geometry.assembly import build_model
from sandbox.tokomak_analysis.main import export_inputs, run_transport

c = replace(Config(), enrichment=(0.60, 0.70, 0.80),
            breeder_packing=0.64, fidelity=2, particles=10000, batches=20).validate()
a = build_model(c, output=Path('sandbox/tokomak_analysis/results/custom'))
export_inputs(a, Path('sandbox/tokomak_analysis/results/custom'))
# a.model is an ordinary openmc.Model: inspect or extend its tallies here.
run_transport(a, Path('sandbox/tokomak_analysis/results/custom'), debug=True)
```

Change `major_radius`, `minor_radius`, shaping, blanket builds, coil positions (`pf_coils` is a list of named radius/z/width/height records), module counts and layer fractions through config. Divertor target/cassette offsets are relative to the plasma major radius and lower plasma boundary. Large changes require fresh geometry checking; validation of scalar ranges cannot guarantee a globally nonintersecting reactor.

`enrichment` has one Li-6 atom fraction per breeder cycle. `breeder_packing` is solid volume fraction: the rest is helium purge gas, not a reduction of solid density without accounting for void. Structural fractions are additional explicit mixture controls. `multiplier` accepts `Be` or `Be12Ti`; the latter is an approximate alternative, not separately transport-qualified here. Steel compositions and winding-pack fractions are documented approximations. Material temperatures and nearest-table tolerance are configurable; the 4.5 K magnet operating label does **not** imply available cryogenic evaluated neutron data.

## Source and power

Default emission is isotropic 14.1 MeV, uniform in physical plasma volume. Native cylindrical proposals use density proportional to R and reject outside the plasma cell. `radial_profile_power > 0` selects the installed native `TokamakSource` profile. Nonzero `poloidal_asymmetry` uses a deterministic, plasma-rejected source bank of configurable size; this is only source generation, not secondary-particle transport sampling. `source_energy_sigma_ev` optionally gives a narrow Gaussian source spectrum.

At `fusion_power_mw=2000`, the assumed 17.6 MeV per DT fusion gives **7.092624×10²⁰ source neutrons/s** and 1602.27 MW source-neutron power. Tallies remain per source neutron until explicitly multiplied by this rate. No plasma burn physics, time dependence or depletion is included.

## Production and other responses

The primary production scores use current native partial cross sections, explicit nuclides and tracklength scoring:

- Li-6: MT105, one triton/reaction.
- Li-7: MT52–82, one triton/reaction for the verified evaluation.
- Be-9: MT700/701, verified triton channels.

The exact three neutron files are hash-pinned to the verified ENDF/B-VIII.0 evidence in `data/`; baseline partial sums are checked against loaded MT205 before transport. Cached MT205 is **not** the primary tally. Unknown libraries fail verification, requiring new channel/yield evidence rather than silently reusing mappings.

`lithium_TBR` sums Li channels; `TBR` sums verified Li+Be channels. This is a **scoped verified TBR**, not a claim that oxygen, steel and every other constituent produce zero tritium. Channel/module contributions are summed within each batch, preserving covariance in their standard errors. `module_breeding.csv` and the inboard/outboard report retain spatial contributions; their batch sum is asserted equal to the global production score.

Other tallies include volume-integrated flux, energy spectra, Cartesian neutron flux, world leakage, signed central-shield surface current, (n,2n) event-rate indicator, VV/magnet spectra and >0.1 MeV exposure. The (n,2n) tally alone is not a complete neutron balance. Heating is `heating-local` in neutron-only mode; `coupled_photons=true` uses collision-estimator `heating` for coupled neutron/photon deposition. Do not equate those approximations without checking the physics.

Main outputs:

- `results.json`: per-source production, batch errors, per-second scaling and timing.
- `batch_tritium.csv`, `module_batch_production.npz`, `module_breeding.csv`: correlated production data.
- `component_heating.csv`: integrated component energy deposition/source and MW at configured power.
- `integrated_tally_totals.csv`, `tally_manifest.json`, named `.npz`: raw tally integrals/errors, bin metadata and units.
- `neutron_flux.vti`: bin-volume-normalised mesh flux per source and per second.
- `statepoint.*.h5`: standard OpenMC results, independently accessible via `openmc.StatePoint`.
- `components.json`, `cell_paths.json`, `build_metrics.json`: IDs, hierarchy and counts.
- `nuclear_data_provenance.json`: code/library/data identities and configuration context.

Raw volume-integrated flux has units cm/source; divide by region volume for flux density cm⁻²/source. Spectral arrays are energy-bin integrals, not PDF values or automatically per-eV densities. Tally errors are one batch standard error; fine mesh/magnet zeros at these statistics are unresolved.

## ParaView and plots

Open `results/detailed/visualization/reactor.vtm` in ParaView. Keep its companion `reactor/` directory beside it. The VTM contains separate component VTU blocks; enable/disable them in the multiblock inspector. `reactor.vti` contains the full image grid; individual `.vtp` files contain sampled component surfaces. Fields include material/cell/component/module/sector IDs and inboard/outboard classification; names/colours are in JSON metadata.

```bash
# In an environment with ParaView's pvpython and this repository on PYTHONPATH:
pvpython -m sandbox.tokomak_analysis.visualization.paraview \
  sandbox/tokomak_analysis/results/detailed/visualization/reactor.vtm \
  --cutaway --screenshot tokamak_paraview.png
```

The export queries the authoritative CSG through `openmc.lib.find_cell` at voxel centres. **It is not an exact boundary mesh.** Full-reactor voxels are about 39.6×39.6×31.3 cm, so thin armour may disappear. The separate `first_wall_zoom` example uses 0.1×0.75×0.5 cm voxels. Increase resolution or restrict `vtk_origin`/`vtk_extent` for other thin structures. Metadata records resolution and limitations. VTK readback is tested; the ParaView GUI was not available here. The environment's OpenGL renderer crashed, so supplied full/cutaway images were rendered with Matplotlib Agg from the same VTK surfaces:

```python
from pathlib import Path
from sandbox.tokomak_analysis.visualization.render import render
p = Path('sandbox/tokomak_analysis/results/detailed/visualization')
render(p/'reactor.vtm', p/'cutaway.png', cutaway=True)
```

Diagnostic XY/XZ/YZ, two toroidal and module/divertor slices also use native CSG cell queries. They were visually inspected. No independently drawn reactor is substituted for geometry validation.

## Scope and next work

Automated and debug-transport checks passed for the supplied baseline; they do not certify every possible config. Re-run them after geometry changes. Production demonstrations are intentionally modest (10k/20k/5k histories); they are not converged shielding or optimisation calculations. Upgrade to validated engineering dimensions, penetrations, realistic coolant routing, thermal/mechanical design and higher-statistics coupled transport before design conclusions. The ordinary `openmc.Model`, cell metadata and reaction-specific tallies can be connected to the existing sensitivity framework, but no tokamak Morris/PCE campaign is claimed here.
