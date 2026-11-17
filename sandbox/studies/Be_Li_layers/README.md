# Be_Li_layers: alternating lithium and beryllium cube

This example uses the shared `analysis_codes` suite with its own editable
`model.py` and `config.py`. The geometry and orientation match the requested
arrangement:

- **Cube:** 20 × 20 × 20 cm, centred at the origin; volume 8,000 cm³.
- **Layers:** alternate 1 cm Li and 1 cm Be, starting with Li at x = −10 cm.
  There are 20 slabs, ten of each material, with no gaps.
- **Layer planes:** y-z; their normal/stacking direction is x.
- **Beam:** a 14 MeV monoenergetic pencil beam along **+z**, parallel to the slabs.
- **Default source:** (0.5, 0, −10.01) cm. It enters the centre of Li layer 10,
  avoiding the interface at x = 0. Layer indices begin at zero.
- **Materials:** natural lithium at 0.534 g/cm³ and Be9 at 1.85 g/cm³, both
  at 293.6 K. These are explicit benchmark choices, not a thermodynamic model.
- **Boundary:** transmitting target surfaces inside a vacuum world box whose
  outer faces are at ±11 cm. The 1 cm surrounding region contains no material.
- **Simulation:** fixed source; no eigenvalue calculation.

With this geometry, an unscattered source neutron remains in its initial Li
slab. Neutrons reach other Li/Be slabs after changing direction through
scattering. The beam does not successively cross Li/Be interfaces while travelling
straight along z; that would require stacking the slabs along z instead.

![Geometry in the x-z plane](results/geometry_xz.png)

## Files and settings

- `model.py`: geometry, material construction, automatic source placement and
  `build_model(particles, batches)`. `layer_layout()` exposes slab bounds.
- `config.py`: all physical and Morris/PCE/Sobol settings for this example.
- `results/`: saved baseline demonstration, layout plot, data provenance,
  transport inputs, checkpoints and dry-run cost report.

Key editable settings:

| Setting | Meaning |
|---|---|
| `CUBE_SIDE_CM` | Side length; currently 20 cm. |
| `LI_LAYER_THICKNESS_CM`, `BE_LAYER_THICKNESS_CM` | Width of each alternating slab; both 1 cm. A final slab is clipped if needed to fill the cube exactly. |
| `LAYER_AXIS` | Stacking direction: x, y or z. Currently x. |
| `LI_DENSITY_G_CM3`, `BE_DENSITY_G_CM3` | Separate material densities. |
| `TEMPERATURE_K` | Both materials' temperature. |
| `LI6_ATOM_FRACTION` | None for natural lithium, otherwise a Li6 fraction in [0,1]. |
| `SOURCE_ENERGY_EV`, `SOURCE_DIRECTION` | Source energy and beam direction. |
| `BEAM_ENTRY_LAYER` | Zero-based slab used for an automatic source parallel to the layers. None selects a Li slab near the centre. Set 11, for example, to enter the adjacent Be slab. |
| `SOURCE_POSITION_CM` | Explicit source point override; None enables automatic placement. |
| `SOURCE_STANDOFF_CM`, `WORLD_MARGIN_CM` | Entry-point standoff and outer vacuum-world margin. |

Automatic placement supports beams along any signed coordinate axis. Oblique
beams require an explicit source position. Make sure that position lies inside
the simulated world and points toward the intended target region.

## Tritium-production QoI

The configured total is **Li/Be-produced tritium**, using the verified
ENDF/B-VIII.0 partial reaction scores:

```python
PRODUCTION_CHANNELS = {
    'Li6': [105],
    'Li7': list(range(52, 83)),
    'Be9': [700, 701],
}
```

Rates are in tritons per source neutron. The total retains within-batch
correlations. Cached MT205 and redundant Be9 MT105 are not the primary scores.
To study lithium-produced tritium only, remove Be9 from `PRODUCTION_CHANNELS`
and change `RESPONSE_LABEL`; Be remains in transport and perturbation discovery
with `TARGET_NUCLIDES=None`. No exclusion list is needed.

## Baseline demonstration already run

The saved demonstration used 1,000 particles per batch, ten batches, seed
20261103: **10,000 source histories**. See
[baseline_summary.json](results/baseline_summary.json) for all constituent means
and standard errors. It is a small transport demonstration, not a completed or
converged Morris/PCE/Sobol campaign.

An extra tally records **Li6 MT105 reactions per source neutron by slab**:

![Li6 production by slab](results/baseline_Li6_mt105_by_layer.png)

Its layer sum was checked against the total Li6 response. This figure is not
all-constituent tritium: Li7 and Be9 contributions are reported separately in
the baseline summary. The CSV includes one standard error per layer. Standard
errors should not be added as if the bins were independent.

The baseline XML files, including that extra spatial tally, are in
`results/baseline_inputs/`. They can be run directly with OpenMC:

```bash
cd /sandbox/openmc
OMP_NUM_THREADS=1 /opt/openmc-venv/bin/python - <<'PY'
import openmc
openmc.run(cwd='Be_Li_layers/results/baseline_inputs')
PY
```

The exported settings suppress statepoints; enable a desired statepoint in the
XML/model settings if using this direct command for new post-processing. The
shared sensitivity workflow builds its own response tallies and does not
currently fit spatial Sobol models for every layer.

## Run the shared analysis

From `/sandbox/openmc`, with this fork's Python/shared-library environment:

```bash
export OPENMC_CROSS_SECTIONS=/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# Native discovery and cost report, without transport:
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers --dry-run

# Full configured Morris -> per-family selection -> high-fidelity PCE/Sobol:
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers --run

# Resume successful deterministic checkpoints:
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers --resume

# Refit/replot saved high-fidelity data without new transport:
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers --reanalyse
```

The full configured sensitivity campaign has **completed**. See
[the main report](results/summary.md), [all Morris/Sobol plots](results/sensitivity_plots.md),
and [final execution checks](results/full_campaign_checks.json). `SCREEN_PARTICLES * SCREEN_BATCHES`
controls Morris histories; `HIGH_PARTICLES * HIGH_BATCHES` controls PCE transport
histories. `TOP_PER_FAMILY` selects separately in XS, energy and angle.
`PCE_BACKEND` accepts `numpy` or optional `openturns`. `WORKERS` and `HIGH_WORKERS`
control the two stages' process counts. The priors are hypothetical ±20%
coefficient ranges, not evaluated covariance uncertainties.

For changed geometry or transport settings, use a fresh output directory:

```bash
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers --dry-run --output sandbox/studies/Be_Li_layers/results_variant
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers --run --output sandbox/studies/Be_Li_layers/results_variant
```

The [shared description](../analysis_codes/description.md) explains the complete
configuration, supported native mutations, screening limitations, validation
criteria and output files.


## Completed full campaign

Executed with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers --run
```

- 171 eligible parameters; 20 Morris trajectories; 3,444 low-fidelity
  observations including baseline and timing checks.
- Morris convergence criterion met: split-half ranking correlation 0.92825,
  top-10 overlap 90%, median relative interval width 10.01%.
- Selected 30 parameters: 12 XS, eight energy modes, ten angle modes.
- 512 high-fidelity training observations, 96 independent validation observations
  and one high-fidelity baseline, each with 100,000 source histories.
- 96,740,000 campaign histories, plus the earlier 10,000-history baseline;
  within the configured 160-million cap.
- All 45 software tests passed. Native controlled checks, exact restoration,
  prior bounds, PDF normalizations and saved-PCE prediction replay also passed.
- All four responses accepted the configured independent validation checks:
  total R² = 0.996609, Li6 = 0.976951, Li7 = 0.999148, Be9 = 0.994001.
- 247 PNG and 247 PDF figures are saved across the galleries and result folders.

The conditional total-production distribution has mean **0.217198**,
standard deviation **0.009244**, and a 5th–95th percentile interval of
**[0.202324, 0.232035] tritons per source neutron**. The high-fidelity baseline is
0.216739 ± 0.000297 (one standard error).

The leading fitted total-order contributions are Li7 MT56 XS (79.02%),
Be9 MT2 XS (9.13%), Li6 MT105 XS (2.68%) and Li7 MT2 XS (2.20%). These fractions
are conditional on the retained inputs and hypothetical independent coefficient
priors. Weak modes remain flagged unresolved. The selected order-2 hyperbolic
PCE contains only univariate terms, so its first- and total-order indices are
equal; this does not prove there are no physical interactions.

![Independent PCE validation](results/high/pce_validation.png)
![Every retained Sobol parameter](results/high/sobol_indices.png)
![Conditional production distribution](results/high/tritium_distribution.png)

[Complete Morris values](results/low/morris_results.csv) ·
[All Sobol values](results/high/sobol_indices.csv) ·
[Native selected perturbation plots](results/nuclear_data_plots/index.md) ·
[Full log](results/full_campaign.log).
