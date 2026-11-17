# Be_Li_layers_x_beam: beam crossing the layers

This study repeats the `Be_Li_layers` sensitivity workflow with the beam along
**+x**, normal to the alternating material slabs. The original +z study is
preserved in `../Be_Li_layers/`.

- Cube: 20 × 20 × 20 cm, centred at the origin.
- Ten 1 cm Li slabs alternating with ten 1 cm Be slabs, stacked along x.
- Beam: monoenergetic 14 MeV, directed +x, starting at (−10.01, 0, 0) cm.
- The initial straight beam crosses Li, Be, Li, Be, … in sequence.
- Lithium: natural isotopic composition, 0.534 g/cm³.
- Beryllium: Be9, 1.85 g/cm³. Both materials are at 293.6 K.
- Vacuum world: faces at ±11 cm; fixed-source transport.

![Beam through the layers](results/geometry_xz.png)

`model.py` defines the geometry/source; `config.py` contains all study settings.
Only the beam direction, study identifier and high-stage process count differ
from the original configuration. Automatic source positioning puts the beam
outside the negative-x face. Changing direction to +z automatically returns it
to a central lithium slab; use an explicit position if desired.

The QoI includes verified production in Li6 MT105, Li7 MT52–82 and Be9 MT700/701,
in tritons per source neutron, summed within batches. Nuclear-data parameters
are discovered for all three constituents and source-inaccessible reactions
are removed from the perturbation registry.

## Configured analysis

- 20 Morris trajectories, 10,000 histories per low-fidelity observation.
- Top ten parameters independently per family, plus reserved physics checks.
- 512 new high-fidelity training runs and 96 held-out validation runs.
- 100,000 histories per high-fidelity observation; eight spawned workers.
- NumPy weighted/sparse Legendre PCE, training-only order selection over 2 and 3.
- Hypothetical independent ±20% coefficients, with the existing two-mode PDF
  shape budget; these are not evaluated covariance uncertainties.
- Explicit 160-million-history cap; planned campaign cost 96.74 million histories.

From `/sandbox/openmc`:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers_x_beam --run
# Same command resumes matching checkpoints, or use --resume explicitly.
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m sandbox.analysis_codes Be_Li_layers_x_beam --reanalyse
```

Results are under `results/`, including `full_campaign.log`, `summary.md`,
`sensitivity_plots.md`, native data galleries and the high-fidelity datasets.
See [the shared description](../analysis_codes/description.md) for the complete
methodology and configuration reference.


## Completed results and comparison

The configured campaign completed: 3,444 low-fidelity observations and 609
high-fidelity observations (baseline + 512 training + 96 validation). Native
production/restoration checks and all 45 software tests passed. An additional
100,000-history independent serial replay matched all parallel sample means,
standard errors and batch scores bit for bit. Total transport including that
check was **96,840,000 histories**, within the 160-million cap.

**Screening limitation:** at 20 trajectories the top-10 overlap was 100%, but
split-half ranking correlation was 0.89260, below the configured 0.90 threshold.
Median relative interval width was 15.05%. Full Morris convergence was not
established; the reduced analysis is conditional on the selected 28 parameters
(10 XS, eight energy, ten angle), not a converged full-space attribution.

All four PCE responses passed configured independent validation: total R² =
**0.979247**, Li6 = 0.976496, Li7 = 0.998100, Be9 = 0.998562. The selected
order-2 hyperbolic model is additive, so equal first/total indices do not prove
absent physical interactions. Held-out validation does not remove the screening
limitation above.

The conditional total-production distribution has mean **0.161747**, standard
deviation **0.007822**, and 5th–95th percentile interval **[0.148919, 0.174868]**
tritons per source neutron. These values apply to the retained hypothetical
input distributions, not evaluated nuclear-data covariance uncertainties.

The leading conditional total-order fractions are Li7 MT56 XS (44.81%), Be9
MT2 XS (38.91%), and Li6 MT105 XS (11.18%). The selected sets differ between
beam orientations; their normalized Sobol fractions should not be treated as
contributions under an identical full-space prior.

### Unperturbed beam-orientation comparison

Both high-fidelity baselines used 100,000 source histories. The change shown is
a point-estimate ratio; the figure shows one standard error for each estimate.

| Production | Parallel +z | Crossing +x | Relative change |
|---|---:|---:|---:|
| Total | 0.216739 | 0.161106 | −25.67% |
| Li6 | 0.015604 | 0.027914 | +78.89% |
| Li7 | 0.197670 | 0.118752 | −39.92% |
| Be9 | 0.003464 | 0.014440 | +316.85% |

Units: tritons per source neutron. The beam-orientation change also changes the
entry point to the appropriate upstream face, as defined by the model builder.

![Baseline comparison](results/baseline_beam_comparison.png)
![Held-out PCE validation](results/high/pce_validation.png)
![All selected Sobol inputs](results/high/sobol_indices.png)
![Conditional production distribution](results/high/tritium_distribution.png)

[Main report](results/summary.md) ·
[All Morris/Sobol figures](results/sensitivity_plots.md) ·
[Comparison CSV](results/baseline_beam_comparison.csv) ·
[Final checks](results/full_campaign_checks.json) ·
[Parallel/serial replay](results/parallel_serial_replay.json).

### History comparison with the LiPb study

This run matches the previous layered Li/Be study's history settings. The saved
LiPb study used **1,000,000** histories per high-fidelity observation; both layer
studies used **100,000**. Both used 10,000 per Morris observation, 512 training
and 96 validation points. LiPb had 14,144 low-stage evaluations versus 3,444 for
the layers. The +x study used eight high-fidelity workers, whereas the saved
LiPb and original +z studies used one. Concurrency reduces wall time; it does
not increase histories per observation.
