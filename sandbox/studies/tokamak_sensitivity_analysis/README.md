# Tokamak nuclear-data sensitivity study

This directory plugs the existing shared Morris/PCE/Sobol runner into the native CSG model in `tokomak_analysis/`. Geometry and materials remain configured in `tokamak_sensitivity_analysis/config.py` through `TOKAMAK_CONFIG`; run budgets and PCE settings are in the same file. Results and checkpoints go to `tokamak_sensitivity_analysis/results/`.

The response is the evaluation-verified **Li-6 + Li-7 + Be-9 tritium production per source neutron**, using Li-6 MT105, Li-7 MT52–82 and Be-9 MT700/701. The existing production adapter verifies the installed evaluation and native baseline partial sums before screening. Other material production channels are not assumed to be zero. All material nuclides remain eligible for transport perturbation discovery; the production tally channels are explicit and covariance is retained by reducing scores within each batch.

From `/sandbox/openmc`, with the repository Python environment and `OPENMC_CROSS_SECTIONS` configured:

```bash
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --dry-run
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --run
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --resume
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --reanalyse
```

`--dry-run` discovers and writes the estimated evaluation/history cost without launching Morris or PCE transport. It stops if a configured history cap is exceeded. Review `results/dry_run_plan.json`; increase caps deliberately if the planned campaign is acceptable. `--run` executes verification, low-history Morris screening, selects the top configured count in each family, then runs independent high-history Latin-hypercube training/validation and PCE/Sobol analysis. `--resume` continues checkpoints. `--reanalyse` refits saved transport data only.

The initial settings are a pilot: 200 histories per low-fidelity batch, 20 Morris trajectories; 128 training and 32 validation points at 10,000 histories per batch. They are not a converged full reactor sensitivity study. Discovery may find many parameters; the dry-run estimate and caps prevent an unexpectedly large campaign. Reduce fidelity/counts for a pilot or expand them after reviewing the estimate.

The PCE inputs use hypothetical independent ±20% coefficients. Energy/angular shape modes preserve the existing shared two-mode PDF budget, so one active shape coefficient has a smaller maximum local change than a cross-section coefficient. Sobol indices are conditional on selected parameters and this prior; they are not evaluated covariance uncertainties. The tokamak sensitivity ordering has not yet been calculated.

## Current dry-run estimate

The executed dry run found **2,818 parameters** (114 XS, 138 energy and 2,566 angular modes). It estimates 56,384 Morris evaluations / 112,768,000 histories, 161 high-fidelity evaluations / 16,100,000 histories and 14,000 verification histories: **128,882,000 histories total**. At the saved segmented tokamak rate (~686 histories/s including setup and output), this is roughly 52 hours. Treat that as a planning estimate, since throughput will vary. The sensitivity campaign itself has not been launched. The exact plan is saved in `results/dry_run_plan.json`.
