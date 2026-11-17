# OpenMC fork studies

This directory contains project-specific studies and utilities. OpenMC's
runtime implementation and regression tests remain in the repository's normal
`openmc/`, `src/`, `include/`, and `tests/` locations.

## Runtime nuclear-data sensitivity studies

The shared analysis engine is `sandbox/analysis_codes/`. Model definitions and
settings live in `sandbox/studies/<study>/`. From the repository root, for
example:

```bash
OMP_NUM_THREADS=2 python -m sandbox.analysis_codes Li_analysis --dry-run
OMP_NUM_THREADS=2 python -m sandbox.analysis_codes Li_analysis --run
OMP_NUM_THREADS=2 python -m sandbox.analysis_codes Li_analysis --resume
OMP_NUM_THREADS=2 python -m sandbox.analysis_codes Li_analysis --reanalyse
```

Use the same commands with `FLiBe_analysis`, `LiPb_analysis`, `Be_Li_layers`,
`Be_Li_layers_x_beam`, or `tokamak_sensitivity_analysis`. Each study reads its
own `config.py` and `model.py`; generated results stay in its configured local
results directory and are excluded from Git.

The standalone DEMO-inspired CSG application is `sandbox/tokomak_analysis/`.
Run it from the repository root with `python -m sandbox.tokomak_analysis.main
--help`. The scientific report source and analysis examples are in their
corresponding subdirectories.
