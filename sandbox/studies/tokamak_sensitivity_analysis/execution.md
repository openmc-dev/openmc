# Execution commands

Run from the repository root with the modified OpenMC Python environment:

```bash
cd /sandbox/openmc
export OPENMC_CROSS_SECTIONS=/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=2

# Recalculate parameter inventory and cost after editing config/model:
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --dry-run

# Start verification, Morris, PCE training/validation and Sobol analysis:
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --run

# Continue saved checkpoints after interruption:
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --resume

# Refit/plot saved results without transport:
python -m sandbox.analysis_codes tokamak_sensitivity_analysis --reanalyse
```

The executed dry-run plans 128,882,000 histories, including 56,384 Morris evaluations. Review this cost before `--run`. Adjust histories, trajectories, training count, fidelity and caps in `config.py`; then rerun `--dry-run`. Geometry/material/source settings live in `TOKAMAK_CONFIG` in the same file, while detailed constructors remain in `tokomak_analysis/`. Save each changed configuration to its own results directory using the CLI `--output PATH` option.
