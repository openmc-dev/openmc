# FLiBe study inputs and results

The reusable implementation is in [`../analysis_codes/`](../analysis_codes/README.md).
This directory contains the editable study inputs:

- **`model.py`**: fresh OpenMC model builder—materials, geometry and source.
- **`config.py`**: physical settings, Morris/high histories, sample counts,
  seeds, budgets, parameter selection and the PCE backend.
- **`results/`**: all inventories, Morris/PCE/Sobol plots, datasets, checkpoints,
  verification and reports. Historical source files are archived there.

The model is a 20 cm Li₂BeF₄ sphere, density 1.94 g/cm³, 900 K, natural Li,
with a 14 MeV external +z pencil beam and a vacuum world radius of 21 cm.
The response includes verified Li6, Li7, Be9 and F19 tritium-producing partial
reactions, reduced within batches. It does not use cached MT205 as its response.

From `/sandbox/openmc`:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --dry-run
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --run
# Resume interrupted calculations:
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --resume
# Refit/replot saved transport data only:
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --reanalyse
```

Choose `PCE_BACKEND = 'numpy'` or `'openturns'` in `config.py`. Optional dependency:
`/opt/openmc-venv/bin/python -m pip install openturns`. Reanalysis uses the same
saved training/validation data. Backend outputs are preserved separately under
`results/pce_backends/`; `results/high/` shows the latest analysis.

For another model, copy **only `model.py` and `config.py`** into a new study
folder, edit them, and pass that folder to `analysis_codes`. Changing transport
settings requires a fresh results directory; `--output` can provide one.
See the [shared guide](../analysis_codes/README.md) for the public model/response
interface, verified yield limitations and a module-by-module description.

[Summary](results/summary.md) · [Morris/Sobol plots](results/sensitivity_plots.md) ·
[Retained nuclear-data gallery](results/all_retained_plots/index.md) ·
[Selected perturbation gallery](results/nuclear_data_plots/index.md).
