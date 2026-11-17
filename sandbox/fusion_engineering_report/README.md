# Fusion-engineering scientific manuscript

**Read:** `fusion_nuclear_data_sensitivity_report.pdf` (also `main.pdf`).

**Overleaf:** upload `fusion_engineering_overleaf.zip`, select `main.tex` as the main document, and use pdfLaTeX. The ZIP contains all LaTeX sources, bibliography, relative-path figures and tables. It also includes a compact evidence snapshot. No OpenMC installation or nuclear-data library is needed to compile the manuscript.

The report covers the native OpenMC runtime data modifications, cached-MT205 scoring correction, mathematical transport/perturbation/Morris/PCE/Sobol formulation, the five saved benchmark analyses, geometry-dependent rankings, an ordered investigation/measurement review list, and the DEMO-inspired Li₂O tokamak demonstration last. Its limitations are part of the scientific conclusions: hypothetical priors, unequal individual shape/XS amplitudes, three unconverged Morris studies, reduced input sets, additive selected PCE bases, and no tokamak sensitivity campaign.

Author names and affiliations are placeholders for the research team. This is a technical manuscript draft, not an accepted journal paper. Measurement priorities are provisional targets for evaluation/experimental-coverage review, not claims that specific existing measurements are deficient.

## Local compilation

```bash
cd /sandbox/openmc/sandbox/fusion_engineering_report
bash build.sh
```

Requires a standard TeX Live installation with pdfLaTeX, BibTeX and the packages in `main.tex`. `main.bbl` is included in the ZIP for portability, alongside its BibTeX source.

## Reproduce the tables and comparative figures

```bash
cd /sandbox/openmc
/opt/openmc-venv/bin/python sandbox/fusion_engineering_report/prepare_report.py
```

This reads saved results only and never launches transport. It requires NumPy, pandas and Matplotlib. The original study directories must remain available for regeneration; the report's `data/` snapshot separately supports checking its quoted results. The script overwrites only manuscript-generated figures/tables/data copies. It does not overwrite study outputs.

- `data/study_summary.csv`: extracted numerical comparison.
- `data/source_manifest.json`: source artifacts and SHA-256 hashes.
- `data/<study>/`: saved configurations, production verification, full Morris CSV, Sobol CSV, PCE arrays and validation metadata.
- `data/cross_model_sobol.csv`: retained-input comparison; missing means not retained, never an assigned zero.
- `data/tokamak/`: executed model results and engineering-reference records.
- `data/code_snapshot/`: inspected scientific implementation files, included for audit.
- `build/`: compilation logs and checking artifacts; excluded from the Overleaf archive.

No new transport was run to write the report. Existing simulations and native code were preserved. Source-code snapshots identify the current inspection state; saved run provenance remains authoritative for historical computations.
