# Proposed next update to public repository

Target: https://github.com/othmansoubai/nanoporous-silicon-thermomechanical-design

This is an **overlay** for the existing public repository, not a replacement for its original data, history, scripts or LAMMPS inputs.

Added (without overwriting the existing 18-row baseline):

- `data/processed/corrected_22_measured_geometries.csv`
- `data/processed/pareto_front_22_mean_values.csv`
- `data/processed/pass7_seed_results_all12.csv`
- `data/processed/pass7_summary_all4.csv`
- `data/jobs/pass7_geometry_manifest.csv`
- `analysis/make_public_figures.py`
- four reproducible vector PDF charts under `figures/`
- `docs/PASS7_VALIDATION.md`
- `docs/FIGURE_REPRODUCIBILITY.md`
- `docs/PUBLIC_UPDATE_NOTES.md`

**Do not commit:** personal MARWAN directories, compute-node SLURM logs, the full 21 MB raw Pass 7 archive, LAMMPS inputs with unverified third-party redistribution rights, or manuscripts that have not received coauthor approval. The public four-geometry CSVs do not replace the author-held hashes of the original pre-MD prediction freeze.

**No license or publication DOI is asserted.** Reuse permissions, citation metadata and exact target journal must be approved by the authors. README update remains a proposed addition until a versioned Git commit actually contains these files.

**Important chronology:** commit `0de9bd5` is an earlier corrected 18-geometry baseline. Once this overlay is merged to `main`, obtain the new commit SHA and only then update the article's Data and Code Availability text to cite that exact new version; do not reuse the previous 18-geometry-only commit to claim Pass 7 availability.
