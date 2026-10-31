# Figure source provenance — 31 October 2026 supervisor-review synchronization

**Scope and accessibility.** This document describes the 14 figure placements in the author's 34-page manuscript revision and distinguishes **(a) figures regenerated from this public repository** from **(b) figures reconstructed using an author-held private research package**. The fact that the manuscript figure lineage has been restored does **not** mean every input, script, or raw molecular-dynamics trajectory is publicly distributed. Figure-generation is also distinct from independently rerunning the full MD simulations or retraining the historical Gaussian-process cross-validation experiments.

The figure-by-figure source mapping is recorded in [`FIGURE_SOURCE_REGISTRY.csv`](FIGURE_SOURCE_REGISTRY.csv). These notes update provenance documentation only; **no simulation, surrogate prediction, or numerical data were changed**.

## Reproduce the figures supported by the public checkout

From the repository root, install the packages in `requirements.txt`, then run:

```bash
python analysis/validate_public_update.py
python analysis/make_public_figures.py --only all
```

The public script has four named functions, each writing a vector PDF to `figures/`:

| Public script slot | Output | Data |
| --- | --- | --- |
| `geometry()` | `fig_pass10_matched_footprints.pdf` | `data/jobs/pass7_geometry_manifest.csv` |
| `orientation_seeds()` | `fig_pass10_pass7_per_seed.pdf` | `data/processed/pass7_seed_results_all12.csv` and `pass7_summary_all4.csv` |
| `pareto22()` | `fig_pareto22_public.pdf` | `data/processed/corrected_22_measured_geometries.csv` |
| `nemd_sensitivity()` | `fig_pass10_nemd_window_sensitivity.pdf` | `data/processed/thermal_runs_51_with_QC.csv` |

**Important:** `fig_pareto22_public.pdf` is a fourth public-data chart; the manuscript currently embeds the separately regenerated **`fig_pass7_pareto22.pdf`**, which is traced to the author-held source package. They represent the same mean-based 22-case Pareto calculation but are separate figure files and visual designs. Do not claim that the public script produces `fig_pass7_pareto22.pdf`.

## Figures requiring the author-held research package

The manuscript's other figure assets were checked/reconstructed on 31 October 2026 using retained geometry builders, corrected numerical summaries, selected raw NEMD/tensile traces, or manuscript-validated tables. Those private generators are **not** part of this lightweight GitHub checkout. In the author-held package:

- `python reproduction/rebuild_figure1.py` reconstructs `fig1_4panel_slices.pdf` using the historical silicon structure-builder and explicit pore coordinates; this is a structural illustration, not a new MD trajectory.
- `python reproduction/regenerate_eight.py` regenerates seven charts/schematics: `fig_pareto.pdf`, `fig_ml_validation.pdf`, `fig_inverse_selection.pdf`, `fig6_inv1_planview.pdf`, `fig_validation.pdf`, `fig_pass7_gp_md.pdf`, and `fig_pass7_pareto22.pdf`. The Pareto classifications are recalculated from sample means (9/17, 9/18, and 10/22, respectively).
- `python build_figures_pass10_private_raw.py` reconstructs thermal profiles, mechanical stress–strain, and stagger/stiffness graphics from selected private MD source records: `fig_pass10_measured_thermal_profiles.pdf`, `fig_pass10_mechanical_stress_strain.pdf`, and `fig_pass10_stagger_stiffness.pdf`.

**Scientific qualifications:** `fig_ml_validation.pdf` reproduces the validated LOOCV/LOFO **tabulated scores**, not independently rerun GP training folds. `fig_pass7_gp_md.pdf` plots frozen pre-MD surrogate predictions versus subsequently measured properties, not an independent refit. Schematic rectangles are pore-footprint projections, not simulated temperature or stress fields. Error bars identified as sample SD are not confidence intervals. The original raw archives and full frozen-prediction provenance are author-held, not uploaded here. These distinctions are relevant to the Methods and Data/Code Availability sections of the manuscript.

## Data integrity and publication status

The core numerical public datasets remain the 22 measured geometry means, 4 prospective orientation summaries and 12 seed-level rows, and original 51 thermal reduction/QC rows. This change only corrects the **figure documentation** after an additional private-source audit; it does not add previously unavailable raw MD records or alter the already published scientific findings.

For a future journal submission, the authors will decide which raw records and supplementary figure generators can be deposited with clear redistribution/licensing terms and a persistent identifier. The manuscript remains under coauthor review.
