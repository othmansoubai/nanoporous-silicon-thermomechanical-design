# Nanoporous silicon thermomechanical design

**Current research dataset: 22 measured geometries, with prospective surrogate validation.** This public repository contains selected numerical data and code for molecular-dynamics-based geometry control of thermal transport and directional stiffness in nanoporous silicon, alongside exploratory Gaussian-process (GP) modeling and multi-objective design analysis.

The combined study comprises **17 original training geometries** (including pristine bulk silicon), the historically inverse-designed **INV1** candidate, and **four new orientation-validation structures**. There are three independently seeded thermal and mechanical MD runs per geometry: **66 NEMD runs and 66 tensile runs** overall.

**Transport and loading directions differ:** thermal conductivity is measured in the **z** direction (κ_z) and elastic modulus in the **x** direction (E_x). The data reflect finite simulation cells and the classical Stillinger–Weber interatomic potential; they are not experimentally validated or length-converged bulk silicon properties.

## Main results

- **INV1 (historical inverse-design case):** corrected κ_z = **2.5352 ± 0.0430 W/(m·K)** and E_x = **68.394 ± 0.152 GPa**, mean ± sample SD from three seeds. INV1 meets **neither** original hard simultaneous target (κ_z ≤ 2.4 W/(m·K) and E_x ≥ 70 GPa).
- **Prospective GP property validation:** Four new pore geometries were evaluated in two equal-porosity rotation pairs *after* GP predictions were frozen. The wide structures show approximately **41.0%** and **53.2%** higher measured E_x than their tall counterparts. Corresponding κ_z differences in the sample means are **−0.42%** and **−5.25%**. With three seeds, these comparisons do **not** establish formal conductivity equivalence.
- **Predictive limits:** The frozen GP correctly predicts the direction of the mechanical orientation contrast but underestimates its size. The four new geometries validate *property prediction*, not a newly GP-optimized inverse-design structure. INV1 is historical; comparing it to the updated four-descriptor GP is retrospective.
- **Multi-objective summary:** the 22 measured means contain **10 mean-value non-dominated structures**, with Pareto classification based on means rather than statistically propagated dominance probabilities.

## Data: current version and provenance

| File | Meaning |
| --- | --- |
| [`data/processed/corrected_22_measured_geometries.csv`](data/processed/corrected_22_measured_geometries.csv) | **Current** 22-structure mean κ_z and E_x dataset |
| [`data/processed/pareto_front_22_mean_values.csv`](data/processed/pareto_front_22_mean_values.csv) | **Current** 10-member measured mean-value front |
| [`data/processed/pass7_seed_results_all12.csv`](data/processed/pass7_seed_results_all12.csv) | Individual results for four new geometries × three seeds per property |
| [`data/processed/pass7_summary_all4.csv`](data/processed/pass7_summary_all4.csv) | Means, sample SDs, frozen GP predictions and signed relative errors |
| [`data/jobs/pass7_geometry_manifest.csv`](data/jobs/pass7_geometry_manifest.csv) | Construction coordinates for the four prospective orientation geometries |
| [`data/processed/all_18_corrected.csv`](data/processed/all_18_corrected.csv) | **Historical** 17 training structures plus held-out INV1 |
| [`data/processed/training_17_corrected.csv`](data/processed/training_17_corrected.csv) | Baseline training dataset |
| [`data/processed/thermal_runs_51_with_QC.csv`](data/processed/thermal_runs_51_with_QC.csv) | Original 51 NEMD seed observations and analysis diagnostics |
| [`data/processed/mechanical_runs_54.csv`](data/processed/mechanical_runs_54.csv) | Original 54 tensile seed results including INV1 |
| [`data/processed/pareto_front_18.csv`](data/processed/pareto_front_18.csv) | **Historical** nine-member front for only 18 structures |

Older inputs and analysis scripts are kept to preserve scientific provenance. The repository does **not** claim that earlier 18-geometry scripts already regenerate the updated 22-geometry figures. [`SHA256SUMS_RELEASE.txt`](SHA256SUMS_RELEASE.txt) records the historical repository's previous release content and should not be treated as a checksum manifest for this new 22-geometry revision.

## Reproduce public-data results and figures

Recommended Python 3.11 with packages in [`requirements.txt`](requirements.txt):

```bash
python -m pip install -r requirements.txt
python analysis/validate_public_update.py
python analysis/make_public_figures.py --only all
```

Each new public-data figure has its own function/slot inside [`analysis/make_public_figures.py`](analysis/make_public_figures.py):

1. `geometry()` → `fig_pass10_matched_footprints.pdf`
2. `orientation_seeds()` → `fig_pass10_pass7_per_seed.pdf`
3. `pareto22()` → `fig_pareto22_public.pdf`
4. `nemd_sensitivity()` → `fig_pass10_nemd_window_sensitivity.pdf`

See [`docs/FIGURE_REPRODUCIBILITY.md`](docs/FIGURE_REPRODUCIBILITY.md) and [`docs/PASS7_VALIDATION.md`](docs/PASS7_VALIDATION.md) for mappings and scientific qualifications. Earlier sources and their analysis remain in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Raw simulation data and limitations

The original numerical-reduction archive and full Pass 7 raw MD trajectories, frozen pre-MD predictions, SLURM logs, and unpublished manuscript files are maintained by the authors **outside this public GitHub repository**. Consequently, plots dependent on those raw files cannot be fully regenerated from this lightweight public repository alone. Public checksums confirm the integrity of the two Pass 7 numerical CSVs, but do not imply all raw MD inputs are publicly distributed.

Three-seed sample standard deviations quantify seed-to-seed variation; they are not confidence bounds or a complete uncertainty budget. The finite MD cells, system length, limited geometry families, and preliminary GP transfer limit the generality of conclusions.

## Citation, authors and licensing

**Manuscript status:** under author review, not yet accepted. Until a journal publication or curated DOI is available, reference the **exact Git commit** used for the results. No DOI has been assigned to this revision.

Othman Soubai, Younes Abouelhanoune, Mohammed Taibi — LSA Laboratory, ENSAH, Abdelmalek Essaâdi University, Morocco. Corresponding author: othman.soubai@etu.uae.ac.ma.

No new license has been selected for the original code and numerical data; public access does not itself grant broad redistribution rights. The redistribution terms for third-party `inputs/Si.sw` still require separate confirmation.
