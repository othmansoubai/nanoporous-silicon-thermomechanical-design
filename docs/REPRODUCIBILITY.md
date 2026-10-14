# Reproducibility details

## Raw materials and provenance

The 17-geometry κ dataset derives from 51 independently executed LAMMPS 20-bin thermal runs; the mechanical data consist of 51 training traces and 3 INV1 traces. Reproducibility checks on the numerical source material confirmed 394/394 SHA-256 file hashes and numerical reproduction of 51 training κ values and 54 E values. The separately prepared raw-data asset intentionally omits scheduler logs and private paths.

`data/processed/all_18_corrected.csv`: includes 17 training structures and one `held_out_candidate`. Sample SD across three seeds, not confidence bounds. The pristine bulk is included among 17 training geometries but excluded from the 16-porous-model regressions.

The thermal per-seed CSV includes QC columns for the standard temperature window C, alternative window D, energy R², fit R², branch mismatch and temporal sensitivity; no QC flags were reported among the 51 runs. The release does not establish thermal-length convergence or statistical independence of MD trajectories.

## Input source of truth

`inputs/in.nemd_porous_20bin`, `inputs/in.nemd_bulk_20bin_nswap100`, `inputs/in.nemd_inv1`, `inputs/in.uniaxial_tensile_parametric`, and `inputs/Si.sw` record the simulation inputs examined for reproducibility. They use LAMMPS to create the lattice and delete pore atoms. The separate `geometry` Python builder exists for cross-checking and visualization. SLURM submission environment, paths and private cluster logs have intentionally **not** been included in this repository.

For INV1, raw logs on the source HPC reported 24,576 created and 4,733 deleted atoms; the `0%` porosity printout in `in.nemd_inv1` is a variable timing bug, not bulk silicon simulation. Input deck configuration and 20-bin profile outputs were checked for correspondence. The publicly shareable raw data package has only numerical traces; provenance logs should remain in private author archives pending privacy review.

## Scripts

- `analysis/reproduce_raw.py` recomputes all 18 thermal and tensile mean/SD pairs from the 54 raw thermal and 54 raw mechanical files, using the same frozen windows (thermal 1-based bins 3–9 and 13–19; energy fit 100,000–500,000; temperature block averaging 110,000–500,000; mechanical strain 0.0005–0.005).
- `analysis/reproduce_gp.py` recomputes porous16 physical-four GP LOOCV and leave-one-family-out validation, fitting preprocessing independently in each fold. GP kernel: Constant×Matérn ν=2.5 + WhiteKernel; hyperparameter bounds and random seed frozen. Exact optimizer outputs may vary modestly across package versions; reference results are supplied in the prediction CSV. The fitted GP should not be used to claim a prospectively selected INV1, nor should a posterior model comparison be treated as independent test-set accuracy.
- `analysis/reproduce_pareto.py` reports the non-dominated set on the **means** of the 18 verified designs, minimizing κ_z and maximizing E_x.
- `analysis/plot_corrected.py` builds a measured Pareto plot from the 18-row measured dataset.

## Check expected numbers

| Outcome | Expected |
|---|---|
| 18-run mean Pareto set | 9 designs |
| INV1 dominated designs | 6 |
| κ physical GP LOOCV MAE / R² | 0.193635 / 0.959119 |
| E physical GP LOOCV MAE / R² | 1.329916 / 0.977821 |
| κ physical GP family-holdout MAE / R² | 0.188204 / 0.967316 |
| E physical GP family-holdout MAE / R² | 3.055824 / 0.868871 |
| INV1 retrospective physical GP prediction | κ_z=2.556641; E_x=71.663546 |
| Max discrepancy raw-to-processed (this test environment) | 2.84×10⁻¹⁴, floating-point tolerance |

A five-feature baseline model comparison is included as a supplementary benchmark table (`data/processed/ml_legacy5_gp_rf_comparison.csv`). Its input descriptors for INV1 are not uniquely physically defined. Accordingly, no prediction from that baseline is treated as an independent prospective validation of INV1.
