# Machine-learning-guided thermomechanical design of nanoporous silicon

Code and data for a study of geometry-controlled thermal conductivity and
mechanical stiffness in nanoporous silicon, combining non-equilibrium molecular
dynamics, uniaxial tensile simulations, Gaussian-process regression, and
multi-objective design analysis.

The work examines **17 training geometries**, including nonporous silicon,
and a separately evaluated inverse-design candidate (**INV1**). The numerical
results presented here were independently reproduced from the simulation
outputs. Thermal conductivity is measured along **z** (κ_z), while Young's
modulus is obtained under loading along **x** (E_x); these are different
directions and should not be treated as equivalent material responses.

**Key result.** INV1 (37.5% stagger, pore aspect ratio 2.5) has
**κ_z = 2.5352 ± 0.0430 W m⁻¹ K⁻¹** and **E_x = 68.394 ± 0.152 GPa**
(three-seed means ± sample standard deviations). The geometry belongs to the
nine-member measured Pareto front over the 18 evaluated structures and
dominates six original designs. It **approaches but does not satisfy** the
original simultaneous target limits κ_z ≤ 2.4 W m⁻¹ K⁻¹ and E_x ≥ 70 GPa.

**Authors:** Othman Soubai, Younes Abouelhanoune, and Mohammed Taibi —
LSA Laboratory, National School of Applied Sciences of Al Hoceima (ENSAH),
Abdelmalek Essaâdi University, Morocco.

> **Repository status:** This public repository provides research data,
> analysis code, molecular-dynamics input files, and reproducibility notes
> for the study described below. Public availability does not imply
> journal acceptance or the assignment of a publication DOI.

---

## What is here

```text
data/processed/  training and combined datasets, QC summaries,
                 descriptors, surrogate-model assessments and Pareto front
data/jobs/       numerical geometry and simulation-configuration manifests
analysis/        scripts for GP validation, Pareto reconstruction,
                 plotting and raw-data reduction
geometry/        independently verified atom-count and geometry builders
inputs/          LAMMPS input decks and the Stillinger–Weber potential file
figures/         measured thermomechanical Pareto front
docs/            reproducibility notes and release checklist
```

### `data/processed/`

| File | Contents |
|---|---|
| `training_17_corrected.csv` | 17 training geometries with κ_z and E_x |
| `all_18_corrected.csv` | 17 training geometries plus INV1, kept identifiable as held out |
| `physical_descriptors_18.csv` | geometry-derived descriptors defined independently of measured responses |
| `thermal_runs_51_with_QC.csv` | per-seed NEMD results and quality checks |
| `mechanical_runs_54.csv` | per-seed tensile results, including INV1 |
| `physical4_family_holdout_predictions.csv` | predictions when each geometry family is excluded in turn |
| `pareto_front_18.csv` | measured non-dominated set, nine structures |

**Important:** INV1 was excluded from fitting and cross-validation of the
17-geometry training dataset. It is a separately evaluated design candidate;
the four-descriptor GP comparison is **retrospective**, not evidence
that this specific surrogate prospectively selected the architecture.

### Numerical methods

- **Thermal:** Müller–Plathe reverse NEMD with the Stillinger–Weber
  potential, a 20-slab temperature profile along *z*, exchange interval of
  100 steps, and 500 ps production. Thermal conductivities are based on
  late-time energy transfer and temperature-gradient fitting. There are
  51 production runs for the 17 training geometries.
- **Mechanical:** initial tensile slope under loading along *x*, fitted over
  engineering strain 0.0005–0.005, with three seeds per geometry.
- **Models:** the main physically interpretable GP uses porosity,
  projected obstruction, pore-orientation descriptor and shift magnitude.
  Its leave-one-family-out MAE is **0.1882 W m⁻¹ K⁻¹** for κ_z and
  **3.0558 GPa** for E_x (16 porous training structures). Model selection
  was exploratory, and transferability outside these geometries is limited.

## Reproduce the analysis

Use Python with NumPy, SciPy, pandas, scikit-learn and Matplotlib:

```bash
python -m pip install -r requirements.txt
python analysis/reproduce_pareto.py
python analysis/reproduce_gp.py
python analysis/plot_corrected.py
```

To reproduce the original numerical reductions from raw MD outputs, obtain
the separately prepared **raw-data ZIP** and extract its `data/raw/` directory
inside this repository's `data/` directory. Then run:

```bash
python analysis/reproduce_raw.py
```

The raw-data package contains temperature profiles, cumulative transfer
energy and stress–strain traces for the 18 designs, with three seeds each.
It is intended for an archival release asset rather than ordinary Git
history. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for
precise fitting windows, GP settings, and known limitations.

## Notes on interpretation

**Atom-counted porosity.** Porosity is measured from atoms removed from the
constructed simulation cell, not estimated from a nominal geometric area.
For INV1 the correct values are 24,576 atoms before deletion and 19,843
after deletion, giving 19.2586% porosity. One archived input prints an
incorrect `0%` because of a logging-variable evaluation defect. The
constructed geometry and simulation outputs are not affected.

**Thermal and mechanical axes differ.** The reported conductivity κ_z and
modulus E_x probe different directions of the nanoporous architecture.

**Model and size limitations.** These results are for finite simulation
cells using the Stillinger–Weber silicon potential. They are not direct
experimental predictions, and do not establish full length convergence.
The small design set and the retrospective comparison for INV1 limit
claims of general-purpose surrogate prediction or prospective optimization.

## Citation

Please cite the corresponding manuscript when published. Until then, cite
the versioned repository release and its commit identifier. A final
publication DOI has not yet been assigned here.

## Contact

Othman Soubai — othman.soubai@etu.uae.ac.ma

LSA Laboratory, ENSAH, Abdelmalek Essaâdi University, Al Hoceima, Morocco.

## Licensing

This repository is publicly accessible, but a license for the authors'
original project code and numerical datasets has not yet been selected.
Public visibility does not itself grant permission to reuse or redistribute
these materials. Licensing terms and redistribution rights for the included
third-party `inputs/Si.sw` parameter file remain under review.
