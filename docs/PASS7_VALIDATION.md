# Pass 7: frozen-prediction orientation validation

**Scope.** Four new porous-silicon geometries in two matched, transverse-rotation pairs; 3 independent thermal NEMD and 3 independent tensile runs per geometry: **12 thermal + 12 mechanical** evaluations. Combined with the 17 original training geometries and historical INV1, the study contains **22 measured geometries**, **66 thermal** and **66 mechanical** seed-level evaluations. These totals are not 132 independent geometries.

**Directional properties.** `kappa_z` is the z-axis, finite-cell, classical Stillinger–Weber silicon conductivity measured by Müller–Plathe NEMD. `E_x` is the small-strain modulus along x from tensile MD, at approximately 300 K. These are different axes. Reported `±` values are sample standard deviations over 3 runs, not measurement calibration or 95% confidence intervals.

## Verified results

| Geometry | Atom-counted porosity | Mean κ_z ± SD (W m⁻¹ K⁻¹) | Mean E_x ± SD (GPa) |
|---|---:|---:|---:|
| ORI35_T_s2 | 28.01% | 1.92438 ± 0.02800 | 42.2142 ± 2.3599 |
| ORI35_W_s2 | 28.01% | 1.91627 ± 0.05835 | 59.5298 ± 1.3454 |
| ORI25_T_s1 | 19.18% | 3.86323 ± 0.11606 | 46.2700 ± 2.3911 |
| ORI25_W_s1 | 19.18% | 3.66055 ± 0.21010 | 70.8799 ± 4.0126 |

Wide minus tall: ΔE_x = +17.316 GPa (+41.02%) at 28.01% porosity and +24.610 GPa (+53.19%) at 19.18%; Δκ_z = −0.00811 W m⁻¹ K⁻¹ (−0.42%) and −0.20268 W m⁻¹ K⁻¹ (−5.25%), respectively. Thermal differences are unresolved at n=3; this **does not demonstrate equivalence** within an a priori margin. Two pair designs differ in more than porosity; do not infer a porosity-only law from these contrasts.

## Frozen GP prospective prediction

The four pre-MD GP predictions were frozen on MARWAN and verified before running Pass 7. The model was **not refitted** using these four structures before the comparisons.

| Geometry | Pre-MD GP κ_z | Pre-MD GP E_x |
|---|---:|---:|
| ORI35_T_s2 | 1.844835 | 45.830157 |
| ORI35_W_s2 | 1.874637 | 57.398456 |
| ORI25_T_s1 | 3.579609 | 55.146336 |
| ORI25_W_s1 | 3.513095 | 74.716047 |

Prospective mean absolute error over the 4 geometry means: **κ_z = 0.13806195 W/(m K)**, **E_x = 4.6149492 GPa**. This evaluates property prediction on *controlled unseen structures*. It is not evidence that the revised GP prospectively selected a new inverse-optimal design. Historical INV1 selection was performed using an earlier descriptor/model formulation and must be described separately.

## Archive integrity

The source Pass 7 numeric CSVs were checked against the MARWAN-recorded SHA-256 hashes:

```
380077627a44a320e4dd365b82729d6561a108a13240a8a4c84ee4b066025313  pass7_seed_results_all12.csv
72c343d6bcda838626f592588e660eb3f52b671930ce4e70897c09f3b7e278c9  pass7_summary_all4.csv
```

The complete original Pass 7 raw LAMMPS output bundle `PASS7_ORIENTATION_COMPLETE.tar.gz` is **author-held**, has not been added to this public repository, and must not be asserted publicly released. The original pre-MD freeze and checksums remain in that author-held bundle. The published summaries include the frozen predictions, but the archived evidence for their creation order is distinct from the CSVs.

## Numerical caveats

All 12 thermal runs achieved strong fitted energy transfer (reported R² 0.999973–0.999997), temperature-gradient branch R² ≥ 0.993339, and branch mismatches ≤ 6.23%. Mechanical moduli were fitted with the same predefined engineering-strain window 0.0005–0.005. Tall-pore raw fits were noisier (minimum reported R² 0.7845 for ORI35 and 0.8219 for ORI25), although alternative fit windows were reasonably stable. These checks do not establish transport-length convergence or all microscopic phonon-scattering mechanisms.

Pareto classification over 22 structures uses mean κ_z and E_x, yielding ten non-dominated structures, **without** propagation of small-sample uncertainty into dominance probabilities.
