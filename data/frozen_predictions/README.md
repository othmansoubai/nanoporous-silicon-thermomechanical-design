# Frozen predictions

Gaussian-process (GP) predictions written to file and checksummed **before**
the molecular-dynamics (MD) simulations of the corresponding cells. The files
are byte-identical copies of the records created on the computing cluster at
freeze time. Each sub-folder verifies against its own original checksum file:

```bash
cd data/frozen_predictions/orientation_cells && sha256sum -c frozen_predictions.sha256
cd data/frozen_predictions/design_cells && sha256sum -c frozen_predictions.sha256
```

| Folder | Cells | Frozen (UTC, from the metadata file) | Model |
|---|---|---|---|
| `orientation_cells/` | ORI35_T_s2, ORI35_W_s2, ORI25_T_s1, ORI25_W_s1 | 2026-10-15T19:58:25.180512+00:00 | GP trained on the 16 porous training cells (repository commit 0de9bd5) |
| `design_cells/` | DES1_5x3_s3, DES2_5x4_s2 | 2026-11-04T01:40:36.169333+00:00 | GP trained on the 21 measured porous cells; batch of two selected by expected hypervolume improvement (repository commit 7f7ad6f) |

Files in each folder:

- `frozen_predictions_preMD.csv`: descriptors, GP mean and predictive standard deviation of kappa_z and E_x
- `geometry_manifest.csv`: pore coordinates (unit cells) and expected atom counts
- `frozen_predictions_metadata.json`: time, code and package versions, training set, settings
- `frozen_predictions.sha256`: checksums written at freeze time
- `design_cells/candidates_all.csv` and `design_cells/design_search_and_freeze.py`: all enumerated candidates and the selection script (both listed in the checksum file)
- `orientation_cells/freeze_predictions.py`: the script that produced the orientation predictions (copied from the same folder; not listed in its checksum file)

The date of the commit that added this folder is its public timestamp.
Measured results of the orientation cells: `data/processed/pass7_summary_all4.csv`
and `pass7_seed_results_all12.csv`.
