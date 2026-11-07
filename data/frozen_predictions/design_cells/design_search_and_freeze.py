#!/usr/bin/env python3
# =============================================================================
# design_search_and_freeze.py
# -----------------------------------------------------------------------------
# Prospective design campaign for the nanoporous-silicon kappa_z / E_x paper.
#
# What this script does, in order:
#   1. SELF-TEST A  - re-counts atoms for every geometry already simulated
#                     (16 porous training cells, INV1, 4 orientation cells)
#                     with the same inclusive-box rule LAMMPS uses, and checks
#                     the counts against the published porosities.
#   2. SELF-TEST B  - recomputes the four coordinate-derived GP features from
#                     the pore coordinates and checks them against
#                     data/processed/physical_descriptors_18.csv.
#   3. SELF-TEST C  - refits the published four-feature GP on the original
#                     16 porous cells and checks that it reproduces the four
#                     predictions frozen before the orientation MD (Pass 7).
#                     This proves the Python environment gives the same GP.
#   4. DESIGN MODEL - refits the same GP on ALL 21 measured porous cells
#                     (16 training + INV1 + 4 orientation cells).
#   5. SEARCH       - enumerates every buildable pair of identical rectangular
#                     footprints with a diagonal stagger in the 8 x 8 cell,
#                     keeps only candidates inside the training feature range
#                     and not already simulated (nor their x<->y transpose).
#   6. DIAGNOSTIC   - reports the best GP probability, over all candidates, of
#                     reaching the ORIGINAL target box (kappa <= 2.4, E >= 70).
#   7. SELECTION    - picks a batch of two candidates by Monte-Carlo EXPECTED
#                     HYPERVOLUME IMPROVEMENT (EHVI) of the measured 22-cell
#                     mean Pareto front (minimise kappa_z, maximise E_x). The
#                     second pick is chosen after adding the first pick's GP
#                     mean to the front ("Kriging believer"), so the two
#                     candidates target different parts of the front.
#   8. FREEZE       - only with --freeze: writes the 2-row geometry manifest,
#                     the frozen predictions, metadata and SHA-256 checksums.
#                     Refuses to overwrite an existing freeze.
#
# Run WITHOUT --freeze first (self-tests + preview, nothing frozen).
# Run WITH --freeze once, BEFORE any new MD is submitted.
#
# Usage (on MARWAN, from ~/phase3_mech/design_campaign):
#   python design_search_and_freeze.py ~/nanoporous-silicon-thermomechanical-design
#   python design_search_and_freeze.py ~/nanoporous-silicon-thermomechanical-design --freeze
# =============================================================================

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import norm

# -----------------------------------------------------------------------------
# DESIGN SETTINGS - part of the frozen record (the script's SHA-256 is stored
# with the predictions). Fix them before running --freeze; never change them
# after any new MD result exists.
# -----------------------------------------------------------------------------
# Original target box of the earlier design workflow (diagnostic only).
ORIGINAL_BOX = {"kappa_max": 2.4, "E_min": 70.0}
# Reference point of the hypervolume = the design region of interest.
# Thermal-insulation regime: kappa_z <= 0.40 x kappa_pristine (finite cell),
# i.e. at least a 60 % reduction. Without this cap, the hypervolume is
# dominated by the stiff, highly conducting end of the front, which is not
# the application. E reference 40 GPa lies below every measured cell
# (lowest measured E_x = 42.11 GPa for A_S4.5_d2).
REF_KAPPA_FRACTION = 0.40
REF_E = 40.0
N_BATCH = 2            # number of candidates sent to MD
N_MC = 20000           # Monte-Carlo samples per candidate for EHVI
MC_SEED = 20261104     # fixed RNG seed -> identical EHVI on every rerun

# Repository commit whose data and GP code this campaign is built on.
REQUIRED_COMMIT = "7f7ad6ffc6a0fbaae38655408556a61bddc1dd2c"

# Lattice constants of the production cell (8 x 8 x 48 conventional cells).
NX, NY, NZ = 8, 8, 48
N0 = 8 * NX * NY * NZ                       # 24 576 atoms in the pristine cell
# Fractional coordinates of the 8 atoms in the diamond-cubic conventional cell.
BASIS = np.array([(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0),
                  (.25, .25, .25), (.25, .75, .75), (.75, .25, .75), (.75, .75, .25)])

# Candidate search space: identical A and B footprints dx x dy (unit cells),
# B shifted diagonally by s cells, as for every measured staggered family.
DX_RANGE = range(1, 7)
DY_RANGE = range(1, 7)
S_RANGE = range(0, 5)

FEATS = ["phi_actual", "projected_obstruction", "orientation_chi", "shift_mag_uc"]
CAMPAIGN = Path(__file__).resolve().parent


# =============================================================================
# Geometry helpers
# =============================================================================
def lattice_positions():
    """All atom positions of the pristine cell in lattice units, shape (N0, 3)."""
    ix, iy, iz = np.meshgrid(np.arange(NX), np.arange(NY), np.arange(NZ), indexing="ij")
    cells = np.stack([ix.ravel(), iy.ravel(), iz.ravel()], axis=1)
    pos = (cells[:, None, :] + BASIS[None, :, :]).reshape(-1, 3)
    assert pos.shape[0] == N0
    return pos


POS = lattice_positions()


def remaining_atoms(b):
    """Atoms left after carving, mirroring the LAMMPS input exactly.

    b = (Axlo, Axhi, Aylo, Ayhi, Bxlo, Bxhi, Bylo, Byhi) in lattice units.
    LAMMPS 'region block' is INCLUSIVE on all faces, pore A occupies the
    z-slices [8k, 8k+4] and pore B the slices [8k+4, 8k+8], k = 0..5.
    """
    x, y, z = POS[:, 0], POS[:, 1], POS[:, 2]
    in_A_xy = (x >= b[0]) & (x <= b[1]) & (y >= b[2]) & (y <= b[3])
    in_B_xy = (x >= b[4]) & (x <= b[5]) & (y >= b[6]) & (y <= b[7])
    in_A_z = np.zeros_like(z, dtype=bool)
    in_B_z = np.zeros_like(z, dtype=bool)
    for k in range(0, NZ, 8):
        in_A_z |= (z >= k) & (z <= k + 4)
        in_B_z |= (z >= k + 4) & (z <= k + 8)
    deleted = (in_A_xy & in_A_z) | (in_B_xy & in_B_z)
    return int(N0 - deleted.sum())


def features(b, n_remaining):
    """Four coordinate-derived GP features, same convention as the release.

    phi is rounded to 2 decimals, exactly as in physical_descriptors_18.csv
    and in the Pass-7 freeze, so that new and old cells are encoded alike.
    """
    dx, dy = b[1] - b[0], b[3] - b[2]
    assert abs(dx - (b[5] - b[4])) < 1e-12 and abs(dy - (b[7] - b[6])) < 1e-12
    ox = max(0.0, min(b[1], b[5]) - max(b[0], b[4]))       # x-overlap of A and B
    oy = max(0.0, min(b[3], b[7]) - max(b[2], b[6]))       # y-overlap of A and B
    obstruction = (2 * dx * dy - ox * oy) / (NX * NY)      # |A u B| / 64
    chi = (dx - dy) / (dx + dy) if (dx + dy) > 0 else 0.0  # in-plane orientation
    shift = math.hypot(b[4] - b[0], b[6] - b[2])           # displacement of B
    phi = round(100.0 * (N0 - n_remaining) / N0, 2)
    return [phi, obstruction, chi, shift]


def footprint_key(b):
    """Canonical (dx, dy, sx, sy) key used to detect already-simulated cells."""
    return (round(b[1] - b[0], 3), round(b[3] - b[2], 3),
            round(b[4] - b[0], 3), round(b[6] - b[2], 3))


def transpose_key(k):
    """Key of the x<->y transposed geometry (same kappa_z by symmetry)."""
    return (k[1], k[0], k[3], k[2])


# =============================================================================
# Pareto / hypervolume helpers (minimise kappa, maximise E)
# =============================================================================
def pareto_front(points):
    """Non-dominated subset of (kappa, E) pairs, sorted by increasing kappa."""
    pts = np.asarray(points, float)
    keep = []
    for i, (k, e) in enumerate(pts):
        dominated = np.any((pts[:, 0] <= k) & (pts[:, 1] >= e)
                           & ((pts[:, 0] < k) | (pts[:, 1] > e)))
        if not dominated:
            keep.append((k, e))
    return np.array(sorted(keep))


def clip_to_roi(front, ref_k, ref_e):
    """Keep only front points inside the region of interest (k < ref_k, E > ref_e)."""
    m = (front[:, 0] < ref_k) & (front[:, 1] > ref_e)
    return front[m]


def hv_improvement(front, k, e, ref_k, ref_e):
    """Hypervolume gained by adding points (k, e) to a 2-D front, vectorised.

    Objectives are converted to 'minimise both' with f1 = kappa, f2 = -E.
    The front's attained boundary h(x) is a staircase: at abscissa x it is the
    smallest f2 among front points with f1 <= x (or the reference -ref_e).
    A new point (a, b) adds the area between b and h(x) for x in [a, ref_k]:
        HVI = integral_a^ref_k max(0, h(x) - b) dx.
    `front` must already be clipped to the region of interest.
    """
    a = np.minimum(np.asarray(k, float), ref_k)      # samples beyond ref add nothing
    b = -np.asarray(e, float)
    f1 = front[:, 0]
    f2 = -front[:, 1]
    xs = np.concatenate([f1, [ref_k]])               # staircase breakpoints
    hv = np.zeros_like(a)
    for j in range(len(xs) - 1):
        # on [xs[j], xs[j+1]] the boundary height is min(f2[:j+1])
        h = min(f2[: j + 1].min(), -ref_e)
        lo = np.maximum(a, xs[j])
        width = np.clip(xs[j + 1] - lo, 0, None)
        hv += width * np.clip(h - b, 0, None)
    # left of the first front point the boundary is the reference itself
    width0 = np.clip((xs[0] if len(f1) else ref_k) - a, 0, None)
    hv += width0 * np.clip(-ref_e - b, 0, None)
    return hv


# =============================================================================
# Data loading
# =============================================================================
def read_csv(path):
    with Path(path).open(newline="", encoding="utf8") as f:
        return list(csv.DictReader(f))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repo", help="path to nanoporous-silicon-thermomechanical-design")
    ap.add_argument("--freeze", action="store_true",
                    help="write the frozen manifest and predictions (run once, before MD)")
    args = ap.parse_args()
    repo = Path(args.repo).expanduser().resolve()

    # ---------------------------------------------------------------- commit
    head = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    if head != REQUIRED_COMMIT:
        sys.exit(f"STOP: repository HEAD is {head}, expected {REQUIRED_COMMIT}.")
    print(f"[ok] repository at commit {head[:7]}")

    # Use the published GP definition itself, so the configuration is identical.
    sys.path.insert(0, str(repo / "analysis"))
    from reproduce_gp import model                      # noqa: E402
    from sklearn.exceptions import ConvergenceWarning   # noqa: E402
    import sklearn                                      # noqa: E402
    warnings.simplefilter("ignore", ConvergenceWarning)

    # ---------------------------------------------------------- input tables
    jobs = read_csv(repo / "data/jobs/thermal_jobs_51.csv")
    desc = {r["label"]: r for r in read_csv(repo / "data/processed/physical_descriptors_18.csv")}
    meas = {r["label"]: r for r in read_csv(repo / "data/processed/all_18_corrected.csv")}
    p7_manifest = read_csv(repo / "data/jobs/pass7_geometry_manifest.csv")
    p7_summary = {r["label"]: r for r in read_csv(repo / "data/processed/pass7_summary_all4.csv")}

    # Pore coordinates of every measured porous geometry.
    bounds = {}
    for r in jobs:
        if r["label"] == "bulk_si":
            continue
        bounds[r["label"]] = tuple(float(r[k]) for k in
                                   ("pAxlo", "pAxhi", "pAylo", "pAyhi", "pBxlo", "pBxhi", "pBylo", "pByhi"))
    # INV1 coordinates as stated in the manuscript: A = x[0,5] y[0,2], B = x[3,8] y[3,5].
    bounds["INV1_AR2.5_s37"] = (0, 5, 0, 2, 3, 8, 3, 5)
    for r in p7_manifest:
        bounds[r["label"]] = tuple(float(r[k]) for k in
                                   ("pAxlo", "pAxhi", "pAylo", "pAyhi", "pBxlo", "pBxhi", "pBylo", "pByhi"))
    training16 = [l for l in bounds if l in meas and meas[l]["role"] == "training"]
    assert len(training16) == 16, training16
    ori4 = [r["label"] for r in p7_manifest]

    # ------------------------------------------------------------ SELF-TEST A
    counts = {l: remaining_atoms(b) for l, b in bounds.items()}
    for l in training16:
        phi = 100 * (N0 - counts[l]) / N0
        assert abs(phi - float(desc[l]["phi_actual"])) < 0.006, (l, phi, desc[l]["phi_actual"])
    assert counts["INV1_AR2.5_s37"] == 19843, counts["INV1_AR2.5_s37"]
    for r in p7_manifest:
        assert counts[r["label"]] == int(r["expected_remaining_atoms"]), r["label"]
    print("[ok] SELF-TEST A: atom counts reproduce all 21 simulated porous cells")

    # ------------------------------------------------------------ SELF-TEST B
    feat = {l: features(b, counts[l]) for l, b in bounds.items()}
    for l in training16 + ["INV1_AR2.5_s37"]:
        ref = [float(desc[l][f]) for f in FEATS]
        assert abs(feat[l][0] - ref[0]) < 0.006, (l, "phi", feat[l][0], ref[0])
        assert np.allclose(feat[l][1:], ref[1:], atol=1e-9), (l, feat[l], ref)
    print("[ok] SELF-TEST B: coordinate features reproduce physical_descriptors_18.csv")

    # ------------------------------------------------------------ SELF-TEST C
    X16 = np.array([[float(desc[l][f]) for f in FEATS] for l in training16])
    Xori = np.array([feat[l] for l in ori4])
    for key, col in (("kappa_z_mean_W_mK", "GP_kappa"), ("E_x_mean_GPa", "GP_E")):
        y16 = np.array([float(meas[l][key]) for l in training16])
        gp = model().fit(X16, y16)
        pred = gp.predict(Xori)
        frozen = np.array([float(p7_summary[l][col]) for l in ori4])
        assert np.allclose(pred, frozen, rtol=1e-6), (key, pred, frozen)
    print("[ok] SELF-TEST C: GP refit reproduces the 4 predictions frozen before Pass 7")

    # ----------------------------------------------------------- DESIGN MODEL
    # Train on every measured porous cell (16 + INV1 + 4 orientation cells).
    train_labels = training16 + ["INV1_AR2.5_s37"] + ori4
    Xtr = np.array([feat[l] if l in ori4 else [float(desc[l][f]) for f in FEATS]
                    for l in train_labels])
    y = {"kappa": [], "E": []}
    for l in train_labels:
        if l in ori4:
            y["kappa"].append(float(p7_summary[l]["kappa_mean"]))
            y["E"].append(float(p7_summary[l]["E_mean"]))
        else:
            y["kappa"].append(float(meas[l]["kappa_z_mean_W_mK"]))
            y["E"].append(float(meas[l]["E_x_mean_GPa"]))
    gps = {k: model().fit(Xtr, np.array(v)) for k, v in y.items()}
    print(f"[ok] design GP fitted on {len(train_labels)} measured porous cells")

    # ------------------------------------------------------------------ SEARCH
    measured_keys = set()
    for b in bounds.values():
        k = footprint_key(b)
        measured_keys |= {k, transpose_key(k)}
    lo = Xtr.min(axis=0)
    hi = Xtr.max(axis=0)
    chi_max = np.abs(Xtr[:, 2]).max()

    cands = []
    for dx in DX_RANGE:
        for dy in DY_RANGE:
            for s in S_RANGE:
                if s + dx > NX or s + dy > NY:          # B must stay inside the box (no wrap)
                    continue
                b = (0, dx, 0, dy, s, s + dx, s, s + dy)
                n = remaining_atoms(b)
                f = features(b, n)
                in_range = (lo[0] <= f[0] <= hi[0] and lo[1] <= f[1] <= hi[1]
                            and abs(f[2]) <= chi_max + 1e-12 and f[3] <= hi[3] + 1e-12)
                already = footprint_key(b) in measured_keys
                cands.append(dict(dx=dx, dy=dy, s=s, bounds=b, N=n, feat=f,
                                  in_range=in_range, already=already))
    pool = [c for c in cands if c["in_range"] and not c["already"]]
    Xc = np.array([c["feat"] for c in pool])
    mk, sk = gps["kappa"].predict(Xc, return_std=True)
    mE, sE = gps["E"].predict(Xc, return_std=True)
    for i, c in enumerate(pool):
        c.update(mk=mk[i], sk=sk[i], mE=mE[i], sE=sE[i])
        # Probability of meeting the ORIGINAL box, treating the two GPs as
        # independent Gaussians. GP SDs are not calibrated: a ranking score only.
        c["P_box"] = (norm.cdf((ORIGINAL_BOX["kappa_max"] - c["mk"]) / c["sk"])
                      * norm.cdf((c["mE"] - ORIGINAL_BOX["E_min"]) / c["sE"]))
    print(f"[ok] {len(cands)} buildable candidates, {len(pool)} in range and not yet simulated")

    # -------------------------------------------------------------- DIAGNOSTIC
    best_box = max(pool, key=lambda c: c["P_box"])
    print(f"\nDiagnostic - original box kappa_z <= {ORIGINAL_BOX['kappa_max']}, "
          f"E_x >= {ORIGINAL_BOX['E_min']}: best candidate score = {best_box['P_box']:.4f} "
          f"({best_box['dx']}x{best_box['dy']}, s={best_box['s']})")

    # --------------------------------------------------------------- SELECTION
    # Measured 22-cell means (17 training incl. pristine, INV1, 4 orientation).
    m22 = read_csv(repo / "data/processed/corrected_22_measured_geometries.csv")
    assert len(m22) == 22
    kappa_pristine = float(next(r["kappa"] for r in m22 if r["label"] == "bulk_si"))
    ref_k = REF_KAPPA_FRACTION * kappa_pristine
    front_all = pareto_front([(float(r["kappa"]), float(r["E"])) for r in m22])
    front = clip_to_roi(front_all, ref_k, REF_E)
    print(f"[ok] measured mean front: {len(front_all)} non-dominated cells, "
          f"{len(front)} inside the insulating region kappa_z < {ref_k:.3f}")

    # ------------------------------------------------------------ SELF-TEST D
    # Vectorised improvement must equal HV(front + p) - HV(front), where the
    # hypervolume of an ROI-clipped front sorted by kappa (E rising along it) is
    #   HV = sum_i (k_{i+1} - k_i) * (E_i - REF_E),   k_{n+1} = ref_k.
    def hv_exact(pts):
        f = clip_to_roi(pareto_front(pts), ref_k, REF_E)
        if len(f) == 0:
            return 0.0
        k_next = np.append(f[1:, 0], ref_k)
        return float(np.sum((k_next - f[:, 0]) * (f[:, 1] - REF_E)))
    test_rng = np.random.default_rng(1)
    for _ in range(500):
        k, e = test_rng.uniform(0.5, 11.0), test_rng.uniform(40.5, 99.0)
        direct = hv_exact(np.vstack([front, [k, e]])) - hv_exact(front)
        assert abs(direct - hv_improvement(front, [k], [e], ref_k, REF_E)[0]) < 1e-9, (k, e)
    print("[ok] SELF-TEST D: vectorised hypervolume improvement matches direct computation")

    rng = np.random.default_rng(MC_SEED)
    z = rng.standard_normal((N_MC, 2))               # common random numbers
    chosen, taken = [], set()
    for pick in range(N_BATCH):
        for c in pool:
            ks = c["mk"] + c["sk"] * z[:, 0]
            es = c["mE"] + c["sE"] * z[:, 1]
            c[f"EHVI_{pick + 1}"] = hv_improvement(front, ks, es, ref_k, REF_E).mean()
        ranked = sorted((c for c in pool if footprint_key(c["bounds"]) not in taken),
                        key=lambda c: -c[f"EHVI_{pick + 1}"])
        best = ranked[0]
        chosen.append(best)
        kb = footprint_key(best["bounds"])
        taken |= {kb, transpose_key(kb)}
        print(f"\nPick {pick + 1} (top 5 by EHVI, first selected)")
        print("   dx dy  s   phi%   O_proj   chi    shift | kappa_GP +/- sd | E_GP +/- sd |  EHVI")
        for c in ranked[:5]:
            f = c["feat"]
            print(f"   {c['dx']:2d} {c['dy']:2d} {c['s']:2d}  {f[0]:5.2f}  {f[1]:6.4f}  {f[2]:+.3f}  "
                  f"{f[3]:5.3f} |  {c['mk']:5.3f} +/- {c['sk']:4.3f} | {c['mE']:5.1f} +/- "
                  f"{c['sE']:3.1f} | {c[f'EHVI_{pick + 1}']:.4f}")
        # Kriging believer: pretend the pick lands on its GP mean, then re-score
        front = clip_to_roi(pareto_front(np.vstack([front, [best["mk"], best["mE"]]])),
                            ref_k, REF_E)

    # Full candidate table, always written (it is a preview, not a frozen record).
    with (CAMPAIGN / "candidates_all.csv").open("w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh)
        ehvi_cols = [f"EHVI_{p + 1}" for p in range(N_BATCH)]
        w.writerow(["dx", "dy", "s", "remaining_atoms", *FEATS, "in_training_range",
                    "already_simulated", "kappa_GP", "kappa_GP_sd", "E_GP", "E_GP_sd",
                    "P_original_box", *ehvi_cols])
        for c in cands:
            row = [c["dx"], c["dy"], c["s"], c["N"], *c["feat"], c["in_range"], c["already"]]
            row += ([c["mk"], c["sk"], c["mE"], c["sE"], c["P_box"], *[c.get(k, "") for k in ehvi_cols]]
                    if "mk" in c else [""] * (5 + N_BATCH))
            w.writerow(row)
    print("\n[ok] wrote candidates_all.csv (preview table)")

    if not args.freeze:
        print("\nPREVIEW ONLY - nothing frozen. Re-run with --freeze before submitting any MD.")
        return

    # ------------------------------------------------------------------ FREEZE
    out_pred = CAMPAIGN / "frozen_predictions_preMD.csv"
    out_man = CAMPAIGN / "geometry_manifest.csv"
    for p in (out_pred, out_man):
        if p.exists():
            sys.exit(f"STOP: {p.name} already exists - refusing to replace a frozen record.")

    with out_man.open("w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "pAxlo", "pAxhi", "pAylo", "pAyhi", "pBxlo", "pBxhi", "pBylo",
                    "pByhi", "expected_remaining_atoms", "expected_porosity_pct", "phase"])
        for p, c in enumerate(chosen, start=1):
            b = c["bounds"]
            label = f"DES{p}_{c['dx']}x{c['dy']}_s{c['s']}"
            c["label"] = label
            w.writerow([label, *b, c["N"], 100 * (N0 - c["N"]) / N0, f"design_pick{p}"])

    with out_pred.open("w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "pick", *FEATS, "GP_pred_kappa_z_W_mK", "GP_model_sd_kappa_z_W_mK",
                    "GP_pred_E_x_GPa", "GP_model_sd_E_x_GPa", "EHVI_at_selection",
                    "P_original_box"])
        for p, c in enumerate(chosen, start=1):
            w.writerow([c["label"], p, *c["feat"], c["mk"], c["sk"], c["mE"], c["sE"],
                        c[f"EHVI_{p}"], c["P_box"]])

    meta = {
        "utc_generated": datetime.now(timezone.utc).isoformat(),
        "repository_commit": head,
        "sklearn_version": sklearn.__version__,
        "numpy_version": np.__version__,
        "training_labels": train_labels,
        "features": FEATS,
        "original_box_diagnostic": ORIGINAL_BOX,
        "best_original_box_score": best_box["P_box"],
        "hypervolume_reference": {"kappa": ref_k, "kappa_fraction_of_pristine": REF_KAPPA_FRACTION,
                                  "E": REF_E},
        "monte_carlo": {"samples": N_MC, "seed": MC_SEED},
        "selection_rule": "batch of 2 by Monte-Carlo expected hypervolume improvement of the "
                          "measured 22-cell mean front (min kappa_z, max E_x), Kriging-believer "
                          "update between picks; candidates restricted to the training feature "
                          "range and not previously simulated (incl. x<->y transpose)",
        "script_sha256": sha256(__file__),
        "interpretation": "PROSPECTIVE GP-selected design candidates, frozen before any MD",
    }
    out_meta = CAMPAIGN / "frozen_predictions_metadata.json"
    out_meta.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf8")

    files = [out_pred, out_meta, out_man, CAMPAIGN / "candidates_all.csv", Path(__file__).resolve()]
    (CAMPAIGN / "frozen_predictions.sha256").write_text(
        "".join(f"{sha256(p)}  {p.name}\n" for p in files), encoding="utf8")
    print("\nFROZEN. Files written: geometry_manifest.csv, frozen_predictions_preMD.csv,")
    print("frozen_predictions_metadata.json, frozen_predictions.sha256")
    for c in chosen:
        print(f"  {c['label']}: N={c['N']}  kappa_GP={c['mk']:.3f}  E_GP={c['mE']:.2f}")


if __name__ == "__main__":
    main()
