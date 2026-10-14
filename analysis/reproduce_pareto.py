#!/usr/bin/env python3
"""Reproduce expanded measured Pareto front: minimize kappa_z, maximize E_x."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
with (ROOT/'data/processed/all_18_corrected.csv').open(newline='') as f:
    rows=list(csv.DictReader(f))
assert len(rows)==18 and sum(r['role']=='held_out_candidate' for r in rows)==1

def dominates(a,b):
    ak,ae=float(a['kappa_z_mean_W_mK']),float(a['E_x_mean_GPa'])
    bk,be=float(b['kappa_z_mean_W_mK']),float(b['E_x_mean_GPa'])
    return ak<=bk and ae>=be and (ak<bk or ae>be)

front=[r for r in rows if not any(dominates(s,r) for s in rows if s is not r)]
front.sort(key=lambda r:float(r['kappa_z_mean_W_mK']))
assert len(front)==9
print('Measured 18-geometry Pareto front (means only):')
for r in front:
    print('%-24s %9.4f W/(m K) %9.3f GPa'%(r['label'],float(r['kappa_z_mean_W_mK']),float(r['E_x_mean_GPa'])))
inv=next(r for r in rows if r['role']=='held_out_candidate')
dominated=sorted(r['label'] for r in rows if r is not inv and dominates(inv,r))
print('Designs dominated by INV1:',', '.join(dominated))
assert len(dominated)==6
