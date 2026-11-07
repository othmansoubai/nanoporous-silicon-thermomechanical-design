#!/usr/bin/env python3
"""Read-only, source-integrity and cross-table checks for the public overlay."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'data'/'processed'
checks={
'pass7_seed_results_all12.csv':'0ccb3b6d3de9aa6f95836527357b4bfef4ff870b4176683119356e95a017b880',
'pass7_summary_all4.csv':'4ff5c96b120428b74f9f6b9e6fce5c724c9e9d63780cd3d9e00ea8ed21f94b9d',
}
for name,sha in checks.items():
    actual=hashlib.sha256((D/name).read_bytes()).hexdigest()
    assert actual==sha,(name,actual,sha)
    print('SHA256 PASS',name)
seeds=pd.read_csv(D/'pass7_seed_results_all12.csv')
summary=pd.read_csv(D/'pass7_summary_all4.csv')
full=pd.read_csv(D/'corrected_22_measured_geometries.csv')
pareto=pd.read_csv(D/'pareto_front_22_mean_values.csv')
manifest=pd.read_csv(ROOT/'data'/'jobs'/'pass7_geometry_manifest.csv')
assert len(seeds)==12 and len(summary)==4 and len(full)==22 and len(pareto)==10 and len(manifest)==4
assert seeds[['label','seed']].drop_duplicates().shape[0]==12
assert set(seeds.label)==set(summary.label)==set(manifest.label)
assert full.label.is_unique
for _,row in summary.iterrows():
    z=seeds.loc[seeds.label==row.label]
    assert len(z)==3
    for field,dst in [('kappa_z_W_mK','kappa_mean'),('E_x_GPa','E_mean')]:
        assert np.isclose(z[field].mean(),row[dst],rtol=1e-12)
    for field,dst in [('kappa_z_W_mK','kappa_sample_sd'),('E_x_GPa','E_sample_sd')]:
        assert np.isclose(z[field].std(ddof=1),row[dst],rtol=1e-12)
    original=full.set_index('label').loc[row.label]
    assert np.isclose(original.kappa,row.kappa_mean,rtol=1e-12)
    assert np.isclose(original.E,row.E_mean,rtol=1e-12)
    print('THREE-SEED PASS',row.label)
front=[]
for i,a in full.iterrows():
    dom=any(j!=i and b.kappa<=a.kappa and b.E>=a.E and (b.kappa<a.kappa or b.E>a.E) for j,b in full.iterrows())
    if not dom: front.append(a.label)
assert len(front)==10
assert set(front)==set(pareto.label)
print('PARETO PASS: 10 mean-value non-dominated structures')
print('PASS: Public overlay validated; source archives remain author-held')
