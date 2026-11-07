#!/usr/bin/env python3
"""Freeze current physical-four GP predictions BEFORE viewing new MD outputs.
Use on a copy of the verified GitHub repository at the indicated commit.
Only the original 16 porous training structures are fitted; INV1 is excluded.
Running this script does NOT select INV2 and does NOT submit any jobs.
"""
import csv
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

CAMPAIGN = Path(__file__).resolve().parent
REPO = Path(sys.argv[1]).expanduser().resolve() if len(sys.argv)>1 else Path.home()/'nanoporous-silicon-thermomechanical-design'
REQUIRED_COMMIT='0de9bd5559b196a4f62f3afde9cbf94986bd17eb'

head=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip()
if head!=REQUIRED_COMMIT:
    raise SystemExit(f"STOP: repository commit is {head}; planned freeze is {REQUIRED_COMMIT}. Inspect before proceeding.")

sys.path.insert(0, str(REPO/'analysis'))
from reproduce_gp import FEATS, model, load_csv
import numpy as np
import sklearn

assert FEATS==['phi_actual','projected_obstruction','orientation_chi','shift_mag_uc']
actual=load_csv(REPO/'data/processed/all_18_corrected.csv')
desc={r['label']:r for r in load_csv(REPO/'data/processed/physical_descriptors_18.csv')}
train=[r for r in actual if r['role']=='training' and r['label']!='bulk_si']
assert len(train)==16
X=np.array([[float(desc[r['label']][f]) for f in FEATS] for r in train],float)

manifest=list(csv.DictReader((CAMPAIGN/'geometry_manifest.csv').open(newline='',encoding='utf8')))
X_new=[]
for r in manifest:
    dx=float(r['pAxhi'])-float(r['pAxlo']);dy=float(r['pAyhi'])-float(r['pAylo'])
    dx_B=float(r['pBxhi'])-float(r['pBxlo']);dy_B=float(r['pByhi'])-float(r['pBylo'])
    assert abs(dx-dx_B)<1e-12 and abs(dy-dy_B)<1e-12
    shift_x=float(r['pBxlo'])-float(r['pAxlo'])
    shift_y=float(r['pBylo'])-float(r['pAylo'])
    overlap_x=max(0,min(float(r['pAxhi']),float(r['pBxhi']))-max(float(r['pAxlo']),float(r['pBxlo'])))
    overlap_y=max(0,min(float(r['pAyhi']),float(r['pByhi']))-max(float(r['pAylo']),float(r['pBylo'])))
    obstruction=(2*dx*dy-overlap_x*overlap_y)/64
    phi=round(float(r['expected_porosity_pct']),2)  # match release feature convention
    chi=(dx-dy)/(dx+dy)
    X_new.append([phi,obstruction,chi,math.hypot(shift_x,shift_y)])
X_new=np.array(X_new,float)

out=CAMPAIGN/'frozen_predictions_preMD.csv'
if out.exists():
    raise SystemExit(f'File already exists; refusing to replace the pre-MD record: {out}')

pred={}
for key in ['kappa_z_mean_W_mK','E_x_mean_GPa']:
    y=np.array([float(r[key]) for r in train],float)
    m=model()
    m.fit(X,y)
    mu,sd=m.predict(X_new,return_std=True)
    pred[key]=(mu,sd)

with out.open('w',newline='',encoding='utf8') as f:
    w=csv.writer(f)
    w.writerow(['label',*FEATS,'GP_pred_kappa_z_W_mK','GP_model_sd_kappa_z_W_mK','GP_pred_E_x_GPa','GP_model_sd_E_x_GPa'])
    for j,r in enumerate(manifest):
        w.writerow([r['label'],*X_new[j],pred['kappa_z_mean_W_mK'][0][j],pred['kappa_z_mean_W_mK'][1][j],pred['E_x_mean_GPa'][0][j],pred['E_x_mean_GPa'][1][j]])

metadata={'utc_generated':datetime.now(timezone.utc).isoformat(),
          'repository_commit':head, 'sklearn_version':sklearn.__version__,
          'training_n_porous':16,'features':FEATS,
          'manifest_sha256':hashlib.sha256((CAMPAIGN/'geometry_manifest.csv').read_bytes()).hexdigest(),
          'predictions_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
          'interpretation':'PROSPECTIVE predictions for orientation controls before MD; not the historical INV1 design selection'}
(CAMPAIGN/'frozen_predictions_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf8')
(CAMPAIGN/'frozen_predictions.sha256').write_text(
    hashlib.sha256(out.read_bytes()).hexdigest()+'  frozen_predictions_preMD.csv\n'+
    hashlib.sha256((CAMPAIGN/'frozen_predictions_metadata.json').read_bytes()).hexdigest()+'  frozen_predictions_metadata.json\n'+
    hashlib.sha256((CAMPAIGN/'geometry_manifest.csv').read_bytes()).hexdigest()+'  geometry_manifest.csv\n',encoding='utf8')
print('FROZEN: ',out)
for r in csv.DictReader(out.open(newline='',encoding='utf8')):
    print(r['label'], 'kappa=',r['GP_pred_kappa_z_W_mK'],'E=',r['GP_pred_E_x_GPa'])
