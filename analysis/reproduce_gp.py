#!/usr/bin/env python3
"""Physically defined four-feature Gaussian Process validation (porous16 only).

Report LOOCV and leave-one-geometry-family-out metrics. INV1 is never fitted.
These comparisons were explored retrospectively after the corrected MD data.
"""
import csv
import warnings
from pathlib import Path
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut, LeaveOneGroupOut
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import mean_absolute_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
FEATS = ['phi_actual', 'projected_obstruction', 'orientation_chi', 'shift_mag_uc']

def load_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def model():
    kernel = ConstantKernel(1.0, (0.01,100.0)) * Matern(1.0,(0.1,100.0),nu=2.5) + WhiteKernel(0.01,(1e-4,1.0))
    return make_pipeline(StandardScaler(),GaussianProcessRegressor(
        kernel=kernel,normalize_y=True,n_restarts_optimizer=0,random_state=42))

def predict_cv(X,y,groups,mode):
    cv = LeaveOneOut() if mode == 'LOOCV' else LeaveOneGroupOut()
    pred = np.zeros(len(y))
    for tr,te in cv.split(X,y,groups if mode != 'LOOCV' else None):
        estimator = model()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',ConvergenceWarning)
            estimator.fit(X[tr],y[tr])
        pred[te] = estimator.predict(X[te])
    return pred

def main():
    descriptors = load_csv(ROOT/'data/processed/physical_descriptors_18.csv')
    actual = load_csv(ROOT/'data/processed/all_18_corrected.csv')
    desc = {r['label']:r for r in descriptors}
    table = {r['label']:r for r in actual}
    rows = [r for r in actual if r['role']=='training' and r['label']!='bulk_si']
    assert len(rows)==16 and table['INV1_AR2.5_s37']['role']=='held_out_candidate'
    X = np.array([[float(desc[r['label']][f]) for f in FEATS] for r in rows])
    inv = np.array([[float(desc['INV1_AR2.5_s37'][f]) for f in FEATS]])
    groups = np.array([desc[r['label']]['group'] for r in rows])
    print('Physical4 GP, porous training N=',len(rows),' features=',FEATS)
    print('Note: model-family/descriptor comparison was exploratory, not prospectively preregistered.')
    for target in ('kappa_z_mean_W_mK','E_x_mean_GPa'):
        y = np.array([float(r[target]) for r in rows])
        for mode in ('LOOCV','FAMILY-HOLDOUT'):
            pred = predict_cv(X,y,groups,mode)
            print(f'{target} {mode}: MAE={mean_absolute_error(y,pred):.6f}, R2={r2_score(y,pred):.6f}')
        gp=model()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',ConvergenceWarning)
            gp.fit(X,y)
        yhat=float(gp.predict(inv)[0]); measured=float(table['INV1_AR2.5_s37'][target])
        print(f'  INV1 (excluded from training): predicted={yhat:.6f} observed={measured:.6f} error_pct={100*(yhat-measured)/measured:+.2f}')

if __name__=='__main__': main()
