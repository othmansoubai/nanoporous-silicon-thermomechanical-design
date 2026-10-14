#!/usr/bin/env python3
"""Regenerate illustrative plots from verified corrected datasets."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
with (ROOT/'data/processed/all_18_corrected.csv').open(newline='') as f:rows=list(csv.DictReader(f))
with (ROOT/'data/processed/pareto_front_18.csv').open(newline='') as f:front=list(csv.DictReader(f))
fig,ax=plt.subplots(figsize=(8.0,5.4))
train=[r for r in rows if r['role']=='training']
ax.scatter([float(r['kappa_z_mean_W_mK']) for r in train],[float(r['E_x_mean_GPa']) for r in train],s=30,label='Original 17 geometries')
ax.plot([float(r['kappa_z_mean_W_mK']) for r in front],[float(r['E_x_mean_GPa']) for r in front],'-',lw=1.0,label='Measured Pareto front (18)')
inv=next(r for r in rows if r['role']=='held_out_candidate')
ax.errorbar([float(inv['kappa_z_mean_W_mK'])],[float(inv['E_x_mean_GPa'])],xerr=[float(inv['kappa_z_sample_sd_W_mK'])],yerr=[float(inv['E_x_sample_sd_GPa'])],fmt='*',markersize=13,capsize=3,label='INV1 (3 seeds)')
ax.annotate('INV1',(float(inv['kappa_z_mean_W_mK']),float(inv['E_x_mean_GPa'])),xytext=(8,-18),textcoords='offset points')
ax.set(xlabel=r'Thermal conductivity $\kappa_z$ [W m$^{-1}$ K$^{-1}$]',ylabel=r'Young modulus $E_x$ [GPa]',title='Measured thermal–mechanical trade-off')
ax.legend(loc='lower right',fontsize=8);fig.tight_layout()
fig.savefig(OUT/'measured_pareto_corrected.png',dpi=300)
plt.close(fig)
print('Saved',OUT/'measured_pareto_corrected.png')
