#!/usr/bin/env python3
"""Independent reduction of the released raw traces. Requires the separate raw-data package.

Thermal: 500 ps MP-NEMD, last 80% time in 20-bin temperature profiles,
C window 1-based bins 3-9 and 13-19. Mechanical: E_x from OLS slope
of sigma_xx vs engineering strain 0.0005--0.005. Prints means ± sample SD.
This verifies numerical reduction, not the physical convergence of the model.
"""
import csv
import re
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'data/raw'
PROC= ROOT/'data/processed/all_18_corrected.csv'
EV=1.602176634e-19
A=(8*5.431e-10)**2

def profiles(p):
    lines=p.read_text().splitlines()
    blocks=[]; i=0
    while i<len(lines):
        s=lines[i].strip()
        if not s or s.startswith('#'):
            i+=1;continue
        parts=s.split();step=int(parts[0]);n=int(parts[1])
        assert n==20
        vals=np.array([[float(q) for q in lines[j].split()] for j in range(i+1,i+1+n)])
        if step>=110000: blocks.append(vals)
        i+=1+n
    assert len(blocks)==40, (p,len(blocks))
    return np.stack(blocks).mean(axis=0)

def calc_kappa(tempfile,energyfile):
    prof=profiles(tempfile)
    grad=[np.polyfit(prof[sl,1]*1e-10,prof[sl,3],1)[0] for sl in (np.arange(2,9),np.arange(12,19))]
    energy=np.loadtxt(energyfile,comments='#');energy=energy[energy[:,0]>=100000]
    assert len(energy)==41
    P=np.polyfit(energy[:,0]*1e-15,energy[:,1]*EV,1)[0]
    return float(P/(2*A*np.mean(np.abs(grad))))

def calc_E(filepath):
    data=np.loadtxt(filepath,comments='#')
    strain=data[:,2]; stress=data[:,5]
    order=np.argsort(strain);strain=strain[order];stress=stress[order]
    ok=(strain>=0.0005)&(strain<=0.005)
    assert ok.sum()>=10
    return float(np.polyfit(strain[ok],stress[ok],1)[0])

def main():
    if not RAW.is_dir():raise SystemExit('Missing data/raw: extract the separate RAW_DATA zip into the same directory as the GitHub-ready package.')
    with PROC.open(newline='') as f:reported=list(csv.DictReader(f))
    errors=[]
    for r in reported:
        label=r['label'];is_inv=(r['role']=='held_out_candidate')
        thermal_seeds=(12345,45678,78901) if is_inv else (12345,34567,56789)
        mechanical_seeds=thermal_seeds
        k=[];e=[]
        for seed in thermal_seeds:
            stem=f'{label}_seed{seed}'
            k.append(calc_kappa(RAW/'thermal'/f'temp_profile_{stem}.dat',RAW/'thermal'/f'e_transfer_{stem}.dat'))
            e.append(calc_E(RAW/'mechanical'/f'stress_strain_{stem}.dat'))
        vals=(np.mean(k),np.std(k,ddof=1),np.mean(e),np.std(e,ddof=1))
        cols=('kappa_z_mean_W_mK','kappa_z_sample_sd_W_mK','E_x_mean_GPa','E_x_sample_sd_GPa')
        errors.extend(abs(v-float(r[c])) for v,c in zip(vals,cols))
        print('%-24s kappa_z = %.6f ± %.6f; E_x = %.6f ± %.6f'%(label,*vals))
    print('18 geometries, 54 thermal and 54 mechanical traces; maximum difference = %.3g'%max(errors))
    assert max(errors)<1e-8, ('Audit mismatch',max(errors))
    print('RAW-TO-PROCESSED REPRODUCIBILITY: PASS')

if __name__=='__main__':main()
