#!/usr/bin/env python3
"""Rebuild public-data manuscript figures, with one explicit slot per figure.

Run from a checkout after applying the overlay:
  python analysis/make_public_figures.py --only all
  python analysis/make_public_figures.py --only geometry

No LAMMPS, re-training, simulated raw trajectories, or internet required.
The original author-held raw MD figures are NOT claimed reproducible here.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

ROOT=Path(__file__).resolve().parents[1]
PROC=ROOT/'data'/'processed'
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True,parents=True)
BLUE='#1E658C'; ORANGE='#B65A37'; TEAL='#237A7A'; INK='#23313B'
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'pdf.fonttype':42,'savefig.facecolor':'white'})

def export(fig, filename):
    dest=OUT/filename
    fig.savefig(dest,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    print('Generated:',dest)

def geometry():
    """FIGURE SLOT 1: exact Pass7 x-y footprint projections from frozen manifest."""
    rows=pd.read_csv(ROOT/'data'/'jobs'/'pass7_geometry_manifest.csv').set_index('label')
    labels=[['ORI35_T_s2','ORI35_W_s2'],['ORI25_T_s1','ORI25_W_s1']]
    fig,axs=plt.subplots(2,2,figsize=(8.4,7.8),sharex=True,sharey=True)
    for i,pair in enumerate(labels):
        for j,label in enumerate(pair):
            ax=axs[i,j]; d=rows.loc[label]
            ax.set_facecolor('white')
            ax.add_patch(Rectangle((0,0),8,8,facecolor='#F8FAFB',edgecolor=INK,lw=1.15,zorder=0))
            A=(float(d.pAxlo),float(d.pAylo),float(d.pAxhi-d.pAxlo),float(d.pAyhi-d.pAylo))
            B=(float(d.pBxlo),float(d.pBylo),float(d.pBxhi-d.pBxlo),float(d.pByhi-d.pBylo))
            ax.add_patch(Rectangle(A[:2],A[2],A[3],facecolor=BLUE,alpha=.73,edgecolor=BLUE,lw=1.4,zorder=1))
            ax.add_patch(Rectangle(B[:2],B[2],B[3],facecolor=ORANGE,alpha=.64,edgecolor=ORANGE,lw=1.4,zorder=2))
            for txt,rec in [('A',A),('B',B)]:
                xx=rec[0]+rec[2]/2; yy=rec[1]+rec[3]/2
                ax.text(xx,yy,txt,ha='center',va='center',color='white',weight='bold',fontsize=12,zorder=4,
                        path_effects=[])   # same legend symbols used consistently below
            ax.set_aspect('equal');ax.set_xlim(-.3,8.3);ax.set_ylim(-.3,8.3)
            ax.set_xticks(np.arange(0,9,2));ax.set_yticks(np.arange(0,9,2))
            ax.tick_params(direction='out',length=3)
            ax.grid(False)
            ax.set_title(f'{label}  |  {"Tall" if j==0 else "Wide"}\n'+f'porosity {float(d.expected_porosity_pct):.2f}%   ·   N={int(d.expected_remaining_atoms):,}',pad=8,fontweight='semibold',fontsize=9.1)
            if i==1:ax.set_xlabel('x (conventional unit cells)')
            if j==0:ax.set_ylabel('y (conventional unit cells)')
    fig.subplots_adjust(left=.095,right=.98,bottom=.24,top=.92,wspace=.15,hspace=.36)
    fig.legend(handles=[Patch(facecolor=BLUE,alpha=.73,label='Pore A (first z-slice type)'),
                        Patch(facecolor=ORANGE,alpha=.64,label='Pore B (alternating z-slice type)')],
               loc='lower center',bbox_to_anchor=(.5,.154),frameon=False,ncol=2,fontsize=9)
    fig.text(.5,.102,'Mechanical tension: fixed x axis  →     |     Thermal transport: z axis, out of plane  ⊙',
             ha='center',color=INK,fontsize=9)
    fig.text(.5,.069,'Rectangles are projections from alternating slices, not two overlapping voids in one plane.',
             ha='center',color='#525B64',fontsize=8)
    export(fig,'fig_pass10_matched_footprints.pdf')

def orientation_seeds():
    """FIGURE SLOT 2: frozen GP versus all individual prospective MD seeds."""
    seeds=pd.read_csv(PROC/'pass7_seed_results_all12.csv')
    summary=pd.read_csv(PROC/'pass7_summary_all4.csv').set_index('label')
    fig,axs=plt.subplots(2,2,figsize=(9.2,6.3),layout='constrained')
    for i,(number,tag) in enumerate([('35','28.01%'),('25','19.18%')]):
        for j,(target,pred,label) in enumerate([('kappa_z_W_mK','GP_kappa',r'$\kappa_z$ (W m$^{-1}$ K$^{-1}$)'),('E_x_GPa','GP_E',r'$E_x$ (GPa)')]):
            ax=axs[i,j]
            names=[f'ORI{number}_T_s{2 if number=="35" else 1}',f'ORI{number}_W_s{2 if number=="35" else 1}']
            for k,name in enumerate(names):
                values=seeds.loc[seeds.label==name,target].astype(float).to_numpy()
                ax.scatter(k+np.array([-.065,0,.065]),values,c=BLUE if k==0 else ORANGE,s=33,label='MD seed' if k==0 else None,zorder=4)
                ax.errorbar(k,float(values.mean()),yerr=float(values.std(ddof=1)),fmt='_',markersize=19,
                            color=INK,capsize=5,lw=1.15,zorder=5)
                ax.scatter(k,float(summary.loc[name,pred]),marker='D',s=45,color='white',edgecolor=INK,zorder=6)
            ax.set(xlim=(-.4,1.4),xticks=(0,1),xticklabels=('Tall','Wide'),ylabel=label,title=f'{tag} porosity')
            ax.grid(axis='y',alpha=.17)
    fig.text(.5,-.025,'Colored circles: independent MD seeds (n=3); bars: mean ± sample SD; open diamonds: pre-MD frozen GP prediction.',ha='center',fontsize=8.6)
    export(fig,'fig_pass10_pass7_per_seed.pdf')

def pareto22():
    """FIGURE SLOT 3: Pareto points from the corrected 22-geometry sample means."""
    allv=pd.read_csv(PROC/'corrected_22_measured_geometries.csv')
    needed={'label','kappa','E'}
    assert needed.issubset(allv.columns)
    front=[]
    for i,a in allv.iterrows():
        dominated=False
        for j,b in allv.iterrows():
            if i!=j and (b.kappa<=a.kappa) and (b.E>=a.E) and ((b.kappa<a.kappa) or (b.E>a.E)):
                dominated=True;break
        if not dominated:front.append(a['label'])
    assert len(front)==10, f'Unexpected Pareto count: {len(front)}'
    fig,ax=plt.subplots(figsize=(7.1,4.5))
    reg=allv.loc[~allv.label.str.startswith('ORI')]
    prospect=allv.loc[allv.label.str.startswith('ORI')]
    ax.scatter(reg.kappa,reg.E,c='#8E9A9D',s=42,label='Original measured (18)')
    ax.scatter(prospect.kappa,prospect.E,c=ORANGE,s=78,marker='s',label='Pass 7 prospective (4)')
    f=allv.loc[allv.label.isin(front)].sort_values('kappa')
    ax.plot(f.kappa,f.E,color=TEAL,lw=1.8,ls='--',label='Mean-value Pareto envelope')
    ax.set(xlabel=r'Conductivity $\kappa_z$ (W m$^{-1}$ K$^{-1}$)',ylabel=r'Modulus $E_x$ (GPa)')
    ax.grid(alpha=.17);ax.legend(loc='lower right',frameon=False)
    ax.annotate('INV1',xy=tuple(reg[reg.label=='INV1_AR2.5_s37'][['kappa','E']].to_numpy()[0]),xytext=(10,-15),textcoords='offset points',fontsize=8)
    fig.text(.5,-.03,'Dominance evaluated on means only; this plot does not quantify uncertainty in Pareto membership.',ha='center',fontsize=8.5)
    export(fig,'fig_pareto22_public.pdf')

def nemd_sensitivity():
    """FIGURE SLOT 4: archived 51-run gradient/temporal-fit sensitivity."""
    source=PROC/'thermal_runs_51_with_QC.csv'
    if not source.exists():
        raise FileNotFoundError(f'{source} missing. Keep existing public repo thermal_runs_51_with_QC.csv.')
    df=pd.read_csv(source)
    assert len(df)==51
    fig,axes=plt.subplots(1,2,figsize=(9.7,3.4),layout='constrained')
    for ax,field,lab in zip(axes,['C_D_pct','temporal_pct'],['Interior window D versus C (%)','Final 50% versus 80% (%)']):
        v=df[field].astype(float)
        ax.scatter(df.phi_pct,v,c=df.kappa_C,cmap='viridis',s=30,edgecolors='white',linewidths=.25)
        ax.axhline(0,color='#777777',ls='--',lw=.8)
        ax.set(xlabel='Atom-counted porosity (%)',ylabel=lab,title=f'51 runs; median |difference| = {np.median(np.abs(v)):.2f}%')
        ax.grid(alpha=.13)
    fig.text(.5,-.04,'Sensitivity of reduction choices on the same MD trajectories; not an independent uncertainty estimate.',ha='center',fontsize=8.2)
    export(fig,'fig_pass10_nemd_window_sensitivity.pdf')

def main():
    a=argparse.ArgumentParser();a.add_argument('--only',choices=['all','geometry','orientation','pareto','sensitivity'],default='all')
    q=a.parse_args().only
    actions={'geometry':geometry,'orientation':orientation_seeds,'pareto':pareto22,'sensitivity':nemd_sensitivity}
    for k,fn in actions.items():
        if q in ('all',k):fn()

if __name__=='__main__':main()
