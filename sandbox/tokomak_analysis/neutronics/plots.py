"""Reproducible figures from saved batches/tallies; no new transport."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_results(out):
    out=Path(out)
    table=pd.read_csv(out/'module_breeding.csv');batch=np.load(out/'module_batch_production.npz')['values']
    labels=sorted(table.poloidal.unique());means=[];errors=[]
    for pol in labels:
        values=batch[:,table.poloidal.to_numpy()==pol].sum(axis=1)
        means.append(values.mean());errors.append(values.std(ddof=1)/np.sqrt(len(values)))
    fig,ax=plt.subplots(figsize=(10,5));ax.bar(labels,means,yerr=errors,capsize=3)
    ax.set(xlabel='Poloidal stack ID (first half inboard, second half outboard)',ylabel='Verified Li + Be tritons/source',title='Poloidal breeding contributions; one batch standard error')
    fig.tight_layout();fig.savefig(out/'poloidal_breeding.png',dpi=180);fig.savefig(out/'poloidal_breeding.pdf');plt.close(fig)
    totals=pd.read_csv(out/'integrated_tally_totals.csv')
    totals.to_json(out/'integrated_tally_totals.json',orient='records',indent=2)
    records=json.loads((out/'tally_manifest.json').read_text())
    wanted=next(r for r in records if r['name']=='Blanket flux spectrum')
    edges=np.array(wanted['filters'][-1]['bins']).reshape(-1,2)
    # Saved EnergyFilter bins are consecutive edge pairs, native bin integrals.
    data=np.load(out/wanted['file'])['mean'].ravel()
    fig,ax=plt.subplots(figsize=(8,5));ax.step(edges[:,0]/1e6,data,where='post')
    ax.set(xscale='log',yscale='log',xlabel='Neutron energy [MeV]',ylabel='Volume-integrated flux [cm/source / energy bin]',title='Blanket neutron spectrum (group-integrated tally)')
    fig.tight_layout();fig.savefig(out/'blanket_flux_spectrum.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);args=parser.parse_args();plot_results(args.output)
