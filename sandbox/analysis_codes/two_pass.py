import argparse
from dataclasses import replace
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
from .analysis import run_analysis, _label, _savefig
from sandbox.analysis_codes.storage import save_json
from .settings import FULL

def screen_plots(output):
    """Plot every screened input and retain uncertainty in machine-readable data."""
    import matplotlib.pyplot as plt
    table=pd.read_csv(output/'morris_results.csv')
    ranked=table[table.response=='total'].sort_values('mu_star',ascending=False)
    for page,start in enumerate(range(0,len(ranked),30),1):
        block=ranked.iloc[start:start+30].iloc[::-1]
        fig,ax=plt.subplots(figsize=(12,max(5,len(block)*.27)))
        ax.barh(np.arange(len(block)),block.mu_star,xerr=np.array([
            np.maximum(0,block.mu_star-block.mu_star_ci_low),
            np.maximum(0,block.mu_star_ci_high-block.mu_star)]))
        ax.set_yticks(np.arange(len(block)),[_label(p) for p in block.parameter])
        ax.set(xlabel='Morris μ* [tritons/source neutron per unit design coordinate]',
               title=f'Low-statistics screening — all inputs, page {page}')
        _savefig(output,'morris_ranking' if page==1 else f'morris_ranking_{page}')
    plt.figure(figsize=(7,5))
    for family,part in ranked.groupby('family'):
        plt.scatter(part.mu_star,part.sigma,label=family,alpha=.7)
    plt.xlabel('Morris μ*');plt.ylabel('Morris σ');plt.legend()
    _savefig(output,'morris_scatter')
    groups=ranked.groupby('family').agg(observed=('mu_star','sum'),noise=('noise_rms','sum'))
    groups.noise*=np.sqrt(2/np.pi)
    groups.rename(columns={'observed':'Observed sum of μ*','noise':'Approximate zero-effect noise reference'}).plot.bar(figsize=(9,5))
    plt.ylabel('Sum of μ*; parameter-count dependent, not Sobol indices')
    _savefig(output,'morris_families')


def write_sensitivity_index(output):
    """Expose all saved Morris/Sobol figures in one persistent gallery."""
    output=Path(output)
    lines=['# Morris and Sobol sensitivity plots\n',
        'Morris rankings can remain noise limited; high-pass Sobol indices are conditional on retained inputs.\n']
    for stage,title in [('high','Higher-statistics refinement'),('low','Initial lower-statistics screening')]:
        lines.append(f'## {title}\n')
        figures=sorted((output/stage).glob('morris*.png'))
        figures+=sorted((output/stage).glob('sobol*.png'))
        for figure in figures:
            relative=figure.relative_to(output).as_posix()
            lines.append(f'### {figure.stem}\n\n![{figure.stem}]({relative})\n\n[PDF]({relative[:-4]}.pdf)\n')
    (output/'sensitivity_plots.md').write_text('\n'.join(lines))
