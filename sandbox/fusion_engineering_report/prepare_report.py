"""Rebuild manuscript tables/figures from saved data; never launch transport.

Run from anywhere with the repository's Python environment. Outputs use only
relative paths so the resulting LaTeX project is portable to Overleaf.
"""
from pathlib import Path
import hashlib,json,shutil,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent.parent
STUDIES=[('Li_analysis','Li sphere'),('FLiBe_analysis','FLiBe sphere'),('LiPb_analysis','PbLi sphere'),('Be_Li_layers','Layers, parallel'),('Be_Li_layers_x_beam','Layers, perpendicular')]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
for d in ('figures','tables','data'): (ROOT/d).mkdir(exist_ok=True)
manifest=[];summary=[];sobols={};morris={};family={};bundles={}
def record(source,target):
 source=Path(source);target=Path(target);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
 manifest.append({'source':str(source.relative_to(REPO)),'copy':str(target.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
def escape(x):
 return str(x).replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('%',r'\%').replace('&',r'\&').replace('#',r'\#')
def label(row):
 mode=''
 if row.family!='xs': mode=' / '+('P' if row.family=='angle' else 'E')+row.parameter.split(':m')[-1].split(':')[0]
 return f'{row.nuclide} MT{row.mt} {row.family}'+mode
for directory,name in STUDIES:
 result=REPO/'sandbox'/'studies'/directory/'results';tag=directory.lower();high=result/'high';low=result/'low'
 s=pd.read_csv(high/'sobol_indices.csv');s['label']=[label(row) for row in s.itertuples()];sobols[name]=s
 m=pd.read_csv(low/'morris_results.csv');m=m[m.response=='total'].copy();m['label']=[label(row) for row in m.itertuples()];morris[name]=m
 f=pd.read_csv(high/'sobol_families.csv');family[name]=f
 v=json.loads((high/'pce_validation.json').read_text());vc=v['validation'];b=json.loads((high/'baseline.json').read_text());dist=json.loads((high/'tritium_distribution.json').read_text());lo=json.loads((low/'execution.json').read_text());hi=json.loads((high/'execution.json').read_text());threshold=json.loads((low/'threshold_filter.json').read_text());stop=json.loads((low/'morris_stopping.json').read_text());conf=json.loads((high/'configuration.json').read_text());bundle=np.load(high/'analysis_arrays.npz');bundles[name]={k:bundle[k] for k in bundle.files};pce=np.load(high/'pce_model.npz')
 baseline=b['total'];entry=dict(directory=directory,name=name,screened=len(m),retained=len(s),morris_histories=lo['total_histories'],morris_evaluations=lo['transport_evaluations'],morris_converged=stop['converged'],threshold_removed=threshold['parameters_removed'],reactions_removed=threshold['reactions_removed'],high_histories=hi['histories'],high_particles=conf['particles'],high_batches=conf['batches'],pce_order=v['order'],pce_terms=len(pce['powers']),max_interaction=int((pce['powers']>0).sum(axis=1).max()),r2=vc['r2'][0],rmse=vc['rmse'][0],rmse_over_validation_sd=float(vc['rmse'][0]/bundle['val_y'][:,0].std(ddof=1)),noise_mse=vc['noise_mse'][0],accepted=vc['accepted'][0],baseline=baseline['mean'],baseline_se=baseline['standard_error'],**dist)
 summary.append(entry)
 destination=ROOT/'data'/directory
 for source in [high/'sobol_indices.csv',high/'sobol_families.csv',high/'pce_validation.json',high/'tritium_distribution.json',high/'baseline.json',high/'configuration.json',high/'provenance.json',high/'analysis_arrays.npz',high/'pce_model.npz',low/'morris_results.csv',low/'morris_convergence.json',low/'morris_stopping.json',low/'threshold_filter.json',low/'model_details.json',low/'study_config.json',low/'execution.json',high/'execution.json',result/'retained_parameters.csv',result/'removed_reactions.csv',result/'selected_parameters.csv',result/'verification/controlled_tests.csv',result/'provenance/endf_channel_metadata.csv',result/'provenance/production_mapping.json',REPO/'sandbox'/'studies'/directory/'config.py',REPO/'sandbox'/'studies'/directory/'model.py']:
  if source.exists():record(source,destination/(source.parent.name+'_'+source.name if source.name in ['execution.json','provenance.json'] else source.name))
 for filename in ['pce_validation.png','sobol_indices.png','tritium_distribution.png']:
  record(high/filename,ROOT/'figures'/f'{tag}_{filename}')
 # Complete Sobol tables, including unresolved selected parameters.
 lines=[r'\begin{longtable}{p{0.59\linewidth}rrl}',r'\caption{All retained inputs: '+escape(name)+r'. Indices are conditional; unresolved flags follow saved coefficient-noise diagnostics.}\\',r'\toprule Parameter & $S_i$ & $S_{T_i}$ & Resolved\\\midrule\endfirsthead',r'\toprule Parameter & $S_i$ & $S_{T_i}$ & Resolved\\\midrule\endhead']
 for row in s.sort_values('total',ascending=False).itertuples():lines.append(f'{escape(row.label)} & {row.first:.5f} & {row.total:.5f} & '+('yes' if row.resolved else 'no')+r'\\')
 lines +=[r'\bottomrule\end{longtable}'];(ROOT/'tables'/f'{tag}_sobol.tex').write_text('\n'.join(lines))
 # Six-channel comparison tables used in the main text.
 top=[r'\begin{tabular}{llrr}',r'\toprule Nuclide & MT / family / mode & $S_i$ & $S_{T_i}$\\\midrule']
 for row in s.sort_values('total',ascending=False).head(6).itertuples():
  mode='' if row.family=='xs' else (' / '+row.parameter.split(':m')[-1].split(':')[0])
  top.append(f'{row.nuclide} & {row.mt} / {row.family}{mode} & {row.first:.4f} & {row.total:.4f}'+r'\\')
 top.append(r'\bottomrule\end{tabular}')
 (ROOT/'tables'/f'{tag}_top.tex').write_text('\n'.join(top))
 # Complete Morris index is a portable searchable CSV rather than 100+ printed pages.
 m.to_csv(destination/'total_morris_report.csv',index=False)
pd.DataFrame(summary).to_csv(ROOT/'data/study_summary.csv',index=False)
(ROOT/'data/study_summary.json').write_text(json.dumps(summary,indent=2))
# Authoritative native geometry images and physics results.
for source,new in [('sandbox/tokomak_analysis/results/detailed/visualization/xz_full.png','tokamak_xz.png'),('sandbox/tokomak_analysis/results/detailed/visualization/cutaway.png','tokamak_cutaway.png'),('sandbox/tokomak_analysis/results/detailed/visualization/blanket_outboard_detail.png','tokamak_blanket.png'),('sandbox/tokomak_analysis/results/detailed/visualization/divertor_detail.png','tokamak_divertor.png'),('sandbox/tokomak_analysis/results/segmented/poloidal_breeding.png','tokamak_breeding.png'),('sandbox/tokomak_analysis/results/segmented/blanket_flux_spectrum.png','tokamak_spectrum.png'),('sandbox/studies/Be_Li_layers/results/geometry_xz.png','layers_geometry.png')]:record(REPO/source,ROOT/'figures'/new)
for sub in ('global','segmented','detailed','coupled_smoke'):
 for filename in ('results.json','build_metrics.json','config.json'):
  record(REPO/'sandbox/tokomak_analysis/results'/sub/filename,ROOT/'data/tokamak'/sub/filename)
for filename in ('implementation_report.md','references.md'):record(REPO/'sandbox/tokomak_analysis'/filename,ROOT/'data/tokamak'/filename)
# Native readback perturbation plots, with their original metadata (no re-sampling).
for pattern,new in [('Li_analysis/results/nuclear_data_plots/xs_rank1_Li7_mt56_p0_d0_mode0.png','li7_mt56_xs.png'),('Li_analysis/results/nuclear_data_plots/xs_rank2_Li6_mt105_p0_d0_mode0.png','li6_mt105_xs.png'),('Li_analysis/results/nuclear_data_plots/angle_rank1_Li7_mt2_p0_d0_mode2.png','li7_mt2_angle.png'),('Li_analysis/results/nuclear_data_plots/energy_rank1_Li7_mt16_p0_d0_mode1.png','li7_mt16_energy.png')]:
 source=REPO/'sandbox'/'studies'/pattern
 if not source.exists():raise FileNotFoundError(source)
 record(source,ROOT/'figures'/new)
 for metadata in source.parent.glob(source.stem+'*json'):record(metadata,ROOT/'data/native_plots'/metadata.name)
# Comparative figures, explicitly reduced-space quantities.
fig,axes=plt.subplots(1,2,figsize=(12,4.8));names=[s['name'] for s in summary];x=np.arange(len(names))
axes[0].errorbar(x,[s['baseline'] for s in summary],yerr=[s['baseline_se'] for s in summary],fmt='o',capsize=4,label='Native baseline ± 1 SE');axes[0].set(ylabel='Verified tritons per source neutron',title='Responses differ by geometry and material');axes[0].legend(fontsize=8)
axes[1].bar(x,[100*s['standard_deviation']/s['mean'] for s in summary],color='#2b6a8c');axes[1].set(ylabel='Conditional coefficient of variation [%]',title='Independent hypothetical input priors')
for ax in axes:ax.set_xticks(x,names,rotation=30,ha='right')
fig.tight_layout();fig.savefig(ROOT/'figures/response_comparison.pdf');fig.savefig(ROOT/'figures/response_comparison.png');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,5));colors={'xs':'#366b95','energy':'#d79137','angle':'#4b9872'}
for i,name in enumerate(names):
 f=family[name].set_index('family');bottom=0
 for fam in ('xs','angle','energy'):
  value=float(f.loc[fam,'first']);axes[0].bar(i,value,bottom=bottom,color=colors[fam],label=fam if i==0 else None);bottom+=value
axes[0].set_xticks(x,names,rotation=30,ha='right');axes[0].set(ylabel='Closed family Sobol index',title='Conditional variance attribution');axes[0].legend()
keys=[('Li7',56,'xs'),('Li6',105,'xs'),('Be9',2,'xs'),('Li7',2,'xs')];width=.18
for j,(nuclide,mt,fam) in enumerate(keys):
 vals=[]
 for name in names[-2:]:
  a=sobols[name];a=a[(a.nuclide==nuclide)&(a.mt==mt)&(a.family==fam)];vals.append(float(a.total.iloc[0]) if len(a) else np.nan)
 axes[1].bar(np.arange(2)+(j-1.5)*width,vals,width,label=f'{nuclide} MT{mt}')
axes[1].set_xticks([0,1],['Parallel to layers','Across layers']);axes[1].set(ylabel='Conditional total-order Sobol index',title='Same inventory: source/geometry coupling matters');axes[1].legend(fontsize=8)
fig.tight_layout();fig.savefig(ROOT/'figures/geometry_and_families.pdf');fig.savefig(ROOT/'figures/geometry_and_families.png');plt.close(fig)
fig,axes=plt.subplots(3,2,figsize=(12,13));
for ax,(directory,name) in zip(axes.flat,STUDIES):
 a=morris[name].sort_values('mu_star',ascending=False).head(10).iloc[::-1];lo=np.maximum(0,a.mu_star-a.mu_star_ci_low);hi=np.maximum(0,a.mu_star_ci_high-a.mu_star)
 ax.barh(a.label,a.mu_star,xerr=[lo,hi],color=[colors[f] for f in a.family],capsize=2);ax.set(title=name,xlabel=r'$\mu^*$ [tritons/source per unit Morris coordinate]');ax.tick_params(axis='y',labelsize=7)
axes.flat[-1].axis('off');axes.flat[-1].text(0,.9,'Top ten shown for comparison only.\nEvery analysed parameter is retained\nin data/*/total_morris_report.csv.\n\nWhiskers: trajectory-bootstrap 95% intervals.\nScreening unconverged in three studies.',va='top')
fig.tight_layout();fig.savefig(ROOT/'figures/morris_comparison.pdf');fig.savefig(ROOT/'figures/morris_comparison.png');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4.7))
for directory,name in STUDIES:
 d=json.loads((REPO/directory/'results/low/morris_convergence.json').read_text());axes[0].plot([v['trajectories'] for v in d],[v['rho'] for v in d],'-o',ms=3,label=name);axes[1].plot([v['trajectories'] for v in d],[v['top10_overlap'] for v in d],'-o',ms=3,label=name)
axes[0].axhline(.9,color='gray',ls='--');axes[1].axhline(.8,color='gray',ls='--')
axes[0].set(xlabel='Trajectories',ylabel='Split-half Spearman correlation',ylim=(0,1.02));axes[1].set(xlabel='Trajectories',ylabel='Top-ten overlap',ylim=(0,1.02));axes[0].legend(fontsize=8);fig.tight_layout();fig.savefig(ROOT/'figures/morris_convergence.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4.6))
for ax,name in zip(axes,['Layers, parallel','Layers, perpendicular']):
 b=bundles[name];v=b['validation'];y=b['val_y'][:,0];se=b['val_se'][:,0];directory=dict((name,directory) for directory,name in STUDIES)[name];p=np.load(REPO/directory/'results/high/pce_model.npz');
 from numpy.polynomial.legendre import legval
 design=np.ones((len(v),len(p['powers'])))
 for j in range(v.shape[1]):
  for degree in np.unique(p['powers'][:,j]):design[:,p['powers'][:,j]==degree]*=(np.sqrt(2*degree+1)*legval(2*v[:,j]-1,[0]*degree+[1]))[:,None]
 pred=(design@p['coefficients'])[:,0];ax.errorbar(y,pred,xerr=se,fmt='.',alpha=.65,ms=4);limits=[min(y.min(),pred.min()),max(y.max(),pred.max())];ax.plot(limits,limits,'k--',lw=1);ax.set(title=name,xlabel='Native OpenMC [tritons/source]',ylabel='PCE [tritons/source]')
fig.tight_layout();fig.savefig(ROOT/'figures/layer_validation.pdf');plt.close(fig)
# Complete selected-input comparison: NaN means not retained, not zero.
merged={}
for name,a in sobols.items():
 for row in a.itertuples():merged.setdefault(row.parameter,{'parameter':row.parameter,'label':row.label})[name]=row.total
pd.DataFrame(merged.values()).to_csv(ROOT/'data/cross_model_sobol.csv',index=False)
# Auto tables reduce transcription errors.
lines=[r'\begin{tabular}{lrrrrrr}',r'\toprule Study & $d$ & $d_a$ & $p$ & $M$ & $R^2$ & RMSE\\\midrule']
for s in summary:lines.append(f"{escape(s['name'])} & {s['screened']} & {s['retained']} & {s['pce_order']} & {s['pce_terms']} & {s['r2']:.5f} & {s['rmse']:.3g}"+r'\\')
lines.append(r'\bottomrule\end{tabular}');(ROOT/'tables/validation.tex').write_text('\n'.join(lines))
lines=[r'\begin{tabular}{lrrrr}',r'\toprule Study & Baseline $\pm$ SE & PCE mean & PCE SD & 5--95\% interval\\\midrule']
for s in summary:lines.append(f"{escape(s['name'])} & ${s['baseline']:.6f}\\pm{s['baseline_se']:.6f}$ & {s['mean']:.6f} & {s['standard_deviation']:.6f} & [{s['q05']:.6f}, {s['q95']:.6f}]"+r'\\')
lines.append(r'\bottomrule\end{tabular}');(ROOT/'tables/distributions.tex').write_text('\n'.join(lines))
(ROOT/'data/source_manifest.json').write_text(json.dumps(manifest,indent=2));print('Prepared',len(manifest),'source artifacts; no transport executed.')
