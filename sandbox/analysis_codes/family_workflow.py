import argparse

from dataclasses import asdict,replace

from pathlib import Path

import json

import shutil

import time

import numpy as np

import pandas as pd

import matplotlib

import matplotlib.pyplot as plt

from scipy.stats import qmc

import openmc

from sandbox.analysis_codes import pce,sobol

from sandbox.analysis_codes.validation import validate

from sandbox.analysis_codes.registry import Parameter

from sandbox.analysis_codes.storage import save_json

from .analysis import run_analysis,_savefig,_label,_markdown

from .budget import BudgetedRunner

from .settings import FULL
from .context import build_model, responses as baseline_responses, register
FAMILIES=("xs","energy","angle")

def select_per_family(results, registry, top):
    """Keep up to n largest Morris mu-star inputs IN EACH family, even unresolved.

    Monte Carlo uncertainty is carried into the selection table. A low-history
    positive mu-star is not by itself evidence of a statistically resolved effect.
    """
    if top<1:raise ValueError('top must be positive')
    ranked=results[results.response=='total'].copy()
    known=set(registry.parameter_name)
    if set(ranked.parameter)!=known:raise ValueError('Morris inputs and retained registry differ')
    selected=[]
    for family in FAMILIES:
        block=ranked[ranked.family==family].sort_values(['mu_star','parameter'],ascending=[False,True]).head(top).copy()
        block['family_rank']=np.arange(1,len(block)+1)
        block['selection_reason']='Top n within family by Morris mu-star; keep all in high-history PCE; unresolved status retained'
        selected.append(block)
    result=pd.concat(selected,ignore_index=True)
    assert not result.parameter.duplicated().any()
    return result


def high_transport(output, selected, registry, training=512,validation=96, *,
                   particles=10000,batches=10,seed=FULL.seed,history_cap=73000000,
                   model_factory=build_model,response_factory=baseline_responses,spec_factory=register,
                   response_label=None,workers=1):
    """Use the existing Runner for all selected inputs; no second down-selection."""
    output=Path(output);hi=output/'high';hi.mkdir(exist_ok=True)
    by_id={r['parameter_name']:r for r in registry}
    chosen=selected.parameter.tolist()
    params=[]
    for identity in chosen:
        data={key:by_id[identity][key] for key in Parameter.__dataclass_fields__}
        data['interval']=tuple(data['interval']);params.append(Parameter(**data))
    limits=np.array([p.limit for p in params])
    setup=dict(parameter_ids=chosen,training=training,validation=validation,seed=seed,
               particles=particles,batches=batches,response=response_label or 'Li6 MT105 + Li7 MT52..82',
               method='No second Morris or parameter re-selection')
    previous=hi/'configuration.json'
    if previous.exists() and json.loads(previous.read_text())!=setup:
        raise ValueError('Preserving high-history observations: changed settings need a new output')
    save_json(previous,setup)
    openmc.reset_auto_ids();model=model_factory(particles,batches)
    spec=spec_factory(model,response_factory(model));names=list(spec.responses)
    histories=particles*batches
    start=time.perf_counter()
    with BudgetedRunner(spec,hi/'checkpoints',history_cap) as runner:
        runner.bind(params)
        save_json(hi/'provenance.json',runner.provenance)
        def evaluate(unit,stage,seed):
            record=runner.evaluate((2*np.asarray(unit)-1)*limits,stage=stage,seed=seed)
            return (np.array([record['responses'][name][0][0] for name in names]),
                    np.array([record['responses'][name][1][0] for name in names]))
        baseline,base_se=evaluate(np.full(len(params),.5),'family_high_baseline',seed)
        save_json(hi/'baseline.json',{n:dict(mean=baseline[i],standard_error=base_se[i]) for i,n in enumerate(names)})
        def evaluate_design(design,stage,offset):
            if workers==1:
                for j,point in enumerate(design):
                    yield evaluate(point,'family_high_'+stage,seed+offset+j)
            else:
                from concurrent.futures import ProcessPoolExecutor
                import multiprocessing as mp
                from .analysis import _initialize_worker,_worker_evaluate
                metadata=dict(name=spec.name,responses=spec.responses,
                    target_nuclides=spec.target_nuclides,threads=spec.threads,
                    history_budget=history_cap,paired_morris=False,
                    inputs=str((hi/'checkpoints/model_inputs'/runner.provenance['model_sha256']).resolve()))
                tasks=((point,'family_high_'+stage,seed+offset+j) for j,point in enumerate(design))
                with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),
                    initializer=_initialize_worker,initargs=(metadata,hi/'checkpoints',params,limits)) as pool:
                    yield from pool.map(_worker_evaluate,tasks)
        designs={};responses={};errors={}
        for stage,count,offset in [('training',training,300000),('validation',validation,400000)]:
            x=qmc.LatinHypercube(len(params),seed=seed+(31 if stage=='training' else 37)).random(count)
            designs[stage]=x;means=[];se=[];records=[]
            for j,(point,result) in enumerate(zip(x,evaluate_design(x,stage,offset))):
                point_seed=seed+offset+j
                m,s=result
                means.append(m);se.append(s)
                records.append(dict(sample=j,seed=point_seed,histories=histories,
                    **{p.identity:float((2*u-1)*p.limit) for p,u in zip(params,point)},
                    **{n:float(m[k]) for k,n in enumerate(names)},
                    **{n+'_se':float(s[k]) for k,n in enumerate(names)}))
                if j%32==0:print(f'HIGH {stage} {j+1}/{count}',flush=True)
            pd.DataFrame(records).to_csv(hi/('pce_'+stage+'.csv'),index=False)
            responses[stage]=np.asarray(means);errors[stage]=np.asarray(se)
        np.savez_compressed(hi/'analysis_arrays.npz',training=designs['training'],validation=designs['validation'],
            train_y=responses['training'],train_se=errors['training'],val_y=responses['validation'],val_se=errors['validation'],
            names=np.asarray(names),baseline=baseline,baseline_se=base_se)
        save_json(hi/'execution.json',dict(evaluations=1+training+validation,histories=(1+training+validation)*histories,
            reserved_histories=runner.history_budget.used(),elapsed_seconds=time.perf_counter()-start,
            exact_native_restoration=True,workers=workers,threads=getattr(spec,"threads",1)))
    return params


def analyse_saved(output,fit_options=None,response_label="Lithium tritium production [tritons/source neutron]",analysis_seed=FULL.seed, backend="numpy", orders=(2,3), folds=4, distribution_samples=100000, prior_limit=.2, validation_options=None):
    """Fit and plot using saved native transport observations only."""
    output=Path(output);hi=output/'high';bundle=np.load(hi/'analysis_arrays.npz')
    selected=pd.read_csv(output/'selected_parameters.csv');d=len(selected)
    registry_file=output/'retained_parameters.json'
    if registry_file.exists():
        by_id={row['parameter_name']:row for row in json.loads(registry_file.read_text())}
        if any(not np.isclose(by_id[identity]['limit'],prior_limit,rtol=0,atol=1e-14) for identity in selected.parameter):
            raise ValueError('Configured perturbation limit differs from saved transport data; use a new output for a new range')
    x=bundle['training'];y=bundle['train_y'];se=bundle['train_se']
    # Keep all main effects: relative-to-largest sparse pruning can erase a
    # weak family simply because one XS dominates. Ridge and training-only CV
    # control the fit; component resolution is reported separately.
    options=dict(q=.75,interaction=2,ridge=1e-6,sparse=False)
    if fit_options:options.update(fit_options)
    basis_limit=options.pop('basis_limit',800)
    orders=[order for order in orders if len(pce.indices(d,order,q=options['q'],interaction=options['interaction']))<=basis_limit]
    if not orders:raise ValueError("No candidate PCE order fits PCE_BASIS_LIMIT")
    engine=pce
    if backend=="openturns":
        from . import openturns_pce as engine
    elif backend!="numpy":raise ValueError("PCE_BACKEND must be numpy or openturns")
    order,cv=engine.cross_validate(x,y,se,orders=orders,folds=folds,seed=analysis_seed,
        **options)
    model=engine.fit(x,y,se,order=order,**options)
    model.save(hi/'pce_model.npz')
    checks=validate(model,bundle['validation'],bundle['val_y'],bundle['val_se'],**(validation_options or {}))
    save_json(hi/'pce_validation.json',dict(order=order,training_only_cv=cv,validation=checks,
        diagnostics=model.diagnostics,orders_considered=orders,fit_options=options,analysis_seed=analysis_seed,validation_options=validation_options or {}))
    names=bundle['names'].tolist();total=names.index('total')
    groups={family:np.flatnonzero(selected.family.to_numpy()==family).tolist() for family in FAMILIES}
    sensitivity=sobol.analyse(model,groups=groups,accepted=checks['accepted'])
    result=selected[['parameter','family','nuclide','mt','family_rank']].copy()
    result['first']=sensitivity['first'][:,total];result['total']=sensitivity['total'][:,total]
    result['resolved']=sensitivity['total_component_resolved'][:,total]
    result.to_csv(hi/'sobol_indices.csv',index=False)
    families=pd.DataFrame([dict(family=f,first=sensitivity['groups'][f]['first'][total],
        total=sensitivity['groups'][f]['total'][total],resolved=bool(sensitivity['groups'][f]['resolved'][total])) for f in groups])
    families.to_csv(hi/'sobol_families.csv',index=False)
    fig,ax=plt.subplots(figsize=(11,max(6,d*.3)));positions=np.arange(d)
    ax.barh(positions-.18,result['first'],height=.36,label='First order')
    ax.barh(positions+.18,result['total'],height=.36,label='Total order')
    ax.set_yticks(positions,[_label(p)+(' *' if not r else '') for p,r in zip(result.parameter,result.resolved)])
    ax.invert_yaxis();ax.set(xlabel='Conditional Sobol index',title='Every retained input — * unresolved component; no inputs hidden')
    ax.legend();_savefig(hi,'sobol_indices')
    # Zoom each family without renormalising its contribution to global variance.
    # Small unresolved estimates must remain visible and explicitly labelled.
    for family in FAMILIES:
        block=result[result.family==family]
        fig,ax=plt.subplots(figsize=(10,max(4,len(block)*.35)))
        position=np.arange(len(block))
        ax.barh(position-.18,block['first'],height=.36,label='First order')
        ax.barh(position+.18,block['total'],height=.36,label='Total order')
        ax.set_yticks(position,[_label(p)+(' *' if not r else '') for p,r in zip(block.parameter,block.resolved)])
        ax.invert_yaxis();ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,3))
        ax.set(xlabel='Sobol index / total selected-input variance (not renormalised)',
            title=f'{family}: every selected input; * unresolved; magnified family axis')
        ax.legend();_savefig(hi,'sobol_'+family)
    families.set_index('family')[['first','total']].plot.bar(figsize=(8,5));plt.ylabel('Grouped Sobol index (polynomial supports)')
    _savefig(hi,'sobol_families')
    observed=bundle['val_y'][:,total];prediction=checks['prediction'][:,total]
    fig,ax=plt.subplots(figsize=(6,5));ax.errorbar(observed,prediction,xerr=bundle['val_se'][:,total],fmt='.',alpha=.6)
    lo=min(observed.min(),prediction.min());up=max(observed.max(),prediction.max());ax.plot([lo,up],[lo,up],'k--')
    ax.set(xlabel='OpenMC tritons/source neutron',ylabel='PCE prediction',title=f'Independent validation: R²={checks["r2"][total]:.4f}')
    _savefig(hi,'pce_validation')
    distribution={}
    if checks['accepted'][total]:
        samples=model.predict(np.random.default_rng(analysis_seed+91).random((distribution_samples,d)))[:,total]
        q05,q95=np.quantile(samples,[.05,.95])
        distribution=dict(mean=samples.mean(),standard_deviation=samples.std(),q05=q05,q95=q95)
        save_json(hi/'tritium_distribution.json',distribution)
        fig,ax=plt.subplots(figsize=(8,5));ax.hist(samples,bins=60,density=True,alpha=.6,label='Validated PCE')
        ax.plot(observed,np.zeros_like(observed),'k|',label='Independent OpenMC validation')
        for value,label in [(bundle['baseline'][total],'Baseline'),(samples.mean(),'PCE mean'),(q05,'5th percentile'),(q95,'95th percentile')]:ax.axvline(value,label=label)
        ax.set(xlabel=response_label,ylabel='Probability density',title=f'Conditional on retained hypothetical independent ±{100*prior_limit:g}% inputs')
        ax.legend(fontsize=8);_savefig(hi,'tritium_distribution')
    else:
        for suffix in ('png','pdf','json'):(hi/('tritium_distribution.'+suffix)).unlink(missing_ok=True)
    return checks,distribution
