from dataclasses import asdict

from pathlib import Path

import time

import multiprocessing as mp

from concurrent.futures import ProcessPoolExecutor

import numpy as np

import pandas as pd

import matplotlib

import matplotlib.pyplot as plt

from scipy.stats import qmc

from sandbox.analysis_codes.runner import Runner

from sandbox.analysis_codes.discovery import discover

from sandbox.analysis_codes.morris import trajectories, effects, statistics, convergence

from sandbox.analysis_codes import pce, sobol

from sandbox.analysis_codes.validation import validate

from sandbox.analysis_codes.storage import save_json

from .settings import FULL
from .context import build_model, responses as baseline_responses, register
from .fingerprint import implementation_sha
_WORKER = None
_WORKER_LIMIT = None
_WORKER_NAMES = None
_WORKER_BATCHES = False
WARNING = "Independent bounded sensitivity priors; not evaluated covariance uncertainties."

def _initialize_worker(metadata, output, parameters, limit):
    global _WORKER, _WORKER_LIMIT, _WORKER_NAMES, _WORKER_BATCHES
    from multiprocessing.util import Finalize
    import openmc
    from sandbox.analysis_codes.config import ModelSpec
    _WORKER_BATCHES=metadata.pop('paired_morris',False)
    inputs = Path(metadata.pop('inputs'))
    model = openmc.Model.from_xml(
        geometry=inputs/'geometry.xml', materials=inputs/'materials.xml',
        settings=inputs/'settings.xml', tallies=inputs/'tallies.xml')
    cap=metadata.pop('history_budget')
    spec = ModelSpec(model=model, **metadata)
    from .budget import BudgetedRunner
    _WORKER = BudgetedRunner(spec, output, cap).__enter__()
    _WORKER.bind(parameters)
    # Native statepoint defaults can reappear during XML round trips. Isolate
    # every process's output files as well as its native library memory.
    import os
    os.chdir(_WORKER.temporary.name)
    _WORKER_LIMIT = limit
    _WORKER_NAMES = list(spec.responses)
    Finalize(None, _WORKER.__exit__, args=(None,None,None), exitpriority=10)


def _worker_evaluate(task):
    unit, stage, seed = task
    record = _WORKER.evaluate((2*np.asarray(unit)-1)*_WORKER_LIMIT,
                              stage=stage,seed=int(seed))
    result=(np.array([record['responses'][n][0][0] for n in _WORKER_NAMES]),
            np.array([record['responses'][n][1][0] for n in _WORKER_NAMES]))
    if _WORKER_BATCHES:
        batch=np.column_stack([np.asarray(record['batch_responses'][n])[:,0] for n in _WORKER_NAMES])
        return (*result,batch)
    return result


def paired_morris_errors(design,batches):
    """Per-input EE uncertainty retaining observed within-trajectory CRN covariance."""
    r,vertices,b,outputs=batches.shape
    if b<2:raise ValueError('At least two batches required for paired Morris')
    result=np.empty((r,vertices-1,outputs))
    for t in range(r):
        for k,j in enumerate(design.changed[t]):
            diff=(batches[t,k+1]-batches[t,k])/design.steps[t,k]
            result[t,j]=diff.std(axis=0,ddof=1)/np.sqrt(b)
    return result


def _label(identity):
    """Use readable labels; complete stable identities remain in CSV tables."""
    fields=str(identity).split(':')
    if len(fields)<6:
        return str(identity)
    name,mt,family,product,distribution,mode=fields[:6]
    if family=='xs':
        return f'{name} {mt} XS'
    return f'{name} {mt} {family} mode {mode[1:]} ({product}, {distribution})'


def _savefig(output, name):
    plt.tight_layout()
    plt.savefig(output / (name + '.png'), dpi=180)
    plt.savefig(output / (name + '.pdf'))
    plt.close()


def run_analysis(model=None, *, output, settings,
                 parameter_ids=None, screen_only=True,
                 source_energy_max_eV=None, response_factory=baseline_responses,
                 spec_factory=register, paired_morris=False, verification_directory=None,
                 response_mapping=None):
    """Run the low-fidelity Morris stage with explicit study settings."""
    import json
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    requested=cfg=settings
    if cfg.max_histories is None:
        raise ValueError('An explicit Morris history budget is required')
    no_full_budget=False;pilot=False;dry_run=False
    mapping=response_mapping
    if not screen_only:raise ValueError('Use high_transport for high-fidelity PCE designs')
    previous_config=output/'study_config.json'
    if previous_config.exists():
        previous=json.loads(previous_config.read_text())
        if previous.get('paired_morris',False)!=paired_morris:
            raise ValueError('Morris seed/covariance rule changed; choose a separate output')
        if previous.get('source_energy_max_eV')!=source_energy_max_eV:
            raise ValueError('Source threshold filter changed; choose a separate output')
        if previous.get('parameter_ids') != (list(parameter_ids) if parameter_ids is not None else None) or previous.get('screen_only',False)!=screen_only:
            raise ValueError('Parameter subset/stage changed; choose a separate output')
        if not dry_run:
            saved=previous['settings'];current=asdict(cfg)
            for key in ('particles','batches','trajectories','training','validation','seed','limit','active'):
                if saved[key]!=current[key]:
                    raise ValueError('Preserving existing corrected results: use a separate --output for changed settings')
    if model is None:model=build_model(cfg.particles,cfg.batches)
    spec=spec_factory(model,response_factory(model))
    histories=model.settings.particles*model.settings.batches
    full_hist=requested.particles*requested.batches
    if not dry_run:
        proof=(Path(verification_directory) if verification_directory is not None else output/'verification')/'verification.json'
        if not proof.exists() or json.loads(proof.read_text()).get('implementation_sha256')!=implementation_sha():
            print('Running controlled production/restoration verification before UQ',flush=True)
            if verification_directory is not None:raise ValueError('Supplied response verification is missing or stale')
            verify(output/'verification')
        if json.loads(proof.read_text())['status']!='passed':
            raise RuntimeError('Controlled verification failed; no UQ allowed')
    started=time.perf_counter()
    from .budget import BudgetedRunner
    with BudgetedRunner(spec,output/'checkpoints',cfg.max_histories) as runner:
        params,inventory=discover(spec.target_nuclides,limit=cfg.limit)
        if source_energy_max_eV is not None:
            from .thresholds import filter_native_thresholds
            params=filter_native_thresholds(params,spec.target_nuclides,source_energy_max_eV,output)
        if not params:raise ValueError('No parameters remain after source-energy filtering')
        if parameter_ids is not None:
            requested_ids=set(parameter_ids)
            available={p.identity:p for p in params}
            missing=requested_ids-available.keys()
            if missing:raise ValueError(f'Unknown refinement parameters: {missing}')
            params=[available[identity] for identity in parameter_ids]
        runner.bind(params);d=len(params);names=list(spec.responses)
        full_count=4+requested.trajectories*(d+1)+requested.training+requested.validation
        count=4+cfg.trajectories*(d+1)+cfg.training+cfg.validation
        cost=dict(eligible_parameters=d,maximum_selected_inputs=cfg.active,
            morris_trajectories=cfg.trajectories,morris_evaluations=cfg.trajectories*(d+1),
            training=cfg.training,validation=cfg.validation,
            histories_per_evaluation=histories,total_histories=count*histories,
            configured_budget=cfg.max_histories,full_budget_explicitly_configured=requested.max_histories is not None,
            full_initial_evaluations=full_count,full_initial_histories=full_count*full_hist,
            full_maximum_trajectories=requested.max_trajectories,
            full_configured_budget=requested.max_histories,
            candidate_terms={str(o):len(pce.indices(cfg.active,o)) for o in cfg.orders},
            execution_mode='pilot' if pilot or no_full_budget else 'full',
            additional_full_histories_if_unfunded=full_count*full_hist if no_full_budget else 0)
        save_json(output/('dry_run_cost_estimate.json' if dry_run else 'cost_estimate.json'),cost)
        print('COST '+__import__('json').dumps(cost),flush=True)
        if dry_run:return cost
        if cfg.max_histories is None or count*histories>cfg.max_histories:
            save_json(output/'budget_stop.json',dict(reason='Initial design exceeds explicit budget',cost=cost))
            raise RuntimeError('Initial design exceeds budget; no campaign launched. Use --pilot.')
        rows=[dict(**asdict(p),parameter_name=p.identity,eligible=True,included=True,
            perturbation_range=f'[-{p.limit:g},{p.limit:g}]',reason='Safe native perturbation; corrected lithium response') for p in params]
        excluded=[dict(nuclide=r['nuclide'],mt=r['mt'],family=r['family'],
            parameter_name='',eligible=False,included=False,reason=r.get('reason','Unsupported'),
            representation=r.get('representation','')) for r in inventory if not r.get('supported')]
        pd.DataFrame(rows+excluded).to_csv(output/'parameter_inventory.csv',index=False)
        save_json(output/'native_inventory.json',inventory)
        save_json(output/'study_config.json',dict(settings=asdict(cfg),requested_full_settings=asdict(requested),
            execution_mode=cost['execution_mode'],response_version='reaction_sum_lithium_v1',
            parameter_ids=list(parameter_ids) if parameter_ids is not None else None,screen_only=screen_only,
            source_energy_max_eV=source_energy_max_eV,paired_morris=paired_morris,
            response_order=names,mapping=mapping,provenance=runner.provenance,
            full_budget_missing=no_full_budget,prior=WARNING))
        save_json(output/'full_study_plan.json',dict(cost,morris_evaluations=cost['morris_evaluations']))
        from lxml import etree
        save_json(output/'model_details.json',dict(materials=[dict(name=m.name,
            density=m.density,density_units=m.density_units,temperature_K=m.temperature,
            nuclides=list(m.get_nuclides())) for m in model.geometry.get_all_materials().values()],
            sources=[etree.tostring(source.to_xml_element()).decode() for source in model.settings.source]))
        def evaluate(unit,stage,seed,*,with_timing=False):
            record=runner.evaluate((2*np.asarray(unit)-1)*cfg.limit,stage=stage,seed=int(seed))
            pair=(np.array([record['responses'][n][0][0] for n in names]),np.array([record['responses'][n][1][0] for n in names]))
            return (*pair,float(record['wall_seconds'])) if with_timing else pair
        baseline,base_se=evaluate(np.full(d,.5),'corrected_baseline',cfg.seed)
        save_json(output/'baseline.json',{n:dict(mean=baseline[i],standard_error=base_se[i]) for i,n in enumerate(names)})
        # A short timed pilot precedes the expensive design, with all data restored.
        timing=[]
        for j in range(3):
            point=np.full(d,.5);point[j%d]=.75
            _,_,duration=evaluate(point,'timing_pilot',cfg.seed+500+j,with_timing=True);timing.append(duration)
        cost.update(pilot_seconds_per_evaluation=float(np.median(timing)),
            estimated_full_serial_hours=float(np.median(timing)*(full_hist/histories)*full_count/3600),
            estimated_full_parallel_hours=float(np.median(timing)*(full_hist/histories)*(requested.trajectories*(d+1)/max(cfg.workers,1)+requested.training+requested.validation+4)/3600),
            timing_pilot_histories=3*histories)
        # Include three timing observations in the hard execution budget.
        if count*histories>cfg.max_histories:raise RuntimeError('Timing pilot leaves insufficient budget')
        save_json(output/'cost_estimate.json',cost)
        print('TIMED COST '+__import__('json').dumps(cost),flush=True)
        target=cfg.trajectories;maximum=cfg.max_trajectories
        plan=trajectories(d,maximum,levels=cfg.levels,seed=cfg.seed)
        reuse=False
        save_json(output/'design_reuse.json',dict(previous_design_coordinates_match=reuse,
            old_response_values_reused=False,transport_seed_rule=('seed+1000+trajectory*(d+1)' if paired_morris else 'seed+1000+trajectory*(d+1)+vertex')))
        means=np.empty((maximum,d+1,len(names)));errors=np.empty_like(means)
        batches=np.empty((maximum,d+1,cfg.batches,len(names))) if paired_morris else None
        done=0
        with ProcessPoolExecutor(max_workers=cfg.workers,mp_context=mp.get_context('spawn'),
                initializer=_initialize_worker,initargs=(dict(name=spec.name,responses=spec.responses,
                target_nuclides=spec.target_nuclides,threads=spec.threads,
                history_budget=cfg.max_histories,paired_morris=paired_morris,inputs=str((output/'checkpoints'/'model_inputs'/runner.provenance['model_sha256']).resolve())),
                output/'checkpoints',params,cfg.limit)) as executor:
            while done<target:
                t=done
                tasks=((plan.points[t,j],'corrected_morris',cfg.seed+1000+t*(d+1)+(0 if paired_morris else j)) for j in range(d+1))
                for j,result in enumerate(executor.map(_worker_evaluate,tasks)):
                    m,se=result[:2];means[t,j]=m;errors[t,j]=se
                    if paired_morris:batches[t,j]=result[2]
                    if j%50==0:
                        print(f'Morris {t*(d+1)+j+1}/{target*(d+1)}',flush=True)
                        save_json(output/'progress.json',dict(stage='morris',trajectory=t+1,vertex=j+1,target=target))
                done+=1
                np.savez_compressed(output/'morris_arrays.npz',points=plan.points[:done],changed=plan.changed[:done],
                    steps=plan.steps[:done],means=means[:done],errors=errors[:done],
                    **({'batches':batches[:done]} if paired_morris else {}))
                current=type(plan)(plan.points[:done],plan.changed[:done],plan.steps[:done],plan.levels,plan.seed)
                ee,ese=effects(current,means[:done],errors[:done])
                if paired_morris:ese=paired_morris_errors(current,batches[:done])
                if done>=2:
                    stat=statistics(ee,ese,bootstraps=cfg.progress_bootstraps,seed=cfg.bootstrap_seed if cfg.bootstrap_seed is not None else cfg.seed,confidence_level=cfg.confidence_level,noise_multiplier=cfg.noise_multiplier)
                    _write_morris(output,params,names,stat)
                    conv=convergence(ee,ese,minimum=cfg.minimum,rank_correlation=cfg.rank_correlation,top_overlap=cfg.top_overlap,max_relative_ci_width=cfg.max_relative_ci_width,bootstraps=cfg.convergence_bootstraps,seed=cfg.bootstrap_seed if cfg.bootstrap_seed is not None else cfg.seed,confidence_level=cfg.confidence_level,noise_multiplier=cfg.noise_multiplier);save_json(output/'morris_convergence.json',conv)
                print(f'Completed trajectory {done}/{target}',flush=True)
                if done==target:
                    converged=bool(conv[-1]['stopping_criterion_met']) if done>=2 else False
                    prospective=target+cfg.trajectory_increment
                    required=(1+prospective*(d+1)+cfg.training+cfg.validation+3)*histories
                    if not converged and prospective<=maximum and required<=cfg.max_histories:
                        target=prospective
                    else:
                        save_json(output/'morris_stopping.json',dict(trajectories=done,
                            converged=converged,reason=('converged' if converged else ('maximum configured trajectories'
                                if prospective>maximum else 'history budget after reserving PCE/validation')),
                            next_increment_histories=required,configured_budget=cfg.max_histories))
        if paired_morris:
            _,independent_error=effects(current,means[:done],errors[:done])
            save_json(output/'paired_noise_diagnostics.json',dict(
                median_paired_to_unpaired_se=float(np.median(ese/np.maximum(independent_error,1e-300))),
                median_paired_to_unpaired_se_total=float(np.median(ese[:,:,names.index('total')]/np.maximum(independent_error[:,:,names.index('total')],1e-300))),
                rule='same seed within trajectory; independent across trajectories; actual batch-difference SE',
                total_parameter_mean_paired_se=ese[:,:,names.index('total')].mean(axis=0)))
        stat=statistics(ee,ese,bootstraps=cfg.bootstraps,seed=cfg.bootstrap_seed if cfg.bootstrap_seed is not None else cfg.seed,confidence_level=cfg.confidence_level,noise_multiplier=cfg.noise_multiplier);_write_morris(output,params,names,stat)
        if screen_only:
            save_json(output/'execution.json',dict(elapsed_seconds=time.perf_counter()-started,
                transport_evaluations=4+done*(d+1),total_histories=(4+done*(d+1))*histories,
                mode='low_statistics_screen',restoration='Exact native rollback every evaluation'))
            return pd.read_csv(output/'morris_results.csv')


def _write_morris(output,params,names,stat):
    rows=[]
    for i,p in enumerate(params):
        for j,n in enumerate(names):
            row=dict(parameter=p.identity,nuclide=p.nuclide,mt=p.mt,family=p.family,response=n)
            row.update({k:float(v[i,j]) for k,v in stat.items()
                        if isinstance(v,np.ndarray) and v.shape==(len(params),len(names))})
            rows.append(row)
    pd.DataFrame(rows).to_csv(output/'morris_results.csv',index=False)


def _markdown(frame):
    """Render small tables without pandas' optional tabulate dependency."""
    headers = [str(x) for x in frame.columns]
    rows = ['| ' + ' | '.join(headers) + ' |',
            '| ' + ' | '.join(['---']*len(headers)) + ' |']
    for values in frame.itertuples(index=False,name=None):
        rows.append('| ' + ' | '.join((f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x)).replace('|','/') for x in values) + ' |')
    return '\n'.join(rows)
