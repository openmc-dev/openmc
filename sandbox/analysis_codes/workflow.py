"""Study-independent orchestration of discovery, Morris, native runs and PCE."""
from pathlib import Path
from dataclasses import replace
import ast
import json
import numpy as np
import pandas as pd
import openmc
from . import context, pce
from .storage import save_json
from .fingerprint import implementation_sha
from .catalogue import catalogue
from .analysis import run_analysis
from .settings import FULL
from .family_workflow import select_per_family, high_transport, analyse_saved
from .two_pass import screen_plots, write_sensitivity_index
from .native_plots import plot_ranked_native
from .production import prepare, verify, sha


def guard_request(study):
    """Protect existing observations from changed transport inputs.

    PCE options are deliberately excluded: they consume saved transport data.
    The previous completed FLiBe study can be adopted after explicit input checks.
    Partially completed legacy campaigns require their original implementation.
    """
    c=study.config;output=study.output;output.mkdir(parents=True,exist_ok=True)
    excluded={key for key in dir(c) if key.startswith('PCE_')}
    excluded.add('ANALYSIS_SEED')
    request=dict(config={key:getattr(c,key) for key in dir(c)
                         if key.isupper() and key not in excluded},
                 model_sha256=sha(study.directory/'model.py'),
                 engine_sha256=implementation_sha())
    request=json.loads(json.dumps(request))
    file=output/'analysis_request.json'
    if file.exists():
        previous=json.loads(file.read_text())
        if previous!=request:
            raise ValueError('Changed model/settings/engine require a new --output; --reanalyse uses saved data')
        return
    legacy=output/'request.json'
    if legacy.exists():
        old=json.loads(legacy.read_text())
        if not all(request['config'].get(k)==v for k,v in old['config'].items()
                   if k not in excluded):
            raise ValueError('Legacy configuration changed: use a new --output')
        archived=output/'source_snapshot/before_analysis_codes/model.py'
        def builder(file):
            return ast.dump(next(n for n in ast.parse(file.read_text()).body
                if isinstance(n,ast.FunctionDef) and n.name=='build_model'),include_attributes=False)
        if not archived.exists() or builder(archived)!=builder(study.directory/'model.py'):
            raise ValueError('Legacy model changed: use a new --output')
        if not all((output/f).exists() for f in ('low/execution.json','high/execution.json','high/analysis_arrays.npz')):
            raise ValueError('Partially complete legacy run: resume with original archived engine or a new --output')
        save_json(output/'legacy_adoption.json',dict(
            original_request_sha256=sha(legacy),model_builder_ast_identical=True,
            original_transport_configuration_identical=True,
            transport_results_reused=True,transport_provenance_preserved=True))
    save_json(file,request)


def fit_saved(study):
    """Select the configured backend and analyse identical saved observations."""
    c=study.config
    return analyse_saved(study.output,fit_options=dict(
        q=c.PCE_HYPERBOLIC_Q,interaction=c.PCE_MAX_INTERACTION,
        ridge=c.PCE_RIDGE,sparse=c.PCE_SPARSE,
        keep_main_effects=c.PCE_KEEP_MAIN_EFFECTS,basis_limit=c.PCE_BASIS_LIMIT),
        response_label=c.RESPONSE_LABEL,analysis_seed=c.ANALYSIS_SEED,
        backend=c.PCE_BACKEND,orders=c.PCE_ORDERS,folds=c.PCE_CV_FOLDS,prior_limit=c.PERTURBATION_LIMIT,
        distribution_samples=c.PCE_DISTRIBUTION_SAMPLES,validation_options=dict(min_r2=getattr(c,"PCE_MIN_R2",.9),max_noise_mse=getattr(c,"PCE_MAX_NOISE_MSE",2.),min_signal_to_noise=getattr(c,"PCE_MIN_SIGNAL_TO_NOISE",3.)))


def report(study,checks,distribution):
    """Write a compact study-specific guide to saved numerical results."""
    out=study.output;c=study.config
    selected=pd.read_csv(out/'selected_parameters.csv')
    conv=json.loads((out/'low/morris_convergence.json').read_text())[-1]
    write_sensitivity_index(out)
    text=f'''# {c.STUDY_NAME} sensitivity analysis

Model definition: `../model.py`; all run settings: `../config.py`.
Response: {c.RESPONSE_LABEL}. Reaction contributions are summed **within batches**
before estimating uncertainties; cached MT205 is not the primary response.

Morris: {c.MORRIS_TRAJECTORIES} trajectories; {c.SCREEN_PARTICLES*c.SCREEN_BATCHES:,}
histories per evaluation. Convergence criterion met: **{conv['stopping_criterion_met']}**.
Unresolved effects are retained with flags; screening does not prove they are zero.

Selected parameters per family: **{selected.groupby('family').size().to_dict()}**.
High-fidelity: {c.TRAINING_SAMPLES} training, {c.VALIDATION_SAMPLES} independent
validation points, {c.HIGH_PARTICLES*c.HIGH_BATCHES:,} histories per evaluation.
PCE backend: **{c.PCE_BACKEND}**, orders {c.PCE_ORDERS}, training-only CV.
Total response held-out R²: **{checks['r2'][0]:.6g}**;
RMSE: {checks['rmse'][0]:.6g}; accepted: **{checks['accepted'][0]}**.

[All Morris/Sobol figures](sensitivity_plots.md) ·
[All retained native nuclear-data plots](all_retained_plots/index.md) ·
[Selected native perturbations](nuclear_data_plots/index.md) ·
[Removed above-source reactions](removed_reactions.csv) ·
[Every parameter](retained_parameters.csv) ·
[Production verification](verification/controlled_tests.csv).

![PCE validation](high/pce_validation.png)
![Every selected Sobol index](high/sobol_indices.png)
![Family indices](high/sobol_families.png)

Response distribution: `{distribution}`.
{'' if not distribution else '![Production probability distribution](high/tritium_distribution.png)'}

Independent hypothetical ±{100*c.PERTURBATION_LIMIT:g}% coefficient priors are **not evaluated covariance**.
XS scale directly; PDF modes share a two-mode bound (one coefficient alone has
at most ±{50*c.PERTURBATION_LIMIT:g}% local PDF change). Native support/interpolation is preserved.
Sobol results and the response distribution are conditional on selected inputs.
Weak components remain marked unresolved. A validated total response does not
establish full-space screening convergence or absence of interactions.

```bash
python -m sandbox.analysis_codes {study.directory} --dry-run
python -m sandbox.analysis_codes {study.directory} --run
python -m sandbox.analysis_codes {study.directory} --resume
python -m sandbox.analysis_codes {study.directory} --reanalyse
```

Completed legacy transport results, original reports and provenance remain intact.
PCE backend reports are also preserved under `pce_backends/<backend>/`.
'''
    (out/'summary.md').write_text(text)


def run(study,*,dry_run=False,reanalyse=False):
    """Execute only the stages required for the selected study and mode."""
    c=study.config;out=study.output;out.mkdir(parents=True,exist_ok=True)
    for key in ('SCREEN_PARTICLES','HIGH_PARTICLES','WORKERS','TRANSPORT_THREADS','TRAINING_SAMPLES','VALIDATION_SAMPLES','TOP_PER_FAMILY'):
        if getattr(c,key)<1:raise ValueError(f'{key} must be positive')
    if getattr(c,'HIGH_WORKERS',1)<1:
        raise ValueError('HIGH_WORKERS must be positive')
    if c.MORRIS_TRAJECTORIES<2 or min(c.SCREEN_BATCHES,c.HIGH_BATCHES,c.VERIFICATION_BATCHES)<2:
        raise ValueError('At least two trajectories and two batches are required')
    if c.PCE_BACKEND not in ('numpy','openturns'):
        raise ValueError('PCE_BACKEND must be numpy or openturns')
    if c.PCE_BACKEND=='openturns':
        from .openturns_pce import require_openturns
        require_openturns()
    if not np.isfinite(c.PERTURBATION_LIMIT) or not 0 < c.PERTURBATION_LIMIT < 1:
        raise ValueError("PERTURBATION_LIMIT must be finite and strictly between 0 and 1")
    if not reanalyse:
        guard_request(study)
        if not hasattr(study.model,'responses'):
            prepare(out)
        hooks=dict(model_factory=context.build_model,response_factory=context.responses,
                   spec_factory=context.register)
        completed=(out/'legacy_adoption.json').exists()
        if completed:
            retained=pd.read_csv(out/'retained_parameters.csv')
            source=json.loads((out/'threshold_filter.json').read_text())['source_max_eV']
        elif dry_run:
            from tempfile import TemporaryDirectory
            with TemporaryDirectory() as temporary:
                retained,source=catalogue(temporary,make_plots=False,**hooks)
        else:
            retained,source=catalogue(out,make_plots=True,**hooks)
        d=len(retained)
        low=(4+c.MORRIS_TRAJECTORIES*(d+1))*c.SCREEN_PARTICLES*c.SCREEN_BATCHES
        high=(1+c.TRAINING_SAMPLES+c.VALIDATION_SAMPLES)*c.HIGH_PARTICLES*c.HIGH_BATCHES
        verification=(2+sum(n in c.PRODUCTION_CHANNELS for n in ('Li6','Li7','Be9','Be9','F19','Li7')))*c.VERIFICATION_PARTICLES*c.VERIFICATION_BATCHES
        plan=dict(parameters=d,parameters_by_family=retained.groupby('family').size().to_dict(),
            morris_evaluations=4+c.MORRIS_TRAJECTORIES*(d+1),morris_histories=low,
            high_evaluations=1+c.TRAINING_SAMPLES+c.VALIDATION_SAMPLES,high_histories=high,
            verification_histories=verification,total_histories=low+high+verification,
            history_cap=c.HISTORY_CAP,pce_backend=c.PCE_BACKEND,perturbation_limit=c.PERTURBATION_LIMIT)
        if low>c.LOW_HISTORY_CAP or high>c.HIGH_HISTORY_CAP or verification>c.VERIFICATION_HISTORY_CAP or low+high+verification>c.HISTORY_CAP or c.LOW_HISTORY_CAP+c.HIGH_HISTORY_CAP+c.VERIFICATION_HISTORY_CAP>c.HISTORY_CAP:
            raise ValueError('Planned campaign exceeds configured history budgets')
        save_json(out/('dry_run_plan.json' if dry_run else 'analysis_plan.json'),plan)
        print(json.dumps(plan,indent=2),flush=True)
        if dry_run:return plan
        if not completed:
            if hasattr(study.model,'responses'):
                if not hasattr(study.model,'verify_response'):
                    raise ValueError('Custom model.responses requires model.verify_response(output) before transport')
                study.model.verify_response(out/'verification')
            else:
                proof=out/'verification/verification.json'
                if not proof.exists() or json.loads(proof.read_text()).get('implementation_sha256')!=implementation_sha():
                    verify(out/'verification')
            cfg=replace(FULL,particles=c.SCREEN_PARTICLES,batches=c.SCREEN_BATCHES,
                trajectories=c.MORRIS_TRAJECTORIES,max_trajectories=c.MORRIS_TRAJECTORIES,
                levels=c.MORRIS_LEVELS,seed=c.SEED,active=c.TOP_PER_FAMILY,limit=c.PERTURBATION_LIMIT,orders=c.PCE_ORDERS,
                training=0,validation=0,max_histories=c.LOW_HISTORY_CAP,workers=c.WORKERS,
                minimum=getattr(c,'MORRIS_MIN_TRAJECTORIES',20),rank_correlation=getattr(c,'MORRIS_RANK_CORRELATION',.9),
                top_overlap=getattr(c,'MORRIS_TOP_OVERLAP',.8),max_relative_ci_width=getattr(c,'MORRIS_MAX_RELATIVE_CI_WIDTH',.5),
                bootstraps=getattr(c,'MORRIS_BOOTSTRAPS',1000),progress_bootstraps=getattr(c,'MORRIS_PROGRESS_BOOTSTRAPS',300),
                convergence_bootstraps=getattr(c,'MORRIS_CONVERGENCE_BOOTSTRAPS',200),confidence_level=getattr(c,'MORRIS_CONFIDENCE_LEVEL',.95),
                noise_multiplier=getattr(c,'MORRIS_NOISE_MULTIPLIER',2.),bootstrap_seed=getattr(c,'MORRIS_BOOTSTRAP_SEED',None))
            openmc.reset_auto_ids();model=context.build_model(cfg.particles,cfg.batches)
            run_analysis(model=model,output=out/'low',settings=cfg,screen_only=True,
                source_energy_max_eV=source,response_factory=context.responses,
                spec_factory=context.register,paired_morris=c.PAIRED_MORRIS,
                verification_directory=out/'verification',response_mapping=json.loads((out/'provenance/production_mapping.json').read_text()) if (out/'provenance/production_mapping.json').exists() else {})
            morris=pd.read_csv(out/'low/morris_results.csv')
            selected=select_per_family(morris,retained,c.TOP_PER_FAMILY)
            checks=morris[morris.response=='total'];extras=[]
            for n,mt,family,mode in c.RESERVED_PARAMETERS:
                ids=retained[(retained.nuclide==n)&(retained.mt==mt)&(retained.family==family)&(retained['mode']==mode)].parameter_name
                for identity in ids:
                    if identity not in set(selected.parameter):
                        row=checks[checks.parameter==identity].copy();row['family_rank']=0
                        row['selection_reason']='Configured reserved physics check';extras.append(row)
            if extras:selected=pd.concat([selected,*extras],ignore_index=True)
            selected.to_csv(out/'selected_parameters.csv',index=False)
            checks=checks.copy();checks['selected']=checks.parameter.isin(selected.parameter)
            checks.to_csv(out/'parameter_decisions.csv',index=False)
            high_transport(out,selected,json.loads((out/'retained_parameters.json').read_text()),
                training=c.TRAINING_SAMPLES,validation=c.VALIDATION_SAMPLES,
                particles=c.HIGH_PARTICLES,batches=c.HIGH_BATCHES,seed=c.SEED,
                history_cap=c.HIGH_HISTORY_CAP,response_label=c.RESPONSE_LABEL,workers=getattr(c,"HIGH_WORKERS",1),**hooks)
            plot_ranked_native(out,parameter_list=selected,per_family=len(selected),**hooks)
    screen_plots(out/'low')
    checks,distribution=fit_saved(study)
    report(study,checks,distribution)
    # Preserve both backends' derived artifacts while retaining simple current paths.
    import shutil
    archive=out/'pce_backends'/c.PCE_BACKEND;archive.mkdir(parents=True,exist_ok=True)
    for pattern in ('pce_*','sobol_*','tritium_distribution.*'):
        for file in (out/'high').glob(pattern):
            if file.is_file():shutil.copy2(file,archive/file.name)
    shutil.copy2(out/'summary.md',archive/'summary.md')
    save_json(archive/'configuration.json',{k:getattr(c,k) for k in dir(c) if k.startswith('PCE_')})
    print(out/'summary.md',flush=True)
    return checks,distribution
