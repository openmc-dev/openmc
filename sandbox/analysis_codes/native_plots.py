"""Plot ranked native data before/after the existing transactional perturbations."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openmc.lib as lib
from openmc.lib.perturbation import _areas
from sandbox.analysis_codes.discovery import discover
from sandbox.analysis_codes.storage import save_json
from .context import build_model, responses as baseline_responses, register
from .budget import BudgetedRunner
from .analysis import _savefig,_label


def native_table(law, target=None):
    """Select a native incident table; never interpolate into a new PDF law."""
    if target is None:
        from .context import settings
        target=settings().SOURCE_ENERGY_EV
    j=int(np.argmin(abs(np.asarray(law['incident_energy'])-target)))
    return float(law['incident_energy'][j]),law['tables'][j]


def sampled_pdf_changed(baseline, candidate):
    """Ignore the unused terminal ordinate for a native histogram PDF."""
    n=int(baseline.get('n_discrete',0))
    stop=-1 if int(baseline['interpolation'])==1 else None
    return bool(np.any(np.asarray(candidate['pdf'])[n:stop] != np.asarray(baseline['pdf'])[n:stop]))


def plot_ranked_native(output,low=None,high=None,per_family=10,parameter_list=None,
                       model_factory=build_model,response_factory=baseline_responses,spec_factory=register):
    """Read back configured endpoint candidates through native APIs without running transport.

    Energy/angular rankings may be unresolved. Their labels report that status;
    plotting a PDF shape is not evidence for a significant transport response.
    """
    output=Path(output);destination=output/'nuclear_data_plots'
    destination.mkdir(exist_ok=True)
    if parameter_list is None:
        lowrank=pd.read_csv(Path(low)/'morris_results.csv')
        highrank=pd.read_csv(Path(high)/'morris_results.csv')
    else:
        lowrank=parameter_list.copy()
        if 'parameter' not in lowrank:lowrank['parameter']=lowrank['parameter_name']
        lowrank['response']='total'
        for column,value in [('resolved',False),('signal_to_noise',0.),('noise_corrected_rms',0.)]:
            if column not in lowrank:lowrank[column]=value
        highrank=lowrank
    from .context import settings
    limit=settings().PERTURBATION_LIMIT
    labels=['Baseline',f'Coefficient −{limit:g}',f'Coefficient +{limit:g}']
    gallery=['# Ranked native nuclear data\n',
        f'Baseline and native readback after individual coefficient changes of −{limit:g} and +{limit:g}. '
        'Shapes share the existing two-mode budget; all other coefficients are zero. '
        'Unresolved rankings are provisional. CSV/JSON preserve the numerical data.\n']
    model=model_factory()
    with BudgetedRunner(spec_factory(model,response_factory(model)),
                        output/'plot_checkpoints',1) as runner:
        params,_=discover(runner.spec.target_nuclides,limit=limit);runner.bind(params)
        lookup={p.identity:(i,p) for i,p in enumerate(params)}
        diagnostics=[]
        for family in ('xs','energy','angle'):
            ranked=highrank[(highrank.response=='total')&(highrank.family==family)]
            stage='High-fidelity PCE selected input; low-history Morris ranking' if parameter_list is not None else 'high-statistics reduced-space Morris'
            if ranked.empty:
                ranked=lowrank[(lowrank.response=='total')&(lowrank.family==family)]
                stage='low-statistics full-space Morris'
            ranked=(ranked.sort_values('family_rank') if parameter_list is not None and 'family_rank' in ranked
                    else ranked.sort_values('noise_corrected_rms',ascending=False)).head(per_family).copy()
            ranked['plot_ranking_source']=stage
            if family=='angle' and parameter_list is None:
                # Coarse histogram tables can make a mode null. Include native
                # elastic P1 references, clearly separated from refined ranks.
                elastic=lowrank[(lowrank.response=='total')&(lowrank.family=='angle')&
                    (lowrank.mt==2)&lowrank.parameter.str.contains(':m1:')].copy()
                elastic['plot_ranking_source']='low-statistics full-space Morris; elastic reference'
                ranked=pd.concat([ranked,elastic]).drop_duplicates('parameter')
            for rank,row in enumerate(ranked.itertuples(),1):
                stage=row.plot_ranking_source
                i,p=lookup[row.parameter];original=runner.state.read(p)
                laws=[original]
                # Values come from C++ after actual native setters, not merely
                # from the Python candidate construction.
                for delta in (-p.limit,p.limit):
                    values=np.zeros(len(params));values[i]=delta
                    candidate=runner.state.candidates(values)
                    with lib.nuclear_data_perturbation(**candidate):
                        runner.state.verify_candidates(candidate)
                        laws.append(runner.state.read(p))
                    runner.state.verify()
                basename=f'{family}_rank{rank}_{p.nuclide}_mt{p.mt}_p{p.product}_d{p.distribution}_mode{p.mode}'
                fig,ax=plt.subplots(figsize=(9,5))
                metadata=dict(parameter=p.identity,ranking_stage=stage,resolved=bool(row.resolved),
                    signal_to_noise=float(row.signal_to_noise),delta=[0.,-p.limit,p.limit],
                    native_readback=True,restored_exactly=True)
                if family=='xs':
                    grid=runner.state.grids[p.target]
                    for law,label in zip(laws,labels):
                        ax.loglog(grid/1e6,law,label=label)
                    ax.set(xlabel='Incident neutron energy [MeV]',ylabel='Microscopic cross section [barn]')
                    pd.DataFrame(dict(incident_energy_eV=grid,baseline=laws[0],minus=laws[1],plus=laws[2])).to_csv(destination/(basename+'.csv'),index=False)
                    # A second panel makes the fractional change readable.
                    fig2,ax2=plt.subplots(figsize=(9,4))
                    mask=laws[0]>0
                    for law,label in zip(laws[1:],labels[1:]):
                        ax2.semilogx(grid[mask]/1e6,law[mask]/laws[0][mask],label=label)
                    ax2.set(xlabel='Incident neutron energy [MeV]',ylabel='Perturbed / baseline cross section',title=_label(p.identity))
                    ax2.legend();plt.figure(fig2.number);_savefig(destination,basename+'_ratio')
                    plt.figure(fig.number)
                else:
                    norms=[];saved=[]
                    for law,label in zip(laws,labels):
                        incident,table=native_table(law)
                        key='energy' if family=='energy' else 'x'
                        x=np.asarray(table[key]);pdf=np.asarray(table['pdf'])
                        n=int(table.get('n_discrete',0));interp=int(table['interpolation'])
                        if interp not in (1,2):raise ValueError(f'Unsupported plot interpolation {interp}')
                        norm=float(pdf[:n].sum()+_areas(x[n:],pdf[n:],interp).sum())
                        norms.append(norm)
                        scale=1e6 if family=='energy' else 1.
                        ax.plot(x[n:]/scale,pdf[n:]*scale,label=label,
                                drawstyle='steps-post' if interp==1 else 'default')
                        if n:
                            # Discrete probabilities have a separate axis;
                            # do not display masses as continuous densities.
                            massax=ax.twinx()
                            massax.stem(x[:n]/scale,pdf[:n],basefmt=' ')
                            massax.set_ylabel('Discrete probability mass')
                        saved.append(dict(incident_energy_eV=incident,table=table))
                    if family=='energy':
                        ax.set(xlabel='Outgoing neutron energy [MeV]',ylabel='Probability density [1/MeV]')
                    else:
                        ax.set(xlabel='Scattering cosine, μ = cos(θ)',ylabel='Probability density, p(μ)',xlim=(-1,1))
                    products=lib.get_reaction_products(p.nuclide,p.mt)
                    product=next(q for q in products if q['product_index']==p.product)
                    baseline_table=saved[0]['table']
                    changed=any(sampled_pdf_changed(baseline_table,t['table']) for t in saved[1:])
                    metadata['sampled_continuous_pdf_changed']=bool(changed)
                    if not changed:
                        metadata['limitation']='Native grid does not resolve this mode: the sampled continuous PDF is unchanged; an unused histogram endpoint may differ.'
                    metadata.update(incident_energy_eV=incident,interpolation=interp,
                        integrals=norms,reference_frame=product['reference_frame'],
                        emitted_particle=product['particle'],representation=original['type'])
                    save_json(destination/(basename+'_tables.json'),saved)
                ax.set_title(_label(p.identity)+'\n'+stage+
                             ('; resolved' if row.resolved else '; unresolved')+
                             (f"; Eᵢₙ={metadata['incident_energy_eV']/1e6:.4g} MeV; {metadata['reference_frame']}" if family!='xs' else ''))
                ax.legend();ax.grid(alpha=.2)
                _savefig(destination,basename)
                save_json(destination/(basename+'_metadata.json'),metadata)
                diagnostics.append(metadata)
                gallery.append(f'## {_label(p.identity)}\n\n'
                    f"{stage}; {'resolved' if row.resolved else 'unresolved'}; S/N={row.signal_to_noise:.3g}.\n\n"
                    f'![Native {family}]({basename}.png)\n')
                if 'limitation' in metadata:gallery.append('**Native representation limitation:** '+metadata['limitation']+'\n')
                if family=='xs':gallery.append(f'![Cross-section ratio]({basename}_ratio.png)\n')
        save_json(destination/'diagnostics.json',diagnostics)
    (destination/'index.md').write_text('\n'.join(gallery))
