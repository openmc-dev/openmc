"""Complete runtime reaction catalogue and plots before sensitivity selection."""
from dataclasses import asdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openmc
import openmc.lib as lib
from sandbox.analysis_codes.runner import Runner
from sandbox.analysis_codes.discovery import discover
from sandbox.analysis_codes.perturbations import shape
from sandbox.analysis_codes.storage import save_json
from openmc.lib.perturbation import _areas
from .context import build_model, responses as baseline_responses, register
from .thresholds import discrete_source_max_energy,filter_native_thresholds
from .native_plots import native_table,sampled_pdf_changed
from .analysis import _savefig,_markdown


def catalogue(output,make_plots=True,model_factory=build_model,response_factory=baseline_responses,spec_factory=register):
    """List every loaded XS and associated neutron law, retaining explicit status.

    Plots use native API data only. Deterministic outgoing-energy laws are
    probability masses; unsupported joint PDFs are recorded, never marginalised.
    """
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    plots=output/'all_retained_plots';plots.mkdir(exist_ok=True)
    openmc.reset_auto_ids()
    model=model_factory();source=discrete_source_max_energy(model)
    gallery=['# All source-retained reactions and neutron distributions\n',
             'This catalogue precedes Morris selection. Every source-retained cross section and associated law has an entry. Read-only laws are not mutated.\n']
    statuses=[]
    with Runner(spec_factory(model,response_factory(model)),output/'catalogue_state') as runner:
        from .context import settings
        params,inventory=discover(runner.spec.target_nuclides,limit=settings().PERTURBATION_LIMIT)
        retained=filter_native_thresholds(params,runner.spec.target_nuclides,source,output)
        runner.bind(retained)
        save_json(output/'catalogue_provenance.json',runner.provenance)
        save_json(output/'native_inventory.json',inventory)
        parameter_rows=[]
        for p in retained:
            active=True
            if p.family!='xs':
                baseline=runner.state.read(p)
                candidate=shape(baseline,[p],np.array([p.limit]))
                active=any(sampled_pdf_changed(a,b) for a,b in zip(baseline['tables'],candidate['tables']))
            parameter_rows.append(dict(**asdict(p),parameter_name=p.identity,native_mode_changes_pdf=active))
        pd.DataFrame(parameter_rows).to_csv(output/'retained_parameters.csv',index=False)
        save_json(output/'retained_parameters.json',parameter_rows)
        thresholds=pd.read_csv(output/'reaction_thresholds.csv')
        registry=[]
        for row in thresholds.itertuples():
            info=lib.get_reaction_info(row.nuclide,row.mt)
            below=row.threshold_eV<=source
            supports=any(p.nuclide==row.nuclide and p.mt==row.mt and p.family=='xs' for p in retained)
            item=dict(nuclide=row.nuclide,mt=row.mt,reaction=row.reaction,
                q_value_eV=info['q_value'],threshold_eV=row.threshold_eV,
                threshold_MeV=row.threshold_MeV,source_retained=below,
                mutable_xs=supports,redundant=info['redundant'],reference_frame=info['reference_frame'],
                data_units='eV*barn' if row.mt in (301,444,901) else 'barn',
                status=('Removed: threshold above source' if not below else
                        ('Retained: mutable' if supports else 'Retained: inspection only; dependent/unsupported mutation')))
            registry.append(item)
            if not below:continue
            base=f'{row.nuclide}_mt{row.mt}'
            gallery.append(f'## {row.nuclide} MT{row.mt} {row.reaction}\n\n{item["status"]}; threshold {row.threshold_MeV:.7g} MeV.\n')
            index=lib.nuclides[row.nuclide]._index
            grid=lib.nuclide_energy_grid(index);xs=lib.get_nuclide_xs(index,row.mt)
            e=grid[len(grid)-len(xs):]
            pd.DataFrame({'energy_eV':e,('energy_xs_eV_barn' if row.mt in (301,444,901) else 'xs_barn'):xs}).to_csv(plots/(base+'_xs.csv'),index=False)
            if make_plots:
                fig,ax=plt.subplots(figsize=(8,5));mask=(e>0)&(xs>0)
                ax.loglog(e[mask]/1e6,xs[mask]);ax.axvline(source/1e6,color='k',ls='--',label='Source maximum')
                ax.set(xlabel='Incident neutron energy [MeV]',ylabel=('Energy-production cross section [eV·barn]' if row.mt in (301,444,901) else 'Microscopic cross section [barn]'),title=f'{row.nuclide} MT{row.mt} {row.reaction}\n{item["status"]}')
                ax.legend();_savefig(plots,base+'_xs')
                gallery.append(f'![Cross section]({base}_xs.png)\n')
            neutron_products=[product for product in lib.get_reaction_products(row.nuclide,row.mt) if product['particle']=='neutron']
            if not neutron_products:
                gallery.append('No secondary-neutron product/law is stored for this reaction in the runtime data; there is no independent neutron energy/angular PDF to plot.\n')
            for product in neutron_products:
                pi=product['product_index']
                for descriptor in product['distributions']:
                    di=descriptor['distribution_index']
                    law=lib.get_emission_distribution(row.nuclide,row.mt,pi,di)
                    stem=f'{base}_p{pi}_d{di}'
                    save_json(plots/(stem+'_native.json'),law)
                    for family in ('energy','angle'):
                        component=law.get(family) or {}
                        representation=component.get('type',('absent' if law['type']=='uncorrelated' else law['type']))
                        mutable=any(p.nuclide==row.nuclide and p.mt==row.mt and p.product==pi and p.distribution==di and p.family==family for p in retained)
                        status=dict(nuclide=row.nuclide,mt=row.mt,product=pi,distribution=di,
                            family=family,representation=representation,mutable=mutable,
                            reference_frame=product['reference_frame'],plot='',reason='')
                        filename=stem+'_'+family
                        if representation in ('continuous-tabular','tabular-angle'):
                            incident,table=native_table(component,source)
                            key='energy' if family=='energy' else 'x';n=int(table.get('n_discrete',0))
                            x=np.asarray(table[key]);pdf=np.asarray(table['pdf']);interp=int(table['interpolation'])
                            if interp not in (1,2):raise ValueError('Unsupported native plotting interpolation')
                            status.update(selected_incident_energy_eV=incident,interpolation=interp,
                                normalization=float(pdf[:n].sum()+_areas(x[n:],pdf[n:],interp).sum()))
                            if make_plots:
                                fig,ax=plt.subplots(figsize=(8,5));scale=1e6 if family=='energy' else 1.
                                ax.plot(x[n:]/scale,pdf[n:]*scale,drawstyle='steps-post' if interp==1 else 'default')
                                if n:
                                    mass=ax.twinx();mass.stem(x[:n]/scale,pdf[:n]);mass.set_ylabel('Discrete probability mass')
                                ax.set(xlabel='Outgoing neutron energy [MeV]' if family=='energy' else 'Scattering cosine, μ = cos(θ)',
                                    ylabel='Probability density [1/MeV]' if family=='energy' else 'Probability density, p(μ)')
                                if family=='angle':ax.set_xlim(-1,1)
                        elif law['type']=='uncorrelated' and family=='energy' and representation=='level-inelastic':
                            energies=lib.sample_emission_energy(row.nuclide,row.mt,source,1,pi,di,seed=settings().NATIVE_PLOT_SEED if getattr(settings(),"NATIVE_PLOT_SEED",None) is not None else settings().SEED)
                            save_json(plots/(filename+'_mass.json'),dict(incident_energy_eV=source,energy_eV=energies,probability_mass=[1.]))
                            status['reason']='Read-only deterministic energy line; probability mass, not continuous PDF'
                            status['selected_incident_energy_eV']=source
                            if make_plots:
                                fig,ax=plt.subplots(figsize=(8,5));ax.stem(energies/1e6,[1.]);ax.set(xlabel='Outgoing neutron energy [MeV]',ylabel='Probability mass',ylim=(0,1.15))
                        elif family=='angle' and representation=='isotropic':
                            status['reason']='Native isotropic representation; no silent tabular conversion'
                            if make_plots:
                                fig,ax=plt.subplots(figsize=(8,5));ax.plot([-1,1],[.5,.5]);ax.set(xlabel='Scattering cosine, μ = cos(θ)',ylabel='Probability density, p(μ)',xlim=(-1,1))
                        else:
                            status['reason']='No independent PDF plotted: conditional/correlated, unsupported, or absent native component; see saved native dictionary'
                            gallery.append(f'**{family}: {representation}.** {status["reason"]}.\n')
                            statuses.append(status);continue
                        if make_plots:
                            ein=status.get('selected_incident_energy_eV',source)/1e6
                            ax.set_title(f'{row.nuclide} MT{row.mt} neutron {family}: {representation}\nEᵢₙ={ein:.6g} MeV; {product["reference_frame"]}; '+('mutable' if mutable else 'read only'))
                            ax.grid(alpha=.2);_savefig(plots,filename)
                            status['plot']=filename+'.png'
                            gallery.append(f'![{family} {representation}]({filename}.png)\n')
                        statuses.append(status)
        runner.state.verify()
    reaction_table=pd.DataFrame(registry)
    reaction_table.to_csv(output/'reaction_registry.csv',index=False)
    for name,subset in [('retained_cross_sections',reaction_table[reaction_table.source_retained]),
                        ('perturbable_cross_sections',reaction_table[reaction_table.source_retained & reaction_table.mutable_xs])]:
        subset.to_csv(output/(name+'.csv'),index=False)
        save_json(output/(name+'.json'),subset.to_dict('records'))
    save_json(output/'reaction_registry.json',registry)
    pd.DataFrame(statuses).to_csv(output/'distribution_registry.csv',index=False)
    (plots/'index.md').write_text('\n'.join(gallery))
    return pd.DataFrame(parameter_rows),source
