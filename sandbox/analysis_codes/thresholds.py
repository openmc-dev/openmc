"""Source-energy filtering of perturbation targets using loaded native data.

This removes parameters, not transport reactions or production tally scores.
The native threshold is the first grid energy for the loaded reaction/temperature.
"""
from collections import Counter
from dataclasses import asdict
from pathlib import Path
import numpy as np
import pandas as pd
import openmc
import openmc.lib as lib
from sandbox.analysis_codes.storage import save_json


def discrete_source_max_energy(model):
    """Return the largest positive-probability discrete source energy in eV.

    Continuous/file sources need an explicit finite cutoff; do not infer one
    from the distribution mean or from a handful of random samples.
    """
    maxima=[]
    for source in model.settings.source:
        if source.strength==0:
            continue
        energy=getattr(source,'energy',None)
        if not isinstance(energy,openmc.stats.Discrete):
            raise ValueError('Source-energy filtering requires a discrete source or an explicit maximum energy')
        x=np.asarray(energy.x,float);p=np.asarray(energy.p,float)
        if not np.any(p>0):raise ValueError('Source has no positive-probability energies')
        maxima.append(float(x[p>0].max()))
    if not maxima:raise ValueError('No active source')
    result=max(maxima)
    if not np.isfinite(result) or result<=0:raise ValueError('Source maximum must be finite and positive')
    return result


def partition_parameters(parameters, thresholds, maximum_energy):
    """Keep thresholds at/below the source maximum, for every parameter family."""
    if not np.isfinite(maximum_energy) or maximum_energy<=0:
        raise ValueError('Maximum source energy must be finite and positive')
    kept=[];removed=[]
    for parameter in parameters:
        key=(parameter.nuclide,parameter.mt)
        threshold=thresholds[key]  # Missing metadata is an error, never a guess.
        if not np.isfinite(threshold) or threshold<0:
            raise ValueError(f'Invalid native threshold for {key}')
        (removed if threshold>maximum_energy else kept).append(parameter)
    return kept,removed


def filter_native_thresholds(parameters, nuclides, maximum_energy, output):
    """Write a full native threshold audit and lists of removed reactions/modes."""
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    thresholds={};rows=[]
    counts=Counter((p.nuclide,p.mt) for p in parameters)
    family_counts=Counter((p.nuclide,p.mt,p.family) for p in parameters)
    for nuclide in nuclides:
        index=lib.nuclides[nuclide]._index
        for mt in lib.get_nuclide_mt_numbers(index):
            mt=int(mt);info=lib.get_reaction_info(nuclide,mt)
            threshold=float(lib.get_reaction_threshold_energy(index,mt))
            thresholds[(nuclide,mt)]=threshold
            count=counts[(nuclide,mt)]
            removed=threshold>maximum_energy and count>0
            rows.append(dict(nuclide=nuclide,mt=mt,reaction=info['name'],
                threshold_eV=threshold,threshold_MeV=threshold/1e6,
                source_max_eV=maximum_energy,source_max_MeV=maximum_energy/1e6,
                eligible_parameters_before_filter=count,
                xs_parameters=family_counts[(nuclide,mt,'xs')],
                energy_parameters=family_counts[(nuclide,mt,'energy')],
                angle_parameters=family_counts[(nuclide,mt,'angle')],
                removed_from_perturbations=removed,
                reason=('Native threshold exceeds source maximum' if removed else
                        ('Already outside safe mutation registry' if count==0 else 'Threshold at/below source maximum'))))
    kept,removed=partition_parameters(parameters,thresholds,maximum_energy)
    table=pd.DataFrame(rows).sort_values(['nuclide','mt'])
    table.to_csv(output/'reaction_thresholds.csv',index=False)
    excluded=table[table.removed_from_perturbations]
    excluded.to_csv(output/'removed_reactions.csv',index=False)
    pd.DataFrame([dict(**asdict(p),parameter_name=p.identity,
        threshold_eV=thresholds[(p.nuclide,p.mt)],source_max_eV=maximum_energy,
        reason='Native threshold exceeds source maximum') for p in removed],
        columns=[*asdict(parameters[0]).keys(),'parameter_name','threshold_eV','source_max_eV','reason'] if parameters else None
        ).to_csv(output/'removed_parameters.csv',index=False)
    save_json(output/'threshold_filter.json',dict(source_max_eV=maximum_energy,
        threshold_definition='Loaded native reaction xs threshold grid energy; temperature index 0',
        comparison='Remove only if threshold > source maximum; equality retained',
        parameters_before=len(parameters),parameters_after=len(kept),
        parameters_removed=len(removed),reactions_removed=len(excluded),
        kept_parameter_ids=[p.identity for p in kept],
        transport_reactions_deleted=False,production_tally_channels_deleted=False))
    return kept
