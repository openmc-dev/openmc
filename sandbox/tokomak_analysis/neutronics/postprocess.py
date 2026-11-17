"""Structured results preserve batch covariance and physical source units."""
from pathlib import Path
import json
import numpy as np
import pandas as pd


def summarize_batches(batches,c):
    names=['Li6','Li7','Be9','lithium_total','total']
    values=np.column_stack([batches[:,0],batches[:,1],batches[:,2],batches[:,:2].sum(axis=1),batches.sum(axis=1)])
    means=values.mean(axis=0);se=values.std(axis=0,ddof=1)/np.sqrt(len(values))
    return {name:dict(per_source=float(mean),standard_error=float(error),per_second=float(mean*c.neutron_rate),
                      per_second_standard_error=float(error*c.neutron_rate)) for name,mean,error in zip(names,means,se)}


def write_outputs(a,directory,batches,performance):
    directory=Path(directory);rates=summarize_batches(batches,a.config)
    result=dict(lithium_tritium_production={k:rates[k] for k in ('Li6','Li7','lithium_total')},verified_tritium_production=rates,lithium_TBR=rates['lithium_total'],TBR=rates['total'],neutron_rate_per_second=a.config.neutron_rate,
      fusion_power_mw=a.config.fusion_power_mw,source_neutron_power_mw=a.config.neutron_rate*a.config.source_energy_ev*1.602176634e-19/1e6,
      response_scope='Verified Li6 MT105 + Li7 MT52..82 + Be9 MT700/701; structural/oxygen triton-production side channels are not included or asserted zero',
      heating_mode='coupled neutron-photon' if a.config.coupled_photons else 'neutron heating-local approximation',performance=performance)
    (directory/'results.json').write_text(json.dumps(result,indent=2));pd.DataFrame(batches,columns=['Li6','Li7','Be9']).to_csv(directory/'batch_tritium.csv',index=False)
    return result



def export_native_tallies(a,out):
    import openmc.lib as lib
    manifest=[]
    for tally in a.model.tallies:
        native=lib.tallies[tally.id]
        name=tally.name.lower().replace(' ','_').replace('/','_')
        mean=np.asarray(native.mean);se=np.asarray(native.std_dev)
        np.savez_compressed(Path(out)/(name+'.npz'),mean=mean,standard_error=se)
        units='eV/source' if 'heating' in name else ('reactions/source' if 'production' in name or 'n2n' in name else ('neutrons/source' if 'leakage' in name else 'cm/source (volume-integrated flux)'))
        manifest.append(dict(tally_id=tally.id,name=tally.name,file=name+'.npz',shape=mean.shape,units=units,
            scores=tally.scores,nuclides=tally.nuclides,filters=[dict(type=type(f).__name__,bins=np.asarray(f.bins).tolist()) for f in tally.filters]))
    (Path(out)/'tally_manifest.json').write_text(json.dumps(manifest,indent=2))


def write_module_results(a,out,bins,batches,total_batches):
    import re
    from ..visualization.export_vtk import path_metadata
    # Parent paths retain sector and local module identities through repetition.
    rows=[]
    for col,(cell_id,instance) in enumerate(bins):
        cell=a.registry.cells[int(cell_id)];path=cell.paths[int(instance)]
        ids=[int(i) for i in re.findall(r'c(\d+)',path)]
        data=[a.registry.records[i] for i in ids]
        sector=next(r['sector'] for r in data if 'sector' in r)
        local=next(r['toroidal_local'] for r in data if 'toroidal_local' in r)
        pol=a.registry.records[int(cell_id)]['poloidal']
        mean=batches[:,col].mean();se=batches[:,col].std(ddof=1)/np.sqrt(len(batches))
        rows.append(dict(module_id=(sector*a.config.toroidal_modules_per_sector+local)*a.config.poloidal_modules+pol,
          sector=sector,toroidal_local=local,poloidal=pol,inboard=bool(a.registry.records[int(cell_id)]['inboard']),
          tritons_per_source=mean,standard_error=se,tritons_per_second=mean*a.config.neutron_rate))
    table=pd.DataFrame(rows);table.to_csv(Path(out)/'module_breeding.csv',index=False)
    np.savez_compressed(Path(out)/'module_batch_production.npz',values=batches,bins=np.asarray(bins))
    # This catches incorrect parent-instance tally indexing and lost module bins.
    np.testing.assert_allclose(batches.sum(axis=1),total_batches.sum(axis=1),rtol=1e-9,atol=1e-10)
    sides={}
    for name,mask in [('inboard',table.inboard.to_numpy()),('outboard',~table.inboard.to_numpy())]:
        values=batches[:,mask].sum(axis=1)
        sides[name]=dict(tritons_per_source=float(values.mean()),standard_error=float(values.std(ddof=1)/np.sqrt(len(values))))
    (Path(out)/'inboard_outboard_breeding.json').write_text(json.dumps(sides,indent=2))



def write_integrated_results(a,out,batches):
    rows=[]
    for j,t in enumerate(a.model.tallies):
        values=batches[:,j];mean=values.mean();se=values.std(ddof=1)/np.sqrt(len(values))
        heating='heating' in t.name.lower()
        rows.append(dict(tally=t.name,mean_per_source=mean,standard_error=se,
          source_scaled_mean=mean*a.config.neutron_rate,
          heating_watts=mean*a.config.neutron_rate*1.602176634e-19 if heating else np.nan,
          note='Sum of nonoverlapping filter bins/scores within each batch; flux scores are volume integrals, not pointwise flux'))
    pd.DataFrame(rows).to_csv(Path(out)/'integrated_tally_totals.csv',index=False)
    np.savez_compressed(Path(out)/'integrated_batch_scores.npz',values=batches,names=np.array([t.name for t in a.model.tallies]))



def write_heating(a,out,bins,batches):
    labels=[]
    for cell_id in bins:
        record=a.registry.records[int(cell_id)];role=record['component']
        if role=='blanket':role='first_wall_W' if record.get('kind')=='armour' else 'blanket_'+record.get('kind','other')
        labels.append(role)
    rows=[]
    for name in sorted(set(labels)):
        values=batches[:,np.array(labels)==name].sum(axis=1)
        mean=values.mean();se=values.std(ddof=1)/np.sqrt(len(values))
        rows.append(dict(component=name,eV_per_source=mean,standard_error=se,
            power_mw=mean*a.config.neutron_rate*1.602176634e-25,
            power_standard_error_mw=se*a.config.neutron_rate*1.602176634e-25))
    pd.DataFrame(rows).to_csv(Path(out)/'component_heating.csv',index=False)
    np.savez_compressed(Path(out)/'heating_batch_scores.npz',values=batches,labels=np.array(labels),cell_ids=np.asarray(bins))


def export_response_vtk(a,out):
    """Cartesian flux mesh with explicit bin-volume normalization."""
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk
    if not a.config.enable_mesh:return
    manifest=json.loads((Path(out)/'tally_manifest.json').read_text())
    record=next(r for r in manifest if r['name']=='Spatial neutron flux')
    scores=np.load(Path(out)/record['file'])
    low=np.array([-1800.,-1800.,-1200.]);high=-low;dimension=np.array(a.config.mesh_dimension)
    spacing=(high-low)/dimension;volume=float(np.prod(spacing))
    image=vtk.vtkImageData();image.SetDimensions(*(dimension+1));image.SetOrigin(*low);image.SetSpacing(*spacing)
    for name,values in [('flux_per_source_cm-2',scores['mean'].ravel()/volume),
       ('flux_se_per_source_cm-2',scores['standard_error'].ravel()/volume),
       ('flux_per_second_cm-2_s-1',scores['mean'].ravel()/volume*a.config.neutron_rate)]:
        array=numpy_to_vtk(values.copy(),deep=True);array.SetName(name);image.GetCellData().AddArray(array)
    writer=vtk.vtkXMLImageDataWriter();writer.SetFileName(str(Path(out)/'neutron_flux.vti'));writer.SetInputData(image);writer.Write()
