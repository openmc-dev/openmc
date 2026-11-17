"""Common-random-point native CSG inventory estimates, with sampling errors."""
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
import pandas as pd
import openmc.lib as lib
from .config import Config
from .geometry.assembly import build_model
from .geometry.primitives import polygon_volume


def compare(out,count=100000):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    lower=np.array([-1850.,-1850.,-1200.]);upper=-lower;volume=float(np.prod(upper-lower))
    points=np.random.default_rng(97531).uniform(lower,upper,(count,3))
    samples=[];rows=[]
    for level in (1,2,3):
        a=build_model(replace(Config(),fidelity=level,enable_mesh=False),tallies=False,output=out/f'level{level}')
        materials=a.model.geometry.get_all_materials();lookup={}
        for material in materials.values():
            atoms=material.get_nuclide_atom_densities()
            lookup[material.id]=np.array([atoms.get('Li6',0)*1e24,atoms.get('Li7',0)*1e24,atoms.get('Be9',0)*1e24,material.get_mass_density()])
        values=np.zeros((count,4));plasma=np.zeros(count)
        a.model.init_lib(directory=str(out/f'level{level}'),threads=2,output=False)
        try:
            for j,point in enumerate(points):
                if np.linalg.norm(point)>=a.config.world_radius:continue
                cell,instance=lib.find_cell(point);record=a.registry.records[cell.id]
                if record['material_id']:values[j]=lookup[record['material_id']]
                plasma[j]=record['component']=='plasma'
        finally:a.model.finalize_lib()
        samples.append(values)
        for k,name in enumerate(('Li6_atoms','Li7_atoms','Be9_atoms','total_mass_g')):
            rows.append(dict(fidelity=level,quantity=name,estimate=float(values[:,k].mean()*volume),
                standard_error=float(values[:,k].std(ddof=1)/np.sqrt(count)*volume)))
        estimated=plasma.mean()*volume;exact=polygon_volume(a.contours['plasma'].points)
        if abs(estimated-exact)>5*plasma.std(ddof=1)/np.sqrt(count)*volume:raise ValueError('Plasma inventory estimate inconsistent with exact revolved-polygon volume')
        print('Inventory complete',level,flush=True)
    pd.DataFrame(rows).to_csv(out/'inventories.csv',index=False)
    differences=[]
    for level in (1,3):
        delta=samples[level-1]-samples[1]
        for k,name in enumerate(('Li6_atoms','Li7_atoms','Be9_atoms','total_mass_g')):
            reference=samples[1][:,k].mean()
            differences.append(dict(fidelity=level,quantity=name,relative_change=float(delta[:,k].mean()/reference),
                relative_paired_standard_error=float(delta[:,k].std(ddof=1)/np.sqrt(count)/reference)))
    pd.DataFrame(differences).to_csv(out/'fidelity_inventory_comparison.csv',index=False)
    np.savez_compressed(out/'common_point_inventory_samples.npz',values=np.array(samples),points=points)
    (out/'method.json').write_text(json.dumps(dict(samples=count,bounds_cm=[lower.tolist(),upper.tolist()],
        sampling='uniform physical volume, identical points for all levels, native find_cell',
        limitation='Global active-zone homogenization uses whole-ring weights; bottom truncation and module gaps can change exact inventories. Differences and paired sampling errors are reported rather than assumed zero.'),indent=2))


if __name__=='__main__':compare(Path(__file__).resolve().parent/'results/inventory')
