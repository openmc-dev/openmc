"""Hierarchy-aware point overlap checks and authoritative component queries."""
from pathlib import Path
import json
import numpy as np
import openmc


def walk(universe,point):
    p=np.array(point,float)
    matching=[cell for cell in universe.cells.values() if p in cell]
    if len(matching)!=1:
        raise ValueError(f'Geometry coverage: universe {universe.id}, point {p.tolist()}, matches {[(c.id,c.name) for c in matching]}')
    cell=matching[0]
    if cell.fill_type=='universe':
        if cell.translation is not None:p-=cell.translation
        if cell.rotation is not None:p=cell.rotation_matrix@p
        return [cell]+walk(cell.fill,p)
    return [cell]


def metadata_for_path(a,path):
    result=dict(a.registry.records[path[-1].id]);result['sector']=-1;result['module_id']=-1
    sector=local=poloidal=-1
    for cell in path:
        r=a.registry.records[cell.id]
        if 'sector' in r:sector=r['sector']
        if 'toroidal_local' in r:local=r['toroidal_local']
        if r.get('component')=='blanket_module':poloidal=r['poloidal']
    result['sector']=sector
    if local>=0 and poloidal>=0:
        result['module_id']=(sector*a.config.toroidal_modules_per_sector+local)*a.config.poloidal_modules+poloidal
    result['poloidal']=poloidal
    result['inboard']=-1
    return result


def sample_geometry(a,out,count=10000,native=False):
    rng=np.random.default_rng(a.config.seed+501)
    points=rng.uniform([-1800,-1800,-1200],[1800,1800,1200],(count,3))
    # Exclude positions outside the spherical world; the native boundary ends there.
    points=points[np.linalg.norm(points,axis=1)<a.config.world_radius]
    counts={};failures=[]
    for point in points:
        try:
            if native:
                import openmc.lib as lib
                cell,instance=lib.find_cell(point);key=str(cell.id)
            else:
                path=walk(a.model.geometry.root_universe,point);key=str(path[-1].id)
            counts[key]=counts.get(key,0)+1
        except Exception as error:
            failures.append(dict(point=point.tolist(),error=str(error)))
            if len(failures)>=10:break
    report=dict(points=len(points),native=native,failures=failures,occupied_cells=len(counts),counts=counts)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/('native_point_validation.json' if native else 'point_validation.json')).write_text(json.dumps(report,indent=2))
    if failures:raise ValueError(f'Geometry validation failed: {failures[0]}')
    return report
