"""Identical radial inventory definitions for three blanket fidelity levels."""
import numpy as np
from ..materials import mix
from .primitives import polygon_volume


def radial_layers(c,m,g):
    layers=[];depth=c.plasma_wall_gap
    def add(name,thickness,material,kind,zone=-1):
        nonlocal depth
        inner=g['offset'](depth);depth=depth+thickness;outer=g['offset'](depth)
        layers.append(dict(name=name,inner=inner,outer=outer,material=material,kind=kind,zone=zone))
    for name,size,key,kind in [('W armour',c.armour,'w','armour'),('FW front steel',c.fw_front_steel,'steel','steel'),
      ('FW helium cooling',c.fw_cooling,'fw_cool','cooling'),('FW rear steel',c.fw_rear_steel,'steel','steel')]:add(name,size,m[key],kind)
    active_depth=g['build']-c.fw_total-c.backplate-c.manifold
    for j in range(c.breeder_cycles):
        for label,fraction,key,kind in [('Multiplier',c.multiplier_fraction,'multiplier_struct','multiplier'),
          ('Li2O breeder',c.breeder_fraction,f'breeder_struct{j}','breeder'),
          ('Cooling plate',c.cooling_fraction,'cool','cooling')]:
            add(f'{label} {j}',active_depth*fraction/c.breeder_cycles,m[key],kind,j)
    add('Backplate',c.backplate,m['steel'],'steel')
    add('Inlet/outlet manifold',c.manifold,m['manifold'],'manifold')
    if c.fidelity==1:
        active=[p for p in layers if p['kind'] in ('breeder','multiplier') or p['name'].startswith('Cooling plate')]
        volumes=np.array([polygon_volume(p['outer'].points)-polygon_volume(p['inner'].points) for p in active])
        homogeneous=mix('Global blanket active homogenization',[p['material'] for p in active],(volumes/volumes.sum()).tolist(),c.breeder_temperature)
        middle=dict(name='Global active breeding zone',inner=active[0]['inner'],outer=active[-1]['outer'],material=homogeneous,kind='breeder',zone=-1)
        layers=layers[:4]+[middle]+layers[-2:]
    return layers
