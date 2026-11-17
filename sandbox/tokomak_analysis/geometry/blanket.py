"""Repeated toroidal modules with poloidal gaps, sidewalls, beds and channels."""
import numpy as np
import openmc
from .primitives import sector,poloidal_sector,union
from .breeder_units import radial_layers


def make_blanket(c,m,g,registry):
    lower=+openmc.ZPlane(z0=c.blanket_bottom)
    envelope=(-g['blanket_outer']&~(-g['fw']))&lower
    sector_width=2*np.pi/c.tf_coils
    module_width=sector_width/c.toroidal_modules_per_sector
    # One canonical toroidal module; repetitions preserve exact geometry.
    toroidal_core=sector(-module_width/2,module_width/2,c.toroidal_gap/2+c.module_sidewall)
    toroidal_active=sector(-module_width/2,module_width/2,c.toroidal_gap/2)
    module_universe=openmc.Universe(name='Poloidally segmented blanket module')
    layers=radial_layers(c,m,g);poloidal_regions=[]
    # Engineering poloidal segmentation: independent IB and OB vertical stacks.
    # Plane/cylinder boundaries avoid extremely shallow cones near the equator.
    stacks=c.poloidal_modules//2
    zmin=max(c.blanket_bottom,float(g['blanket_outer'].points[:,1].min()))
    zmax=float(g['blanket_outer'].points[:,1].max())
    if (zmax-zmin)/stacks<=c.poloidal_gap_cm+2*c.poloidal_sidewall_cm:raise ValueError('Poloidal gaps/walls consume the module height')
    edges=np.linspace(zmin,zmax,stacks+1)
    mid=openmc.ZCylinder(r=c.major_radius)
    def band(lo,hi):return +openmc.ZPlane(z0=lo)&-openmc.ZPlane(z0=hi)
    for j in range(c.poloidal_modules):
        row=j%stacks;inboard=j<stacks
        lo=float(edges[row]);hi=float(edges[row+1]);gap=c.poloidal_gap_cm/2;wall=c.poloidal_sidewall_cm
        side=-mid if inboard else +mid
        active=band(lo+gap,hi-gap)&side
        radial_core=-openmc.ZCylinder(r=c.major_radius-wall) if inboard else +openmc.ZCylinder(r=c.major_radius+wall)
        core=band(lo+gap+wall,hi-gap-wall)&radial_core
        poloidal_regions.append(active)
        cells=[]
        for layer in layers:
            region=-layer['outer'] & ~(-layer['inner'])
            kind=layer['kind'];material=layer['material']
            if c.fidelity==3 and kind in ('breeder','multiplier','cooling','manifold'):
                channel_regions=[]
                n=c.grouped_channels if kind in ('cooling','manifold') else c.ribs_per_module
                fraction=(c.fw_helium_fraction if layer['name'].startswith('FW') else c.cooling_helium_fraction) if kind=='cooling' else c.manifold_helium_fraction
                if kind in ('breeder','multiplier'):fraction=c.bed_structure_fraction
                start=lo+gap+wall;span=hi-gap-wall-start
                for k in range(n):
                    left=start+(k+.5)*span/n-fraction*span/n/2
                    channel_regions.append(band(left,left+fraction*span/n))
                strips=union(channel_regions)
                if kind in ('cooling','manifold'):fluid=m['he'];solid=m['steel']
                else:fluid=m['steel'];solid=m[f'breeder{layer["zone"]}'] if kind=='breeder' else m['multiplier']
                cells.append(registry.cell(layer['name']+' channels/ribs',region&strips,fluid,component='blanket',kind=kind,zone=layer['zone'],poloidal=j))
                cells.append(registry.cell(layer['name']+' matrix',region&~strips,solid,component='blanket',kind=kind,zone=layer['zone'],poloidal=j))
            else:cells.append(registry.cell(layer['name'],region,material,component='blanket',kind=kind,zone=layer['zone'],poloidal=j))
        radial_universe=openmc.Universe(cells=cells,name=f'Breeder unit {j}')
        parent=registry.cell(f'Poloidal module {j}',active&core&toroidal_core,radial_universe,component='blanket_module',poloidal=j,inboard=inboard)
        module_universe.add_cell(parent)
        module_universe.add_cell(registry.cell(f'Module sidewall {j}',active&~(core&toroidal_core),m['steel'],component='blanket',kind='sidewall',poloidal=j))
    gaps=[band(float(z)-c.poloidal_gap_cm/2,float(z)+c.poloidal_gap_cm/2) for z in edges]
    module_universe.add_cell(registry.cell('Poloidal intermodule gaps',union(gaps),m['he'],component='blanket_gap'))
    toroidal_universe=openmc.Universe(name='Toroidal blanket module array');active_modules=[]
    for k in range(c.toroidal_modules_per_sector):
        lo=-sector_width/2+k*module_width;hi=lo+module_width;centre=(lo+hi)/2
        region=sector(lo,hi,c.toroidal_gap/2);active_modules.append(region)
        cell=registry.cell(f'Toroidal module {k}',region,module_universe,component='toroidal_module',toroidal_local=k)
        cell.rotation=(0.,0.,np.rad2deg(centre))
        toroidal_universe.add_cell(cell)
    toroidal_universe.add_cell(registry.cell('Toroidal intermodule gaps',~union(active_modules),m['he'],component='blanket_gap'))
    return registry.cell('Blanket envelope',envelope,toroidal_universe,component='blanket_envelope'),layers
