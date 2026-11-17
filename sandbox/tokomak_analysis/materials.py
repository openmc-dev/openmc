"""Documented elemental mass fractions and volume-consistent bed mixing."""
import numpy as np
import openmc


def alloy(name,rho,composition,temp):
    m=openmc.Material(name=name);m.set_density('g/cm3',rho);m.temperature=temp
    for element,fraction in composition.items():m.add_element(element,fraction,percent_type='wo')
    return m


def mix(name,materials,fractions,temp):
    if not np.isclose(sum(fractions),1.):raise ValueError('Volume fractions must sum to one')
    m=openmc.Material.mix_materials(materials,fractions,'vo',name=name);m.temperature=temp
    return m


def make_materials(c):
    steel=alloy('EUROFER97 representative heat',7.8,
      dict(Fe=.89024,Cr=.0893,C=.0012,Mn=.0047,V=.002,W=.0107,Ta=.0014,Ti=.00009,N=.00018,Si=.00006,S=.00004,Ni=.00002,Mo=.000015,Nb=.00002,Cu=.00003),c.structure_temperature)
    ss=alloy('SS316L engineering nominal',8.,dict(Fe=.655,Cr=.17,Ni=.12,Mo=.025,Mn=.02,Si=.0097,C=.0003),c.structure_temperature)
    w=alloy('Tungsten',19.25,{'W':1.},c.structure_temperature)
    be=alloy('Beryllium solid',1.85,{'Be':1.},c.breeder_temperature)
    cu=alloy('Copper',8.96,{'Cu':1.},c.structure_temperature)
    cucrzr=alloy('CuCrZr nominal',8.9,dict(Cu=.99,Cr=.008,Zr=.002),c.structure_temperature)
    he=alloy('Helium coolant ideal gas',c.helium_pressure_mpa*1e6*.004002602/(8.314462618*c.coolant_temperature)/1000,{'He':1.},c.coolant_temperature)
    purge=alloy('Helium purge ideal gas',c.purge_pressure_mpa*1e6*.004002602/(8.314462618*c.breeder_temperature)/1000,{'He':1.},c.breeder_temperature)
    cold=alloy('Helium magnet coolant equivalent',.125,{'He':1.},c.transport_magnet_temperature)
    glass=openmc.Material(name='Silica glass insulation approximation');glass.add_element('Si',1);glass.add_element('O',2);glass.set_density('g/cm3',2.2);glass.temperature=c.transport_magnet_temperature
    nbsn=openmc.Material(name='Nb3Sn');nbsn.add_element('Nb',3);nbsn.add_element('Sn',1);nbsn.set_density('g/cm3',8.4);nbsn.temperature=c.transport_magnet_temperature
    water=openmc.Material(name='Pressurized water equivalent');water.add_nuclide('H1',2);water.add_nuclide('O16',1);water.set_density('g/cm3',.7);water.temperature=c.structure_temperature
    # S(a,b) is intentionally not applied to pressurized hot-water approximation;
    # this limitation is reported. The model is fast-spectrum fusion neutronics.
    wp=mix('Nb3Sn winding pack homogenization',[ss,cu,nbsn,glass,cold],[.35,.35,.1,.1,.1],c.transport_magnet_temperature)
    if c.multiplier=='Be12Ti':
        multiplier=openmc.Material(name='Be12Ti solid');multiplier.add_element('Be',12);multiplier.add_element('Ti',1);multiplier.set_density('g/cm3',2.28);multiplier.temperature=c.breeder_temperature
    else:multiplier=be
    multbed=mix('Multiplier packed bed',[multiplier,purge],[c.multiplier_packing,1-c.multiplier_packing],c.breeder_temperature)
    result=dict(steel=steel,ss=ss,w=w,be=be,cu=cu,cucrzr=cucrzr,he=he,purge=purge,glass=glass,wp=wp,water=water,
      fw_cool=mix('FW steel helium',[steel,he],[1-c.fw_helium_fraction,c.fw_helium_fraction],c.structure_temperature),
      cool=mix('Cooling plate steel helium',[steel,he],[1-c.cooling_helium_fraction,c.cooling_helium_fraction],c.structure_temperature),
      manifold=mix('Blanket manifold steel helium',[steel,he],[1-c.manifold_helium_fraction,c.manifold_helium_fraction],c.structure_temperature),
      vessel_mix=mix('VV shielding steel water',[ss,water],[c.vessel_steel_fraction,1-c.vessel_steel_fraction],c.structure_temperature),
      shield=mix('Central column steel tungsten',[ss,w],[.8,.2],c.structure_temperature),multiplier=multbed)
    result['multiplier_struct']=mix('Multiplier bed with equivalent ribs',[multbed,steel],[1-c.bed_structure_fraction,c.bed_structure_fraction],c.breeder_temperature)
    for i,enrichment in enumerate(c.enrichment):
        solid=openmc.Material(name=f'Li2O solid zone {i} e6={enrichment:g}')
        if enrichment>0:solid.add_nuclide('Li6',2*enrichment)
        if enrichment<1:solid.add_nuclide('Li7',2*(1-enrichment))
        solid.add_element('O',1);solid.set_density('g/cm3',c.li2o_solid_density);solid.temperature=c.breeder_temperature
        bed=mix(f'Li2O packed bed zone {i}',[solid,purge],[c.breeder_packing,1-c.breeder_packing],c.breeder_temperature)
        result[f'breeder{i}']=bed
        result[f'breeder_struct{i}']=mix(f'Li2O bed equivalent ribs {i}',[bed,steel],[1-c.bed_structure_fraction,c.bed_structure_fraction],c.breeder_temperature)
    return result
