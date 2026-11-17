"""Discrete central-solenoid modules and separate casing/column shielding."""
import openmc


def annulus(inner,outer,zlo,zhi):
    region=-openmc.ZCylinder(r=outer)&+openmc.ZPlane(z0=zlo)&-openmc.ZPlane(z0=zhi)
    if inner>0:region &= +openmc.ZCylinder(r=inner)
    return region


def make_column(c,m,registry):
    cells=[];envelopes=[];dz=c.cs_height/c.cs_modules
    for j in range(c.cs_modules):
        lo=-c.cs_height/2+j*dz+c.cs_gap/2;hi=lo+dz-c.cs_gap
        envelope=annulus(c.cs_inner_radius,c.cs_outer_radius,lo,hi)
        envelopes.append(envelope)
        insulated=annulus(c.cs_inner_radius+c.cs_case,c.cs_outer_radius-c.cs_case,lo+c.cs_case,hi-c.cs_case)
        winding=annulus(c.cs_inner_radius+c.cs_case+c.cs_insulation,c.cs_outer_radius-c.cs_case-c.cs_insulation,
                        lo+c.cs_case+c.cs_insulation,hi-c.cs_case-c.cs_insulation)
        for name,region,material in [('CS winding',winding,m['wp']),('CS insulation',insulated&~winding,m['glass']),('CS casing',envelope&~insulated,m['ss'])]:
            cells.append(registry.cell(f'{name} {j}',region,material,component='central_solenoid',module=j))
    cells.append(registry.cell('Central support spine',annulus(0,c.support_spine_radius,-c.support_height/2,c.support_height/2),m['ss'],component='central_support'))
    cells.append(registry.cell('Central outer support',annulus(*c.support_radii,-c.support_height/2,c.support_height/2),m['ss'],component='central_support'))
    cells.append(registry.cell('Inboard column shielding',annulus(*c.column_shield_radii,-c.column_shield_height/2,c.column_shield_height/2),m['shield'],component='inboard_shield'))
    envelopes.extend(cell.region for cell in cells[-3:])
    return cells,envelopes
