"""Plasma-shaped isotropic DT sources with explicit power normalization."""
from pathlib import Path
import numpy as np
import openmc


def radial_coordinate(points,boundary,r0):
    p=boundary-np.array([r0,0.]);q=np.roll(p,-1,axis=0)
    edge=q-p;normals=np.column_stack([edge[:,1],-edge[:,0]])
    distances=np.einsum('ij,ij->i',normals,p)
    return np.max((points-np.array([r0,0.]))@normals.T/distances,axis=1)


def sample_profile(c,boundary,count,seed):
    """Volume-uniform proposals in cylindrical coordinates, profile rejection.

    Used only for the optional asymmetric source-bank mode, not for transport
    collisions or secondary nuclear-data sampling.
    """
    rng=np.random.default_rng(seed);accepted=[];n=0
    rlo,rhi=boundary[:,0].min(),boundary[:,0].max();zlo,zhi=boundary[:,1].min(),boundary[:,1].max()
    while n<count:
        size=max(1024,2*(count-n));r=np.sqrt(rng.uniform(rlo*rlo,rhi*rhi,size));z=rng.uniform(zlo,zhi,size)
        rz=np.column_stack([r,z]);rho=radial_coordinate(rz,boundary,c.major_radius)
        theta=np.arctan2(z/c.elongation,r-c.major_radius)
        probability=np.maximum(0,1-rho*rho)**c.radial_profile_power*(1+c.poloidal_asymmetry*np.cos(theta))/(1+abs(c.poloidal_asymmetry))
        keep=(rho<=1)&(rho>=0)&(rng.random(size)<probability)
        rz=rz[keep][:count-n];phi=rng.uniform(0,2*np.pi,len(rz))
        accepted.append(np.column_stack([rz[:,0]*np.cos(phi),rz[:,0]*np.sin(phi),rz[:,1]]));n+=len(rz)
    return np.concatenate(accepted)


def make_source(c,plasma,boundary,output=None):
    energy=openmc.stats.Discrete([c.source_energy_ev],[1.])
    if c.source_energy_sigma_ev>0:energy=openmc.stats.Normal(c.source_energy_ev,c.source_energy_sigma_ev)
    if c.poloidal_asymmetry:
        if output is None:raise ValueError('An output directory is required for a nonuniform poloidal source bank')
        directory=Path(output);directory.mkdir(parents=True,exist_ok=True)
        points=sample_profile(c,boundary,c.source_bank_size,c.seed+81)
        rng=np.random.default_rng(c.seed+82);mu=rng.uniform(-1,1,len(points));phi=rng.uniform(0,2*np.pi,len(points))
        u=np.column_stack([np.sqrt(1-mu*mu)*np.cos(phi),np.sqrt(1-mu*mu)*np.sin(phi),mu])
        energies=np.full(len(points),c.source_energy_ev)
        if c.source_energy_sigma_ev:energies=rng.normal(c.source_energy_ev,c.source_energy_sigma_ev,len(points))
        if np.any(energies<=0):raise ValueError('Source energy width generated nonphysical energies')
        path=(directory/'plasma_source.h5').resolve()
        openmc.write_source_file([openmc.SourceParticle(r=p,u=v,E=e) for p,v,e in zip(points,u,energies)],path)
        return openmc.FileSource(path)
    constraints={'domains':[plasma],'rejection_strategy':'resample'}
    if c.radial_profile_power:
        grid=np.linspace(0,1,51)
        return openmc.TokamakSource(major_radius=c.major_radius,minor_radius=c.minor_radius,
            elongation=c.elongation,triangularity=c.triangularity,shafranov_shift=0.,
            r_over_a=grid,emission_density=(1-grid*grid)**c.radial_profile_power,energy=energy,constraints=constraints)
    # R-squared sampling ensures uniform physical volume, not uniform R density.
    rlo,rhi=boundary[:,0].min(),boundary[:,0].max()
    space=openmc.stats.CylindricalIndependent(r=openmc.stats.PowerLaw(rlo,rhi,1),
        phi=openmc.stats.Uniform(0,2*np.pi),z=openmc.stats.Uniform(boundary[:,1].min(),boundary[:,1].max()))
    return openmc.IndependentSource(space=space,angle=openmc.stats.Isotropic(),energy=energy,constraints=constraints)
