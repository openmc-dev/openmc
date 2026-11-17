"""Current native partial-reaction scoring; no cached-MT205 primary response."""
import numpy as np
import openmc


def tally(name,filters,scores,nuclides=None,estimator='tracklength'):
    t=openmc.Tally(name=name);t.filters=filters;t.scores=scores;t.estimator=estimator
    if nuclides:t.nuclides=nuclides
    return t


def add_tallies(a):
    c=a.config;records=a.registry.records;cells=a.registry.cells
    blanket=[i for i,r in records.items() if r.get('component')=='blanket_envelope']
    modules=[i for i,r in records.items() if r.get('component')=='blanket_module']
    a.model.geometry.determine_paths()
    instances=[(cell,index) for cell in modules for index in range(len(cells[cell].paths))]
    # Parent-cell filters match native coordinate levels; nuclide scoring uses
    # the actual material leaf and its atom densities.
    total6=tally('Li6 tritium production',[],['105'],['Li6'])
    total7=tally('Li7 tritium channels',[],list(map(str,range(52,83))),['Li7'])
    # All triton-producing partial scores are reduced per batch in the runner.
    total_be=tally('Be9 tritium production',[],['700','701'],['Be9'])
    tallies=[total6,total7,total_be,
      tally('Module Li6 production',[openmc.CellInstanceFilter(instances)],['105'],['Li6']),
      tally('Module Li7 production',[openmc.CellInstanceFilter(instances)],list(map(str,range(52,83))),['Li7']),
      tally('Module Be9 production',[openmc.CellInstanceFilter(instances)],['700','701'],['Be9']),
      tally('Whole model flux',[],['flux']),tally('Neutron multiplication n2n',[],['(n,2n)'])]
    energy=np.geomspace(1e-5,20e6,81)
    tallies.append(tally('Blanket flux spectrum',[openmc.CellFilter(blanket),openmc.EnergyFilter(energy)],['flux']))
    world=next(s for s in a.model.geometry.get_all_surfaces().values() if s.boundary_type=='vacuum')
    column_surface=next(s for s in a.model.geometry.get_all_surfaces().values() if isinstance(s,openmc.ZCylinder) and s.r==c.column_shield_radii[1])
    tallies.append(tally('Central shield outer-surface neutron current',[openmc.SurfaceFilter(column_surface)],['current'],estimator='analog'))
    tallies.append(tally('World neutron leakage',[openmc.SurfaceFilter(world)],['current'],estimator='analog'))
    selected=[i for i,r in records.items() if r.get('component') in ('vessel','TF_winding','PF_winding','central_solenoid')]
    if c.enable_exposure:
        tallies.append(tally('Vessel and magnet flux spectrum',[openmc.CellFilter(selected),openmc.EnergyFilter(energy)],['flux']))
        tallies.append(tally('Fast magnet flux integral above 0.1 MeV',[openmc.CellFilter(selected),openmc.EnergyFilter([1e5,20e6])],['flux']))
    if c.enable_heating:
        targets=[i for i,r in records.items() if r.get('component') in ('blanket','divertor','TF_winding','PF_winding','vessel','central_solenoid','central_support')]
        # Heating-local is neutron-local photon-deposition approximation.
        score='heating' if c.coupled_photons else 'heating-local'
        tallies.append(tally('Component nuclear heating',[openmc.CellFilter(targets)],[score],estimator='collision' if c.coupled_photons else 'tracklength'))
    if c.coupled_photons:
        for t in tallies:t.filters=[openmc.ParticleFilter('neutron'),*t.filters]
        # Include energy deposition by BOTH species for coupled heating.
        if c.enable_heating:
            t=next(t for t in tallies if t.name=='Component nuclear heating')
            t.filters=[f for f in t.filters if not isinstance(f,openmc.ParticleFilter)]
    if c.enable_mesh:
        mesh=openmc.RegularMesh(name='Reactor Cartesian response mesh');mesh.dimension=c.mesh_dimension
        mesh.lower_left=(-1800,-1800,-1200);mesh.upper_right=(1800,1800,1200)
        tallies.append(tally('Spatial neutron flux',[openmc.MeshFilter(mesh),openmc.ParticleFilter('neutron')],['flux']))
    a.model.tallies=openmc.Tallies(tallies)
