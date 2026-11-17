"""Build, validate, transport, plot and export the native CSG tokamak."""
from pathlib import Path
from dataclasses import asdict,replace
import argparse,json,time,resource
import numpy as np
import openmc
import openmc.lib as lib
from .config import Config
from .geometry.assembly import build_model


def export_inputs(a,out):
    out.mkdir(parents=True,exist_ok=True);a.config.save(out/'config.json')
    a.model.geometry.merge_surfaces=False
    a.model.export_to_xml(directory=out)
    (out/'components.json').write_text(json.dumps(list(a.registry.records.values()),indent=2))
    a.model.geometry.determine_paths()
    lookup={}
    for cell in a.model.geometry.get_all_material_cells().values():
        lookup[str(cell.id)]=cell.paths
    (out/'cell_paths.json').write_text(json.dumps(lookup,indent=2))
    (out/'build_metrics.json').write_text(json.dumps(dict(cells=len(a.model.geometry.get_all_cells()),
      physical_leaf_instances=sum(len(cell.paths) for cell in a.model.geometry.get_all_cells().values() if cell.fill_type!='universe'),
      nominal_blanket_modules=a.config.tf_coils*a.config.toroidal_modules_per_sector*a.config.poloidal_modules,
      surfaces=len(a.model.geometry.get_all_surfaces()),materials=len(a.model.geometry.get_all_materials()),construction_seconds=a.construction_seconds),indent=2))


def run_transport(a,out,debug=False):
    from .neutronics.postprocess import write_outputs
    # The same response verifier used by the sensitivity suite verifies native
    # Li reaction channels against pinned evaluation metadata in a separate step.
    export_inputs(a,out);start=time.perf_counter();batches=[]
    a.model.init_lib(directory=str(out),threads=a.config.threads,geometry_debug=debug,output=True)
    try:
        from .neutronics.library import verify_library
        verify_library(a,out)
        first=next(t for t in a.model.tallies if t.name=='Li6 tritium production')
        second=next(t for t in a.model.tallies if t.name=='Li7 tritium channels')
        third=next(t for t in a.model.tallies if t.name=='Be9 tritium production')
        # Check loaded data: the partial sum must reproduce native baseline205.
        for name,mts in [('Li6',[105]),('Li7',range(52,83)),('Be9',[700,701])]:
            idx=lib.nuclides[name]._index;grid=lib.nuclide_energy_grid(idx)
            def full(mt):
                tail=lib.get_nuclide_xs(idx,mt);v=np.zeros_like(grid);v[-len(tail):]=tail;return v
            np.testing.assert_allclose(sum(full(mt) for mt in mts),full(205),rtol=1e-6,atol=1e-12)
        module6=next(t for t in a.model.tallies if t.name=='Module Li6 production')
        module7=next(t for t in a.model.tallies if t.name=='Module Li7 production')
        module_be=next(t for t in a.model.tallies if t.name=='Module Be9 production')
        module_bins=module6.find_filter(openmc.CellInstanceFilter).bins
        module_batches=[];module_previous=np.zeros(len(module_bins))
        integrated_batches=[];integrated_previous=np.zeros(len(a.model.tallies))
        heat=next((t for t in a.model.tallies if t.name=='Component nuclear heating'),None)
        heat_batches=[];heat_previous=np.zeros(len(heat.find_filter(openmc.CellFilter).bins)) if heat else None
        lib.simulation_init();previous=np.zeros(3)
        for _ in range(a.config.batches):
            lib.next_batch()
            # Native SUM is the per-batch accumulated normalized score.
            cumulative=np.array([lib.tallies[first.id].results[:,:,1].sum(),lib.tallies[second.id].results[:,:,1].sum(),lib.tallies[third.id].results[:,:,1].sum()])
            batches.append(cumulative-previous);previous=cumulative
            module_cumulative=(lib.tallies[module6.id].results[:,:,1].reshape(len(module_bins),-1).sum(axis=1)
                               +lib.tallies[module7.id].results[:,:,1].reshape(len(module_bins),-1).sum(axis=1)
                               +lib.tallies[module_be.id].results[:,:,1].reshape(len(module_bins),-1).sum(axis=1))
            module_batches.append(module_cumulative-module_previous);module_previous=module_cumulative
            all_cumulative=np.array([lib.tallies[t.id].results[:,:,1].sum() for t in a.model.tallies])
            integrated_batches.append(all_cumulative-integrated_previous);integrated_previous=all_cumulative
            if heat:
                h=lib.tallies[heat.id].results[:,:,1].sum(axis=1)
                heat_batches.append(h-heat_previous);heat_previous=h
        lib.simulation_finalize()
        lib.statepoint_write(str((out/f'statepoint.{a.config.batches}.h5').resolve()))
        from .neutronics.postprocess import export_native_tallies,write_module_results
        export_native_tallies(a,out)
        write_module_results(a,out,module_bins,np.asarray(module_batches),np.asarray(batches))
        if heat:
            from .neutronics.postprocess import write_heating
            write_heating(a,out,heat.find_filter(openmc.CellFilter).bins,np.asarray(heat_batches))
        from .neutronics.postprocess import export_response_vtk
        export_response_vtk(a,out)
        from .validation import sample_geometry
        sample_geometry(a,out,count=2000,native=True)
    finally:a.model.finalize_lib()
    elapsed=time.perf_counter()-start
    from .neutronics.postprocess import write_integrated_results
    write_integrated_results(a,out,np.asarray(integrated_batches))
    result=write_outputs(a,out,np.asarray(batches),dict(histories=a.config.particles*a.config.batches,
        elapsed_seconds=elapsed,histories_per_second=a.config.particles*a.config.batches/elapsed,
        peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,geometry_debug=debug))
    print(json.dumps(result,indent=2));return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path);parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'results/segmented')
    parser.add_argument('--fidelity',type=int,choices=[1,2,3]);parser.add_argument('--particles',type=int);parser.add_argument('--batches',type=int)
    parser.add_argument('--run',action='store_true');parser.add_argument('--debug',action='store_true');parser.add_argument('--plot',action='store_true')
    parser.add_argument('--vtk',action='store_true');parser.add_argument('--validate',action='store_true')
    args=parser.parse_args();c=Config.load(args.config) if args.config else Config()
    overrides={k:getattr(args,k) for k in ('fidelity','particles','batches') if getattr(args,k) is not None};c=replace(c,**overrides).validate()
    a=build_model(c,output=args.output);export_inputs(a,args.output)
    if args.validate:
        from .validation import sample_geometry
        sample_geometry(a,args.output,count=20000)
    if args.run:run_transport(a,args.output,args.debug)
    if args.plot or args.vtk:
        from .visualization.export_vtk import visualize
        visualize(a,args.output,plots=args.plot,vtk=args.vtk)


if __name__=='__main__':main()
