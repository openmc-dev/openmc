"""Sample authoritative native CSG cells; no independently defined CAD geometry."""
from pathlib import Path
import json,re
import numpy as np
import openmc
import openmc.lib as lib

COLORS={'plasma':(245,145,170),'blanket':(90,170,100),'blanket_gap':(190,230,240),
 'vessel':(125,135,150),'TF_winding':(145,75,180),'PF_winding':(170,100,190),
 'TF_case':(90,95,110),'TF_insulation':(190,160,210),'PF_case':(95,100,120),
 'divertor':(220,170,55),'central_solenoid':(175,95,190),'central_support':(120,120,135),
 'inboard_shield':(100,140,160),'chamber_vacuum':(250,250,250),'external_vacuum':(245,245,245)}


def path_metadata(a):
    a.model.geometry.determine_paths();result={}
    components=sorted({r['component'] for r in a.registry.records.values()})
    codes={name:i+1 for i,name in enumerate(components)}
    for cell in a.model.geometry.get_all_cells().values():
        if cell.fill_type=='universe':continue
        for instance,path in enumerate(cell.paths):
            r=dict(a.registry.records[cell.id]);sector=local=pol=-1
            for token in re.findall(r'c(\d+)',path):
                data=a.registry.records[int(token)]
                sector=data.get('sector',sector);local=data.get('toroidal_local',local)
                if data['component']=='blanket_module':pol=data['poloidal']
            module=(sector*a.config.toroidal_modules_per_sector+local)*a.config.poloidal_modules+pol if local>=0 and pol>=0 else -1
            result[(cell.id,instance)]=(cell.id,r['material_id'],codes[r['component']],module,sector)
    return result,codes


def write_vtk(a,out):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    lookup,codes=path_metadata(a);nx,ny,nz=a.config.vtk_dimension
    extent=np.array(a.config.vtk_extent);spacing=2*extent/[nx,ny,nz]
    origin=np.asarray(a.config.vtk_origin)
    axes=[np.linspace(o-e+d/2,o+e-d/2,n) for o,e,d,n in zip(origin,extent,spacing,[nx,ny,nz])]
    points=np.column_stack([v.ravel(order='F') for v in np.meshgrid(*axes,indexing='ij')])
    values=np.full((len(points),6),-1,dtype=np.int32)
    for j,point in enumerate(points):
        if np.linalg.norm(point)>=a.config.world_radius:continue
        cell,instance=lib.find_cell(point)
        if (cell.id,instance) not in lookup:raise ValueError(f'Unmapped native cell instance {cell.id}:{instance}')
        values[j,:5]=lookup[cell.id,instance]
        values[j,5]=0 if np.hypot(point[0],point[1])<a.config.major_radius else 1
        if j%100000==0:print(f'VTK native samples {j}/{len(points)}',flush=True)
    image=vtk.vtkImageData();image.SetDimensions(nx+1,ny+1,nz+1);image.SetOrigin(*(origin-extent));image.SetSpacing(*spacing)
    for i,name in enumerate(['cell_id','material_id','component_id','module_id','sector_id','inboard_outboard']):
        array=numpy_to_vtk(values[:,i].copy(),deep=True);array.SetName(name);image.GetCellData().AddArray(array)
    image.GetCellData().SetActiveScalars('component_id')
    writer=vtk.vtkXMLImageDataWriter();writer.SetFileName(str(out/'reactor.vti'));writer.SetInputData(image);writer.Write()
    blocks=vtk.vtkMultiBlockDataSet();index=0;inventory=[]
    for name,code in codes.items():
        count=int((values[:,2]==code).sum());inventory.append(dict(component=name,component_id=code,sampled_voxels=count))
        if not count:continue
        threshold=vtk.vtkThreshold();threshold.SetInputData(image)
        threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,'component_id')
        threshold.SetLowerThreshold(code);threshold.SetUpperThreshold(code);threshold.Update()
        block=vtk.vtkUnstructuredGrid();block.DeepCopy(threshold.GetOutput())
        blocks.SetBlock(index,block);blocks.GetMetaData(index).Set(vtk.vtkCompositeDataSet.NAME(),name)
        # Optional extracted voxel surfaces are useful in ParaView without
        # shipping all interior cells. These are stair-step surfaces, not CAD.
        surface=vtk.vtkDataSetSurfaceFilter();surface.SetInputData(block);surface.Update()
        writer=vtk.vtkXMLPolyDataWriter();writer.SetFileName(str(out/(name+'.vtp')));writer.SetInputData(surface.GetOutput());writer.Write()
        index+=1
    writer=vtk.vtkXMLMultiBlockDataWriter();writer.SetFileName(str(out/'reactor.vtm'));writer.SetInputData(blocks);writer.Write()
    reader=vtk.vtkXMLMultiBlockDataReader();reader.SetFileName(str(out/'reactor.vtm'));reader.Update()
    if reader.GetOutput().GetNumberOfBlocks()!=index:raise ValueError('VTK multiblock readback failed')
    report=dict(dimensions=[nx,ny,nz],spacing_cm=spacing.tolist(),extent_cm=extent.tolist(),origin_cm=origin.tolist(),
      blocks=index,vtk_readback_passed=True,classification='native openmc.lib.find_cell, mapped through authoritative CSG instance paths',
      limitation='Cell-centre voxel sampling; sub-voxel layers may vanish. Not an exact boundary mesh or transport geometry.',components=inventory)
    materials=[]
    for material in a.model.geometry.get_all_materials().values():
        name=material.name
        rgb=(130,140,155)
        if 'Li2O' in name:rgb=(75,165,225)
        elif 'Multiplier' in name or 'Beryllium' in name:rgb=(85,185,95)
        elif name=='Tungsten':rgb=(240,195,50)
        elif 'winding' in name:rgb=(145,75,180)
        elif 'insulation' in name:rgb=(190,160,210)
        elif 'water' in name.lower():rgb=(55,175,220)
        elif 'Helium' in name:rgb=(160,225,240)
        materials.append(dict(material_id=material.id,name=name,rgb=rgb,density_g_cm3=material.get_mass_density()))
    (out/'material_colors.json').write_text(json.dumps(materials,indent=2))
    (out/'export_metadata.json').write_text(json.dumps(report,indent=2))
    np.savez_compressed(out/'sampled_classification.npz',values=values,points=points)
    return report


def make_plots(a,out):
    plots=[]
    configurations=[('xy_midplane','xy',(0,0,0),(3800,3800)),('xz_full','xz',(0,0,0),(3800,2400)),
      ('yz_full','yz',(0,0,0),(3800,2400)),('blanket_outboard_detail','xz',(1240,0,150),(240,350)),
      ('blanket_inboard_detail','xz',(540,0,100),(200,350)),('divertor_detail','xz',(830,0,-510),(450,240))]
    for name,basis,origin,width in configurations:
        plot=openmc.Plot(name=name);plot.filename=str((out/(name+'.png')).resolve());plot.basis=basis;plot.origin=origin
        plot.width=width;plot.pixels=(900,700);plot.color_by='cell';plot.show_overlaps=True
        plot.colors={cell:COLORS.get(a.registry.records[cell.id]['component'],(80,80,80)) for cell in a.model.geometry.get_all_cells().values() if cell.fill_type!='universe'}
        # Material-specific colors expose breeder, multiplier and W layers.
        for cell in a.model.geometry.get_all_cells().values():
            if cell.fill_type!='material':continue
            name_m=cell.fill.name
            if 'Li2O' in name_m:plot.colors[cell]=(75,165,225)
            elif 'Multiplier' in name_m or 'Beryllium' in name_m:plot.colors[cell]=(85,185,95)
            elif name_m=='Tungsten':plot.colors[cell]=(240,195,50)
        plots.append(plot)
    a.model.plots=openmc.Plots(plots)


def sampled_plots(a,out):
    """Slice diagnostics sampled from native find_cell, independent of plotter."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    lookup,codes=path_metadata(a)
    palette=[(1.,1.,1.)];cell_colors={}
    for cell in a.model.geometry.get_all_cells().values():
        if cell.fill_type=='universe':continue
        role=a.registry.records[cell.id]['component'];rgb=COLORS.get(role,(80,80,80))
        if role=='blanket' and cell.fill_type=='material':
            material=cell.fill.name
            if 'Li2O' in material:rgb=(75,165,225)
            elif 'Multiplier' in material:rgb=(85,185,95)
            elif material=='Tungsten':rgb=(240,195,50)
            elif 'Helium' in material:rgb=(160,225,240)
            else:rgb=(130,140,155)
        colour=tuple(v/255 for v in rgb)
        if colour not in palette:palette.append(colour)
        cell_colors[cell.id]=palette.index(colour)
    colors=palette
    views=[('xy_midplane','xy',(0,0,0),(3800,3800)),('xz_full','xz',(0,0,0),(3800,2400)),
      ('yz_full','yz',(0,0,0),(3800,2400)),('blanket_outboard_detail','xz',(1240,0,150),(240,350)),
      ('blanket_inboard_detail','xz',(540,0,100),(200,350)),('divertor_detail','xz',(830,0,-510),(450,240))]
    for name,basis,origin,width in views:
        nx,ny=360,280
        u=np.linspace(-width[0]/2,width[0]/2,nx)+origin['xyz'.index(basis[0])]
        v=np.linspace(-width[1]/2,width[1]/2,ny)+origin['xyz'.index(basis[1])]
        values=np.zeros((ny,nx),int)
        for j,y in enumerate(v):
            for i,x in enumerate(u):
                point=np.array(origin,float);point['xyz'.index(basis[0])]=x;point['xyz'.index(basis[1])]=y
                if np.linalg.norm(point)>=a.config.world_radius:continue
                cell,instance=lib.find_cell(point);values[j,i]=cell_colors[cell.id]
        fig,ax=plt.subplots(figsize=(10,7))
        ax.imshow(values,origin='lower',extent=(u[0],u[-1],v[0],v[-1]),cmap=ListedColormap(colors),vmin=0,vmax=len(colors)-1,interpolation='nearest')
        ax.set(xlabel=basis[0]+' [cm]',ylabel=basis[1]+' [cm]',title=name+' — native CSG cell sampling')
        fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=180);plt.close(fig)
        print('Native diagnostic slice',name,flush=True)
    # Additional toroidal slices at 7.5 and 11.25 degrees query the same CSG.
    for angle in (7.5,11.25):
        r=np.linspace(200,1800,320);z=np.linspace(-1150,1150,320);values=np.zeros((len(z),len(r)),int)
        phi=np.deg2rad(angle)
        for j,zz in enumerate(z):
            for i,rr in enumerate(r):
                cell,instance=lib.find_cell((rr*np.cos(phi),rr*np.sin(phi),zz));values[j,i]=cell_colors[cell.id]
        fig,ax=plt.subplots(figsize=(7,9));ax.imshow(values,origin='lower',extent=(r[0],r[-1],z[0],z[-1]),cmap=ListedColormap(colors),vmin=0,vmax=len(colors)-1,interpolation='nearest')
        ax.set(xlabel='R [cm]',ylabel='Z [cm]',title=f'Native toroidal cut at phi={angle:g}°')
        fig.tight_layout();fig.savefig(out/(f'toroidal_{angle:g}.png'),dpi=180);plt.close(fig)


def visualize(a,out,plots=True,vtk=True):
    out=Path(out).resolve()/'visualization';out.mkdir(parents=True,exist_ok=True)
    a.model.init_lib(directory=str(out/'inputs'),threads=a.config.threads,output=False)
    try:
        if plots:sampled_plots(a,out)
        if vtk:return write_vtk(a,out)
    finally:a.model.finalize_lib()
