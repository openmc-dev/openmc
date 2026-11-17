"""Headless Matplotlib rendering of the exported VTK component surfaces.

The environment's OpenGL/EGL VTK RenderWindow crashes, so rendering uses Agg.
VTK still supplies the authoritative sampled surfaces and optional cutaway.
"""
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .export_vtk import COLORS


def render(vtm,output,cutaway=True):
    reader=vtk.vtkXMLMultiBlockDataReader();reader.SetFileName(str(vtm));reader.Update()
    blocks=reader.GetOutput();fig=plt.figure(figsize=(13,10));ax=fig.add_subplot(111,projection='3d')
    for j in range(blocks.GetNumberOfBlocks()):
        name=blocks.GetMetaData(j).Get(vtk.vtkCompositeDataSet.NAME())
        if name in ('external_vacuum','chamber_vacuum','blanket_gap'):continue
        data=blocks.GetBlock(j)
        if cutaway:
            plane=vtk.vtkPlane();plane.SetOrigin(0,0,0);plane.SetNormal(0,1,0)
            clip=vtk.vtkClipDataSet();clip.SetInputData(data);clip.SetClipFunction(plane);clip.Update();data=clip.GetOutput()
        surface=vtk.vtkDataSetSurfaceFilter();surface.SetInputData(data);surface.Update();poly=surface.GetOutput()
        if poly.GetNumberOfPoints()==0:continue
        points=vtk_to_numpy(poly.GetPoints().GetData())
        offsets=vtk_to_numpy(poly.GetPolys().GetOffsetsArray());indices=vtk_to_numpy(poly.GetPolys().GetConnectivityArray())
        faces=[points[indices[offsets[k]:offsets[k+1]]] for k in range(len(offsets)-1)]
        color=tuple(x/255 for x in COLORS.get(name,(110,110,110)))
        artist=Poly3DCollection(faces,facecolors=color,edgecolors='none',alpha=.2 if name=='plasma' else (.35 if name=='vessel' else 1.))
        ax.add_collection3d(artist)
    ax.set(xlim=(-1850,1850),ylim=(-1850,1850),zlim=(-1200,1200),xlabel='X [cm]',ylabel='Y [cm]',zlabel='Z [cm]')
    ax.set_box_aspect((1,1,.65));ax.view_init(elev=23,azim=-55)
    ax.set_title('Native OpenMC CSG: sampled '+('cutaway' if cutaway else 'full assembly')+'\nFull-reactor voxels underresolve thin first-wall layers')
    fig.tight_layout();fig.savefig(output,dpi=150);plt.close(fig)
