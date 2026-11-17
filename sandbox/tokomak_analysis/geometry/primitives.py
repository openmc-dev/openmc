"""Audited native Polygon RZ revolution and disjoint angular sectors."""
import numpy as np
import openmc


def revolved(points):
    points=np.asarray(points,float)
    if np.any(points[:,0]<0):raise ValueError('Revolved polygon radius cannot be negative')
    if len(points)==3:
        # Installed Polygon triangulation uses QJ, which needs four input
        # vertices; a collinear edge midpoint preserves the exact triangle.
        points=np.insert(points,1,(points[0]+points[1])/2,axis=0)
    return openmc.model.Polygon(points,basis='rz')


def sector(lo,hi,inset=0.):
    """Radians; a span below pi. Inset is a normal distance, not an angle."""
    if not 0<hi-lo<np.pi:raise ValueError('Sector must span between 0 and pi')
    p=openmc.Plane(a=-np.sin(lo),b=np.cos(lo),d=inset)
    q=openmc.Plane(a=-np.sin(hi),b=np.cos(hi),d=-inset)
    return +p & -q


def clip(points,normal,constant):
    """Sutherland-Hodgman clipping: normal . (R,Z) >= constant."""
    result=[];points=np.asarray(points)
    for a,b in zip(points,np.roll(points,-1,axis=0)):
        fa=np.dot(normal,a)-constant;fb=np.dot(normal,b)-constant
        if fa>=-1e-10:result.append(a)
        if (fa>0 and fb<0) or (fa<0 and fb>0):result.append(a+fa/(fa-fb)*(b-a))
    return np.asarray(result)


def rz_halfspace(normal,constant):
    """Exact signed RZ half-plane revolved around Z, for physical R>=0."""
    dx,dy=-np.asarray(normal,float);cc=float(constant)
    if abs(dx)<1e-12:
        surf=openmc.ZPlane(z0=-cc/dy)
        return +surf if dy<0 else -surf
    if abs(dy)<1e-12:
        radius=-cc/dx
        if radius<=0:
            # All physical R are on one side; a negative/zero cylinder is
            # neither needed nor legal. Use a sphere union as universal set.
            sph=openmc.Sphere(r=1.)
            whole=(-sph|+sph)
            return whole if dx<0 else ~whole
        surf=openmc.ZCylinder(r=radius)
        return +surf if dx<0 else -surf
    up=dy/dx<0
    surf=openmc.model.ZConeOneSided(z0=-cc/dy,r2=dy*dy/(dx*dx),up=up)
    return +surf if (up and dx-dy<0) or (not up and dx+dy<0) else -surf


def poloidal_sector(r0,lo,hi):
    normal=np.array([-np.sin(lo),np.cos(lo)])
    lower=rz_halfspace(normal,normal[0]*r0)
    normal=np.array([np.sin(hi),-np.cos(hi)])
    return lower&rz_halfspace(normal,normal[0]*r0)


def union(regions):return openmc.Union(list(regions))


def polygon_volume(points):
    """Exact full-revolution volume of a simple RZ polygon (Pappus)."""
    p=np.asarray(points);q=np.roll(p,-1,axis=0)
    cross=p[:,0]*q[:,1]-q[:,0]*p[:,1]
    return abs(2*np.pi*np.sum((p[:,0]+q[:,0])*cross)/6)


class Registry:
    def __init__(self):self.records={};self.cells={}
    def cell(self,name,region,material=None,component=None,**tags):
        cell=openmc.Cell(name=name,region=region,fill=material)
        self.cells[cell.id]=cell
        self.records[cell.id]=dict(cell_id=cell.id,name=name,component=component or name,
            material_id=material.id if isinstance(material,openmc.Material) else 0,**tags)
        return cell


def share_interfaces(geometry):
    """Deduplicate coincident surfaces within each universe coordinate frame.

    Global deduplication can conflate equal local planes in differently rotated
    universes. Shared local interfaces are essential for stable on-surface
    navigation through adjacent thin material layers.
    """
    for universe in geometry.get_all_universes().values():
        canonical={};replace={}
        for cell in universe.cells.values():
            for surface in cell.region.get_surfaces().values():
                key=(surface._type,surface.boundary_type,*[round(surface.coefficients[k],11) for k in surface._coeff_keys])
                if key in canonical:replace[surface.id]=canonical[key]
                else:canonical[key]=surface
        for cell in universe.cells.values():cell.region.remove_redundant_surfaces(replace)
