"""All lengths in centimetres; reference provenance is in references.md."""
from dataclasses import dataclass,asdict,field
from pathlib import Path
import json
import math


@dataclass(frozen=True)
class PFCoil:
    radius: float
    z: float
    width: float=100.
    height: float=65.


@dataclass(frozen=True)
class Config:
    # DEMO-2018-inspired, not a frozen engineering reference CAD configuration.
    major_radius: float=910.
    minor_radius: float=290.
    elongation: float=1.65
    triangularity: float=.33
    contour_points: int=32
    plasma_wall_gap: float=22.5
    tf_coils: int=16
    tf_minor_radius: float=590.
    tf_elongation: float=1.45
    tf_triangularity: float=.15
    tf_radial_build: float=75.
    tf_half_width: float=35.
    tf_case: float=6.
    tf_insulation: float=2.
    pf_coils: tuple[PFCoil,...]=field(default_factory=lambda:(
        PFCoil(430,1030),PFCoil(1050,1120),PFCoil(1750,450),
        PFCoil(1750,-450),PFCoil(1050,-1120),PFCoil(430,-1030)))
    cs_inner_radius: float=65.
    cs_outer_radius: float=165.
    cs_height: float=1080.
    cs_modules: int=6
    cs_gap: float=2.
    support_spine_radius: float=60.
    support_radii: tuple[float,float]=(172.,195.)
    support_height: float=1300.
    column_shield_radii: tuple[float,float]=(205.,225.)
    column_shield_height: float=1200.
    pf_case: float=6.
    # (R-R0, Z-Z_bottom) targets relative to the shaped plasma bottom.
    divertor_inner_target: tuple=((-180.,18.5),(-130.,-66.5))
    divertor_outer_target: tuple=((50.,28.5),(-55.,-66.5))
    divertor_armour: float=.8
    divertor_copper: float=1.6
    divertor_water: float=1.
    divertor_support: float=3.6
    cassette_radial_offsets: tuple[float,float]=(-170.,-10.)
    cassette_vertical_offsets: tuple[float,float]=(-101.5,-73.5)
    cassette_coolant_radial_offsets: tuple[float,float]=(-150.,-30.)
    cassette_coolant_vertical_offsets: tuple[float,float]=(-91.5,-83.5)
    cassette_liner: float=2.
    cs_case: float=5.
    cs_insulation: float=2.
    blanket_inboard: float=75.5
    blanket_outboard: float=98.2
    blanket_bottom: float=-400.
    toroidal_modules_per_sector: int=3
    poloidal_modules: int=16
    toroidal_gap: float=1.
    poloidal_gap_cm: float=1.
    module_sidewall: float=.5
    poloidal_sidewall_cm: float=.5
    armour: float=.2
    fw_front_steel: float=.8
    fw_cooling: float=1.
    fw_rear_steel: float=.2
    fw_helium_fraction: float=.6
    breeder_cycles: int=3
    breeder_fraction: float=.48
    multiplier_fraction: float=.45
    cooling_fraction: float=.07
    backplate: float=3.
    manifold: float=12.
    manifold_helium_fraction: float=.5
    bed_structure_fraction: float=.03
    cooling_helium_fraction: float=.6
    grouped_channels: int=6
    ribs_per_module: int=3
    enrichment: tuple[float,...]=(.6,.6,.6)
    li2o_solid_density: float=2.013
    breeder_packing: float=.64
    multiplier_packing: float=.64
    multiplier: str='Be'
    helium_pressure_mpa: float=8.
    purge_pressure_mpa: float=.2
    breeder_temperature: float=900.
    structure_temperature: float=600.
    coolant_temperature: float=673.15
    transport_magnet_temperature: float=293.6
    operating_magnet_temperature: float=4.5
    vessel_gap: float=10.
    vessel_inner_wall: float=3.6
    vessel_interwall: float=50.4
    vessel_outer_wall: float=6.
    vessel_steel_fraction: float=.6
    thermal_gap: float=10.
    thermal_shield: float=5.
    world_radius: float=2400.
    fidelity: int=2
    fusion_power_mw: float=2000.
    source_energy_ev: float=14.1e6
    source_energy_sigma_ev: float=0.
    radial_profile_power: float=0.
    poloidal_asymmetry: float=0.
    source_bank_size: int=100000
    particles: int=1000
    batches: int=10
    seed: int=20261106
    threads: int=2
    coupled_photons: bool=False
    mesh_dimension: tuple[int,int,int]=(32,32,24)
    vtk_dimension: tuple[int,int,int]=(96,96,80)
    vtk_origin: tuple[float,float,float]=(0.,0.,0.)
    vtk_extent: tuple[float,float,float]=(1900.,1900.,1250.)
    enable_mesh: bool=True
    enable_heating: bool=True
    enable_exposure: bool=True
    temperature_tolerance: float=1000.

    def validate(self):
        if not 0<self.minor_radius<self.major_radius or self.elongation<=0 or abs(self.triangularity)>=1:
            raise ValueError('Invalid plasma dimensions or shaping')
        if self.fidelity not in (1,2,3):raise ValueError('fidelity must be 1, 2 or 3')
        if self.contour_points<16 or self.contour_points%4:raise ValueError('contour_points must be a multiple of four >=16')
        if min(self.tf_coils,self.toroidal_modules_per_sector,self.poloidal_modules,self.breeder_cycles,self.grouped_channels,self.ribs_per_module,self.cs_modules)<1:
            raise ValueError('Component/module counts must be positive')
        if min(self.divertor_armour,self.divertor_copper,self.divertor_water,self.divertor_support,self.cassette_liner,self.pf_case,self.module_sidewall,self.poloidal_sidewall_cm,self.poloidal_gap_cm,self.toroidal_gap)<=0:
            raise ValueError('Structural thicknesses, coolant widths and gaps must be positive')
        if self.cs_height/self.cs_modules-self.cs_gap<=2*(self.cs_case+self.cs_insulation) or self.cs_inner_radius+2*(self.cs_case+self.cs_insulation)>=self.cs_outer_radius:
            raise ValueError('CS case/insulation consumes the winding pack')
        if self.cs_outer_radius>=self.support_radii[0]:raise ValueError('CS outer radius intersects the central support')
        if any(p.width<=2*self.pf_case or p.height<=2*self.pf_case or p.radius<=p.width/2 for p in self.pf_coils):raise ValueError('Invalid PF coil envelope/case thickness')
        if self.poloidal_modules<4 or self.poloidal_modules%2:raise ValueError('At least four poloidal modules required')
        for x in (self.breeder_packing,self.multiplier_packing,self.fw_helium_fraction,self.cooling_helium_fraction,self.bed_structure_fraction,self.manifold_helium_fraction,self.vessel_steel_fraction):
            if not 0<x<1:raise ValueError('Packing/mix fractions must be strictly between zero and one')
        if len(self.enrichment)!=self.breeder_cycles or any(not 0<=e<=1 for e in self.enrichment):raise ValueError('One valid Li6 enrichment per breeder cycle required')
        if not math.isclose(self.breeder_fraction+self.multiplier_fraction+self.cooling_fraction,1.):raise ValueError('Active layer fractions must sum to one')
        if self.multiplier not in ('Be','Be12Ti'):raise ValueError('Supported multiplier choices: Be, Be12Ti')
        if min(self.blanket_inboard,self.blanket_outboard)<=self.fw_total+self.backplate+self.manifold:raise ValueError('No remaining active inboard blanket depth')
        if min(self.armour,self.fw_front_steel,self.fw_cooling,self.fw_rear_steel,self.backplate,self.manifold,self.world_radius,self.fusion_power_mw,self.source_energy_ev)<=0:raise ValueError('Layer sizes, source power and energy must be positive')
        if self.source_energy_sigma_ev<0 or self.source_energy_sigma_ev*8>=self.source_energy_ev:raise ValueError('Source energy width must be nonnegative and small compared with mean energy')
        if abs(self.poloidal_asymmetry)>=1 or self.radial_profile_power<0:raise ValueError('Invalid source profile')
        if self.particles<1 or self.batches<2 or self.threads<1:raise ValueError('Positive histories/threads and >=2 batches required')
        if any(n<2 for n in self.mesh_dimension+self.vtk_dimension):raise ValueError('Mesh dimensions must be >=2')
        if self.toroidal_gap+2*self.module_sidewall>10:raise ValueError('Toroidal module gap/sidewalls are too large')
        return self

    @property
    def fw_total(self):return self.armour+self.fw_front_steel+self.fw_cooling+self.fw_rear_steel

    @property
    def neutron_rate(self):return self.fusion_power_mw*1e6/(17.6e6*1.602176634e-19)

    def save(self,path):Path(path).write_text(json.dumps(asdict(self),indent=2))

    @classmethod
    def load(cls,path):
        data=json.loads(Path(path).read_text())
        if 'pf_coils' in data:data['pf_coils']=tuple(PFCoil(**p) for p in data['pf_coils'])
        for k in ('enrichment','mesh_dimension','vtk_dimension','vtk_extent','vtk_origin','support_radii','column_shield_radii','divertor_inner_target','divertor_outer_target','cassette_radial_offsets','cassette_vertical_offsets','cassette_coolant_radial_offsets','cassette_coolant_vertical_offsets'):
            if k in data:data[k]=tuple(data[k])
        return cls(**data).validate()
