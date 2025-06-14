from collections import OrderedDict
from collections.abc import Callable, Iterable
from io import StringIO
from numbers import Integral, Real
from warnings import warn
import os
import tempfile

import h5py
import numpy as np

from openmc.mixin import EqualityMixin
from openmc.stats import Uniform
import openmc.checkvalue as cv
from . import HDF5_VERSION, HDF5_VERSION_MAJOR
from .ace import Table, get_metadata, get_table, Library
from .angle_distribution import AngleDistribution
from .angle_energy import AngleEnergy
from .correlated import CorrelatedAngleEnergy
from .data import ATOMIC_SYMBOL, EV_PER_MEV
from .endf import Evaluation, SUM_RULES, get_head_record, get_tab1_record,  get_cont_record
from .fission_energy import FissionEnergyRelease
from .function import Tabulated1D
from .njoy import make_ace_photonuclear
from .reaction import Reaction, REACTION_NAME, FISSION_MTS, _get_products, _get_fission_products_endf, _get_photon_products_endf, _get_activation_products
from .product import Product
from .energy_distribution import EnergyDistribution, LevelInelastic, \
    DiscretePhoton
from .function import Tabulated1D, Polynomial
from .kalbach_mann import KalbachMann
from .laboratory import LaboratoryAngleEnergy
from .nbody import NBodyPhaseSpace
from .product import Product
from .uncorrelated import UncorrelatedAngleEnergy

_REACTION_NAME = {key:value.replace("(n,","(gamma,") for key,value in REACTION_NAME.items()}

class PhotonuclearReaction(EqualityMixin):
    def __init__(self, mt):
        self._center_of_mass = True
        self._redundant = False
        self._q_value = 0.
        self._xs = None
        self._products = []
        self._derived_products = []
        self.mt = mt

    def __repr__(self):
        if self.mt in _REACTION_NAME:
            return f"<PhotonuclearReaction: MT={self.mt} {_REACTION_NAME[self.mt]}>"
        else:
            return f"<PhotonuclearReaction: MT={self.mt}>"
    
    @property
    def center_of_mass(self):
        return self._center_of_mass

    @center_of_mass.setter
    def center_of_mass(self, center_of_mass):
        cv.check_type('center of mass', center_of_mass, (bool, np.bool_))
        self._center_of_mass = center_of_mass

    @property
    def redundant(self):
        return self._redundant

    @redundant.setter
    def redundant(self, redundant):
        cv.check_type('redundant', redundant, (bool, np.bool_))
        self._redundant = redundant

    @property
    def q_value(self):
        return self._q_value

    @q_value.setter
    def q_value(self, q_value):
        cv.check_type('Q value', q_value, Real)
        self._q_value = q_value

    @property
    def products(self):
        return self._products

    @products.setter
    def products(self, products):
        cv.check_type('reaction products', products, Iterable, Product)
        self._products = products

    @property
    def derived_products(self):
        return self._derived_products

    @derived_products.setter
    def derived_products(self, derived_products):
        cv.check_type('reaction derived products', derived_products,
                      Iterable, Product)
        self._derived_products = derived_products

    @property
    def xs(self):
        return self._xs

    @xs.setter
    def xs(self, xs):
        cv.check_type("reaction cross section", xs, Callable)
        self._xs = xs
        
    @classmethod
    def from_endf(cls, ev, mt):
        """Generate a reaction from an ENDF evaluation

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        rx : openmc.data.PhotonuclearReaction
            Reaction data

        """
        rx = cls(mt)

        # Integrated cross section
        if (3, mt) in ev.section:
            file_obj = StringIO(ev.section[3, mt])
            get_head_record(file_obj)
            params, rx.xs = get_tab1_record(file_obj)
            rx.q_value = params[1]
            
        # Get fission product yields (nu) as well as delayed neutron energy
        # distributions
        if mt in FISSION_MTS:
            rx.products, rx.derived_products = _get_fission_products_endf(ev)        

        if (6, mt) in ev.section:
            # Product angle-energy distribution
            for product in _get_products(ev, mt):
                if mt in FISSION_MTS and product.particle == 'neutron':
                    rx.products[0].applicability = product.applicability
                    rx.products[0].distribution = product.distribution
                else:
                    rx.products.append(product)
        elif (4, mt) in ev.section or (5, mt) in ev.section:
            # Uncorrelated angle-energy distribution

            # Note that the energy distribution for MT=455 is read in
            # _get_fission_products_endf rather than here
            if (5, mt) in ev.section:
                product = Product('photon')
                file_obj = StringIO(ev.section[5, mt])
                items = get_head_record(file_obj)
                nk = items[4]
                for i in range(nk):
                    params, applicability = get_tab1_record(file_obj)
                    dist = UncorrelatedAngleEnergy()
                    dist.energy = EnergyDistribution.from_endf(file_obj, params)

                    product.applicability.append(applicability)
                    product.distribution.append(dist)
            elif mt == 2:
                # Elastic scattering -- no energy distribution is given since it
                # can be calulcated analytically
                product = Product('photon')
                dist = UncorrelatedAngleEnergy()
                product.distribution.append(dist)
            elif mt >= 50 and mt < 91:
                # Level inelastic scattering -- no energy distribution is given
                # since it can be calculated analytically. Here we determine the
                # necessary parameters to create a LevelInelastic object
                product = Product('neutron')
                dist = UncorrelatedAngleEnergy()

                A = ev.target['mass']
                threshold = abs(rx.q_value)
                mass_ratio = (A-1)/A
                dist.energy = LevelInelastic(threshold, mass_ratio)

                product.distribution.append(dist)

            if (4, mt) in ev.section:
                for dist in product.distribution:
                    dist.angle = AngleDistribution.from_endf(ev, mt)

            if mt in FISSION_MTS and (5, mt) in ev.section:
                # For fission reactions,
                rx.products[0].applicability = product.applicability
                rx.products[0].distribution = product.distribution
            else:
                rx.products.append(product)

        if (8, mt) in ev.section:
            rx.products += _get_activation_products(ev, rx)

        if (12, mt) in ev.section or (13, mt) in ev.section:
            rx.products += _get_photon_products_endf(ev, rx)
        return rx
    
    def to_hdf5(self, group):
        """Write reaction to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['mt'] = self.mt
        if self.mt in _REACTION_NAME:
            group.attrs['label'] = np.bytes_(_REACTION_NAME[self.mt])
        else:
            group.attrs['label'] = np.bytes_(self.mt)
        group.attrs['Q_value'] = self.q_value
        group.attrs['center_of_mass'] = 1 if self.center_of_mass else 0
        group.attrs['redundant'] = 1 if self.redundant else 0
        
        dset = group.create_dataset('xs', data=self.xs.y)
        threshold_idx = getattr(self.xs, '_threshold_idx', 0)
        dset.attrs['threshold_idx'] = threshold_idx
        
        for i, p in enumerate(self.products):
            pgroup = group.create_group(f'product_{i}')
            p.to_hdf5(pgroup)

    @classmethod
    def from_hdf5(cls, group, energy):
        """Generate reaction from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from
        energy : array
            array of energies at which cross sections are tabulated at.

        Returns
        -------
        openmc.data.Reaction
            Reaction data

        """

        mt = group.attrs['mt']
        rx = cls(mt)
        rx.q_value = group.attrs['Q_value']
        rx.center_of_mass = bool(group.attrs['center_of_mass'])
        rx.redundant = bool(group.attrs.get('redundant', False))
        xs = group['xs'][()]
        threshold_idx = group['xs'].attrs['threshold_idx']
        tabulated_xs = Tabulated1D(energy[threshold_idx:], xs)
        tabulated_xs._threshold_idx = threshold_idx
        rx.xs = tabulated_xs

        # Determine number of products
        n_product = 0
        for name in group:
            if name.startswith('product_'):
                n_product += 1

        # Read reaction products
        for i in range(n_product):
            pgroup = group[f'product_{i}']
            rx.products.append(Product.from_hdf5(pgroup))

        return rx

    @classmethod
    def from_ace(cls, ace, i_reaction):
        # Get nuclide energy grid
        n_grid = ace.nxs[3]
        grid = ace.xss[ace.jxs[1] : ace.jxs[1] + n_grid] * EV_PER_MEV

        if i_reaction > 0:
            mt = int(ace.xss[ace.jxs[6] + i_reaction - 1])
            rx = cls(mt)
            
            # Get Q-value of reaction
            rx.q_value = ace.xss[ace.jxs[7] + i_reaction - 1]*EV_PER_MEV

            # ==================================================================
            # CROSS SECTION

            # Get locator for cross-section data
            loc = int(ace.xss[ace.jxs[8] + i_reaction - 1])

            # Determine starting index on energy grid
            threshold_idx = int(ace.xss[ace.jxs[9] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energy = int(ace.xss[ace.jxs[9] + loc])
            energy = grid[threshold_idx : threshold_idx + n_energy]

            # Read reaction cross section
            xs = ace.xss[ace.jxs[9] + loc + 1 : ace.jxs[9] + loc + 1 + n_energy]

            # Fix negatives -- known issue for Y89 in JEFF 3.2
            if np.any(xs < 0.0):
                warn(
                    "Negative cross sections found for MT={} in {}. Setting "
                    "to zero.".format(rx.mt, ace.name)
                )
                xs[xs < 0.0] = 0.0

            tabulated_xs = Tabulated1D(energy, xs)
            tabulated_xs._threshold_idx = threshold_idx
            rx.xs = tabulated_xs
            
            # ==================================================================
            # YIELD AND ANGLE-ENERGY DISTRIBUTION
            
            for i_typ in range(ace.nxs[5]):
                loc = ace.jxs[10]+i_typ*ace.nxs[7]-1
                mts = ace.xss[int(ace.xss[loc+5]):int(ace.xss[loc+5])+int(ace.xss[loc+2])].astype(int)
                try:
                    i_mtr = mts.tolist().index(mt)
                except ValueError:
                    #skip products for other reactions
                    continue
                match int(ace.xss[loc+1]):
                    case 1:
                        particle=Product('neutron')
                    case 2:
                        particle=Product('photon')
                    case _:
                        #TODO: support more product particles
                        # for now skip particle if it is not a neutron or photon
                        warn(f"Unsupported secondary particle type in {ace.name} for MT={mt}. "
                            "This product will be skipped.")
                        continue
                        
                # Determine reference frame        
                ty = int(ace.xss[int(ace.xss[loc+6])+i_mtr])
                rx.center_of_mass = (ty < 0)
                
                # Determine multiplicity
                idx = int(ace.xss[loc+8])+int(ace.xss[int(ace.xss[loc+7])+i_mtr])-1
                match int(ace.xss[idx]):
                    case 6 | 16 | 12:
                        assert int(ace.xss[idx+1]) == mt
                        yield_ = Tabulated1D.from_ace(ace, idx+2)
                    case _:
                        raise NotImplementedError('partial yields not implemented yet')   
                particle.yield_ = yield_
                
                   
                # Determine locator for energy distribution  
                idx = int(ace.xss[int(ace.xss[loc+11])+i_mtr])   
                distribution = AngleEnergy.from_ace(ace, int(ace.xss[loc+12]), idx)
                
                # Determine locator for angular distribution
                idx = int(ace.xss[int(ace.xss[loc+9])+i_mtr])
                if idx==0:
                    # No angular distribution data are given for this reaction,
                    # isotropic scattering is assumed in LAB
                    energy = np.array([particle.yield_.x[0], particle.yield_.x[-1]])
                    mu_isotropic = Uniform(-1., 1.)
                    distribution.angle = AngleDistribution(
                    energy, [mu_isotropic, mu_isotropic])
                elif idx==-1:
                    pass
                else:
                    assert idx>=0
                    distribution.angle = AngleDistribution.from_ace(ace, int(ace.xss[loc+10]), idx)     
                    
                particle.distribution.append(distribution)
                rx.products.append(particle)
        else:
            # elastic cross section
            mt = 2
            rx = cls(mt)

            # Get elastic cross section values
            elastic_xs = ace.xss[ace.jxs[4] : ace.jxs[4] + ace.nxs[3]]

            # Fix negatives -- known issue for Ti46,49,50 in JEFF 3.2
            if np.any(elastic_xs < 0.0):
                warn(
                    "Negative elastic scattering cross section found for {}. "
                    "Setting to zero.".format(ace.name)
                )
                elastic_xs[elastic_xs < 0.0] = 0.0

            tabulated_xs = Tabulated1D(grid, elastic_xs)
            tabulated_xs._threshold_idx = 0
            rx.xs = tabulated_xs

        return rx

class IncidentPhotonuclear(EqualityMixin):
    """photo-nuclear interaction data.

    This class stores data derived from an ENDF-6 format photo-nuclear interaction
    sublibrary. Instances of this class are not normally instantiated by the
    user but rather created using the factory methods
    :meth:`Photonuclear.from_hdf5`, :meth:`Photonuclear.from_ace`, and
    :meth:`Photonuclear.from_endf`.

    Parameters
    ----------
    name : str
        Name of the nuclide using the GND naming convention
    atomic_number : int
        Number of photo-nuclears in the target nucleus
    atomic_number : int
        Number of photo-nuclears in the target nucleus
    mass_number : int
        Number of nucleons in the target nucleus
    metastable : int
        Metastable state of the target nucleus. A value of zero indicates ground
        state.
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide.

    Attributes
    ----------
    atomic_number : int
        Number of photo-nuclears in the target nucleus
    atomic_symbol : str
        Atomic symbol of the nuclide, e.g., 'Zr'
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide.
    fission_energy : None or openmc.data.FissionEnergyRelease
        The energy released by fission, tabulated by component (e.g. prompt
        neutrons or beta particles) and dependent on incident neutron energy        
    mass_number : int
        Number of nucleons in the target nucleus
    metastable : int
        Metastable state of the target nucleus. A value of zero indicates ground
        state.
    name : str
        Name of the nuclide using the GND naming convention
    reactions : collections.OrderedDict
        Contains the cross sections, secondary angle and energy distributions,
        and other associated data for each reaction. The keys are the MT values
        and the values are Reaction objects.

    """

    def __init__(
        self, name, atomic_number, mass_number, metastable, atomic_weight_ratio
    ):
        self.name = name
        self.atomic_number = atomic_number
        self.mass_number = mass_number
        self.reactions = OrderedDict()
        self.energy = []
        self._fission_energy = None
        self.metastable = metastable
        self.atomic_weight_ratio = atomic_weight_ratio

    def __contains__(self, mt):
        return mt in self.reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        else:
            raise KeyError("No reaction with MT={}.".format(mt))

    def __repr__(self):
        return "<IncidentPhotonuclear: {}>".format(self.name)

    def __iter__(self):
        return iter(self.reactions.values())

    @property
    def atomic_number(self):
        return self._atomic_number

    @property
    def name(self):
        return self._name

    @property
    def mass_number(self):
        return self._mass_number

    @property
    def metastable(self):
        return self._metastable

    @property
    def atomic_weight_ratio(self):
        return self._atomic_weight_ratio

    @property
    def atomic_symbol(self):
        return ATOMIC_SYMBOL[self.atomic_number]

    @name.setter
    def name(self, name):
        cv.check_type("name", name, str)
        self._name = name

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type("atomic number", atomic_number, Integral)
        cv.check_greater_than("atomic number", atomic_number, 0, True)
        self._atomic_number = atomic_number

    @mass_number.setter
    def mass_number(self, mass_number):
        cv.check_type("mass number", mass_number, Integral)
        cv.check_greater_than("mass number", mass_number, 0, True)
        self._mass_number = mass_number

    @metastable.setter
    def metastable(self, metastable):
        cv.check_type("metastable", metastable, Integral)
        cv.check_greater_than("metastable", metastable, 0, True)
        self._metastable = metastable

    @atomic_weight_ratio.setter
    def atomic_weight_ratio(self, atomic_weight_ratio):
        cv.check_type("atomic weight ratio", atomic_weight_ratio, Real)
        cv.check_greater_than("atomic weight ratio", atomic_weight_ratio, 0.0)
        self._atomic_weight_ratio = atomic_weight_ratio
        
    @property
    def fission_energy(self):
        return self._fission_energy

    @fission_energy.setter
    def fission_energy(self, fission_energy):
        cv.check_type('fission energy release', fission_energy,
                      FissionEnergyRelease)
        self._fission_energy = fission_energy        

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate incident photo-nuclear data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : openmc.data.endf.Evaluation or str
            ENDF evaluation to read from. If given as a string, it is assumed to
            be the filename for the ENDF file.


        Returns
        -------
        openmc.data.Photonuclear
            photo-nuclear interaction data

        """

        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        atomic_number = ev.target["atomic_number"]
        mass_number = ev.target["mass_number"]
        metastable = ev.target["isomeric_state"]
        atomic_weight_ratio = ev.target["mass"]

        element = ATOMIC_SYMBOL[atomic_number]
        name = "{}{}".format(element, mass_number)

        data = cls(name, atomic_number, mass_number, metastable, atomic_weight_ratio)

        # Read each reaction
        for mf, mt, nc, mod in ev.reaction_list:
            if mf == 3:
                data.reactions[mt] = PhotonuclearReaction.from_endf(ev, mt)
                
        # Read fission energy release (requires that we already know nu for
        # fission)
        if (1, 458) in ev.section:
            data.fission_energy = FissionEnergyRelease.from_endf(ev, data)

        data._evaluation = ev
        return data

    @classmethod
    def from_ace(cls, ace_or_filename, metastable_scheme="nndc"):
        """Generate incident photo-nuclear continuous-energy data from an ACE table

        Parameters
        ----------
        ace_or_filename : openmc.data.ace.Table or str
            ACE table to read from. If the value is a string, it is assumed to
            be the filename for the ACE file.
        metastable_scheme : {'nndc', 'mcnp'}
            Determine how ZAID identifiers are to be interpreted in the case of
            a metastable nuclide. Because the normal ZAID (=1000*Z + A) does not
            encode metastable information, different conventions are used among
            different libraries. In MCNP libraries, the convention is to add 400
            for a metastable nuclide except for Am242m, for which 95242 is
            metastable and 95642 (or 1095242 in newer libraries) is the ground
            state. For NNDC libraries, ZAID is given as 1000*Z + A + 100*m.

        Returns
        -------
        openmc.data.Photonuclear
            Incident photo-nuclear continuous-energy data

        """

        # First obtain the data for the first provided ACE table/file
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # If mass number hasn't been specified, make an educated guess
        zaid, xs = ace.name.split(".")
        if not xs.endswith("u"):
            raise TypeError(
                "{} is not a continuous-energy photo-nuclear ACE table.".format(ace)
            )
        name, element, Z, mass_number, metastable = get_metadata(
            int(zaid), metastable_scheme
        )

        data = cls(name, Z, mass_number, metastable, ace.atomic_weight_ratio)

        # Read energy grid
        n_energy = ace.nxs[3]
        i = ace.jxs[1]
        energy = ace.xss[i : i + n_energy] * EV_PER_MEV
        data.energy = energy

        total_xs = ace.xss[ace.jxs[2] : ace.jxs[2] + n_energy]
        nonelastic_xs = ace.xss[ace.jxs[3] : ace.jxs[3] + n_energy]
        elastic_xs = total_xs-nonelastic_xs 
        heating_number = ace.xss[ace.jxs[5] : ace.jxs[5] +  n_energy] * EV_PER_MEV
        
        # Create redundant reaction for total (MT=1)
        total = PhotonuclearReaction(1)
        total.xs = Tabulated1D(energy, total_xs)
        total.redundant = True
        data.reactions[1] = total
        
        # Create redundant reaction for nonelastic (MT=3)
        nonelastic = PhotonuclearReaction(3)
        nonelastic.xs = Tabulated1D(energy, nonelastic_xs)
        nonelastic.redundant = True
        data.reactions[3] = nonelastic

        # Create redundant reaction for heating (MT=301)
        heating = PhotonuclearReaction(301)
        heating.xs = Tabulated1D(energy, heating_number*total_xs)
        heating.redundant = True
        data.reactions[301] = heating
            
        # Read each reaction
        n_reaction = ace.nxs[4] + 1
        for i in range(n_reaction):
            if i==0:
                if ace.jxs[4]==0:
                    assert np.allclose(total_xs,nonelastic_xs)
                else:
                    raise NotImplementedError('photonuclear elastic scattering is not supported.')
            else:
                rx = PhotonuclearReaction.from_ace(ace, i)
                data.reactions[rx.mt] = rx

        return data

    def export_to_hdf5(self, path, mode="a", libver="earliest"):
        """Export incident photo-nuclear data to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        mode : {'r', 'r+', 'w', 'x', 'a'}
            Mode that is used to open the HDF5 file. This is the second argument
            to the :class:`h5py.File` constructor.
        libver : {'earliest', 'latest'}
            Compatibility mode for the HDF5 file. 'latest' will produce files
            that are less backwards compatible but have performance benefits.

        """

        # Open file and write version
        f = h5py.File(str(path), mode, libver=libver)
        f.attrs["filetype"] = np.bytes_("data_photonuclear")
        if "version" not in f.attrs:
            f.attrs["version"] = np.array(HDF5_VERSION)

        # If data come from ENDF, don't allow exporting to HDF5
        if hasattr(self, '_evaluation'):
            raise NotImplementedError('Cannot export incident photonuclear data that '
                                      'originated from an ENDF file.')
        # Write basic data
        g = f.create_group(self.name)
        g.attrs['Z'] = self.atomic_number
        g.attrs['A'] = self.mass_number
        g.attrs['metastable'] = self.metastable
        g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio

        # Determine union energy grid
        union_grid = np.array([])
        for rx in self:
            union_grid = np.union1d(union_grid, rx.xs.x)
            for product in rx.products:
                union_grid = np.union1d(union_grid, product.yield_.x)
        g.create_dataset("energy", data=union_grid)
        
        # Write reaction data
        rxs_group = g.create_group('reactions')
        for rx in self.reactions.values():
            # Skip writing redundant reaction if it doesn't have neutron
            # production or is a summed transmutation reaction.
            # Also write heating.
            if rx.redundant:
                neutron_rx = any(p.particle == 'neutron' for p in rx.products)
                keep_mts = (301,)
                if not (neutron_rx or rx.mt in keep_mts):
                    continue

            rx_group = rxs_group.create_group(f'reaction_{rx.mt:03}')
            rx.to_hdf5(rx_group)
        
        # Write fission energy release data
        if self.fission_energy is not None:
            fer_group = g.create_group('fission_energy_release')
            self.fission_energy.to_hdf5(fer_group)
        
        f.close()

    @classmethod
    def from_hdf5(cls, group_or_filename):
        """Generate photo-nuclear interaction data from HDF5 group

        Parameters
        ----------
        group_or_filename : h5py.Group or str
            HDF5 group containing interaction data. If given as a string, it is
            assumed to be the filename for the HDF5 file, and the first group is
            used to read from.

        Returns
        -------
        openmc.data.Photonuclear
            photo-nuclear interaction data

        """
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(str(group_or_filename), "r")

            # Make sure version matches
            if "version" in h5file.attrs:
                major, minor = h5file.attrs["version"]
                # For now all versions of HDF5 data can be read
            else:
                raise IOError(
                    "HDF5 data does not indicate a version. Your installation of "
                    "the OpenMC Python API expects version {}.x data.".format(
                        HDF5_VERSION_MAJOR
                    )
                )

            group = list(h5file.values())[0]
            h5file.close()

        name = group.name[1:]
        atomic_number = group.attrs["Z"]
        mass_number = group.attrs["A"]
        metastable = group.attrs["metastable"]
        atomic_weight_ratio = group.attrs["atomic_weight_ratio"]

        data = cls(name, atomic_number, mass_number, metastable, atomic_weight_ratio)

        # Read energy grid
        data.energy = group["energy"][()]

        # Read reaction data
        rxs_group = group["reactions"]
        for name, obj in sorted(rxs_group.items()):
            if name.startswith("reaction_"):
                rx = PhotonuclearReaction.from_hdf5(obj, data.energy)
                data.reactions[rx.mt] = rx

        # Read fission energy release data
        if 'fission_energy_release' in group:
            fer_group = group['fission_energy_release']
            data.fission_energy = FissionEnergyRelease.from_hdf5(fer_group)

        return data

    @classmethod
    def from_njoy(cls, filename, evaluation=None, **kwargs):
        """Generate incident photo-nuclear data by running NJOY.

        Parameters
        ----------
        filename : str
            Path to ENDF file
        evaluation : openmc.data.endf.Evaluation, optional
            If the ENDF file contains multiple material evaluations, this
            argument indicates which evaluation to use.
        **kwargs
            Keyword arguments passed to :func:`openmc.data.njoy.make_ace_photonuclear`

        Returns
        -------
        data : openmc.data.Photonuclear
            Incident photo-nuclear continuous-energy data

        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run NJOY to create an ACE library
            kwargs.setdefault("output_dir", tmpdir)
            for key in ("acer", "pendf"):
                kwargs.setdefault(key, os.path.join(kwargs["output_dir"], key))
            kwargs["evaluation"] = evaluation
            make_ace_photonuclear(filename, **kwargs)

            # Create instance from ACE tables within library
            lib = Library(kwargs["acer"])
            data = cls.from_ace(lib.tables[0])
            
            # Add fission energy release data
            ev = evaluation if evaluation is not None else Evaluation(filename)
            if (1, 458) in ev.section:
                data.fission_energy = FissionEnergyRelease.from_endf(ev, data)

        return data

    def get_reaction_components(self, mt):
        """Determine what reactions make up redundant reaction.

        Parameters
        ----------
        mt : int
            ENDF MT number of the reaction to find components of.

        Returns
        -------
        mts : list of int
            ENDF MT numbers of reactions that make up the redundant reaction and
            have cross sections provided.

        """
        mts = []
        if mt in SUM_RULES:
            for mt_i in SUM_RULES[mt]:
                mts += self.get_reaction_components(mt_i)
        if mts:
            return mts
        else:
            return [mt] if mt in self else []        

    def _get_redundant_reaction(self, mt, mts):
        """Create redundant reaction from its components

        Parameters
        ----------
        mt : int
            MT value of the desired reaction
        mts : iterable of int
            MT values of its components

        Returns
        -------
        openmc.Reaction
            Redundant reaction

        """

        rx = PhotonuclearReaction(mt)
        # Get energy grid
        energy = self.energy
        xss = [self.reactions[mt_i].xs for mt_i in mts]
        idx = min([xs._threshold_idx if hasattr(xs, '_threshold_idx')
                   else 0 for xs in xss])
        rx.xs = Tabulated1D(energy[idx:], Sum(xss)(energy[idx:]))
        rx.xs._threshold_idx = idx

        rx.redundant = True

        return rx
