from __future__ import division, unicode_literals
import sys
from collections import OrderedDict, Iterable, Mapping
from numbers import Integral, Real
from warnings import warn

import numpy as np
import h5py

from .data import ATOMIC_SYMBOL, SUM_RULES
from .ace import Table, get_table
from .fission_energy import FissionEnergyRelease
from .function import Tabulated1D, Sum
from .product import Product
from .reaction import Reaction, _get_photon_products
from .urr import ProbabilityTables
import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin

if sys.version_info[0] >= 3:
    basestring = str


class IncidentNeutron(EqualityMixin):
    """Continuous-energy neutron interaction data.

    Instances of this class are not normally instantiated by the user but rather
    created using the factory methods :meth:`IncidentNeutron.from_hdf5` and
    :meth:`IncidentNeutron.from_ace`.

    Parameters
    ----------
    name : str
        Name of the table
    atomic_number : int
        Number of protons in the nucleus
    mass_number : int
        Number of nucleons in the nucleus
    metastable : int
        Metastable state of the nucleus. A value of zero indicates ground state.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    temperature : float
        Temperature of the target nuclide in MeV.

    Attributes
    ----------
    atomic_number : int
        Number of protons in the nucleus
    atomic_symbol : str
        Atomic symbol of the nuclide, e.g., 'Zr'
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide.
    energy : numpy.ndarray
        The energy values (MeV) at which reaction cross-sections are tabulated.
    fission_energy : None or openmc.data.FissionEnergyRelease
        The energy released by fission, tabulated by component (e.g. prompt
        neutrons or beta particles) and dependent on incident neutron energy
    mass_number : int
        Number of nucleons in the nucleus
    metastable : int
        Metastable state of the nucleus. A value of zero indicates ground state.
    name : str
        ZAID identifier of the table, e.g. 92235.70c.
    reactions : collections.OrderedDict
        Contains the cross sections, secondary angle and energy distributions,
        and other associated data for each reaction. The keys are the MT values
        and the values are Reaction objects.
    summed_reactions : collections.OrderedDict
        Contains summed cross sections, e.g., the total cross section. The keys
        are the MT values and the values are Reaction objects.
    temperature : float
        Temperature of the target nuclide in MeV.
    urr : None or openmc.data.ProbabilityTables
        Unresolved resonance region probability tables

    """

    def __init__(self, name, atomic_number, mass_number, metastable,
                 atomic_weight_ratio, temperature):
        self.name = name
        self.atomic_number = atomic_number
        self.mass_number = mass_number
        self.metastable = metastable
        self.atomic_weight_ratio = atomic_weight_ratio
        self.temperature = temperature

        self._energy = None
        self._fission_energy = None
        self.reactions = OrderedDict()
        self.summed_reactions = OrderedDict()
        self.urr = None

    def __contains__(self, mt):
        return mt in self.reactions or mt in self.summed_reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        elif mt in self.summed_reactions:
            return self.summed_reactions[mt]
        else:
            raise KeyError('No reaction with MT={}.'.format(mt))

    def __repr__(self):
        return "<IncidentNeutron: {}>".format(self.name)

    def __iter__(self):
        return iter(self.reactions.values())

    @property
    def name(self):
        return self._name

    @property
    def atomic_number(self):
        return self._atomic_number

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
    def energy(self):
        return self._energy

    @property
    def fission_energy(self):
        return self._fission_energy

    @property
    def temperature(self):
        return self._temperature

    @property
    def reactions(self):
        return self._reactions

    @property
    def summed_reactions(self):
        return self._summed_reactions

    @property
    def urr(self):
        return self._urr

    @name.setter
    def name(self, name):
        cv.check_type('name', name, basestring)
        self._name = name

    @property
    def atomic_symbol(self):
        return atomic_symbol[self.atomic_number]

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type('atomic number', atomic_number, Integral)
        cv.check_greater_than('atomic number', atomic_number, 0)
        self._atomic_number = atomic_number

    @mass_number.setter
    def mass_number(self, mass_number):
        cv.check_type('mass number', mass_number, Integral)
        cv.check_greater_than('mass number', mass_number, 0, True)
        self._mass_number = mass_number

    @metastable.setter
    def metastable(self, metastable):
        cv.check_type('metastable', metastable, Integral)
        cv.check_greater_than('metastable', metastable, 0, True)
        self._metastable = metastable

    @atomic_weight_ratio.setter
    def atomic_weight_ratio(self, atomic_weight_ratio):
        cv.check_type('atomic weight ratio', atomic_weight_ratio, Real)
        cv.check_greater_than('atomic weight ratio', atomic_weight_ratio, 0.0)
        self._atomic_weight_ratio = atomic_weight_ratio

    @temperature.setter
    def temperature(self, temperature):
        cv.check_type('temperature', temperature, Real)
        cv.check_greater_than('temperature', temperature, 0.0, True)
        self._temperature = temperature

    @energy.setter
    def energy(self, energy):
        cv.check_type('energy grid', energy, Iterable, Real)
        self._energy = energy

    @fission_energy.setter
    def fission_energy(self, fission_energy):
        cv.check_type('fission energy release', fission_energy,
                      FissionEnergyRelease)
        self._fission_energy = fission_energy

    @reactions.setter
    def reactions(self, reactions):
        cv.check_type('reactions', reactions, Mapping)
        self._reactions = reactions

    @summed_reactions.setter
    def summed_reactions(self, summed_reactions):
        cv.check_type('summed reactions', summed_reactions, Mapping)
        self._summed_reactions = summed_reactions

    @urr.setter
    def urr(self, urr):
        cv.check_type('probability tables', urr,
                      (ProbabilityTables, type(None)))
        self._urr = urr

    def get_reaction_components(self, mt):
        """Determine what reactions make up summed reaction.

        Parameters
        ----------
        mt : int
            ENDF MT number of the reaction to find components of.

        Returns
        -------
        mts : list of int
            ENDF MT numbers of reactions that make up the summed reaction and
            have cross sections provided.

        """
        if mt in self.reactions:
            return [mt]
        elif mt in SUM_RULES:
            mts = SUM_RULES[mt]
        complete = False
        while not complete:
            new_mts = []
            complete = True
            for i, mt_i in enumerate(mts):
                if mt_i in self.reactions:
                    new_mts.append(mt_i)
                elif mt_i in SUM_RULES:
                    new_mts += SUM_RULES[mt_i]
                    complete = False
            mts = new_mts
        return mts

    def export_to_hdf5(self, path, mode='a'):
        """Export table to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        mode : {'r', r+', 'w', 'x', 'a'}
            Mode that is used to open the HDF5 file. This is the second argument
            to the :class:`h5py.File` constructor.

        """

        f = h5py.File(path, mode, libver='latest')

        # Write basic data
        g = f.create_group(self.name)
        g.attrs['Z'] = self.atomic_number
        g.attrs['A'] = self.mass_number
        g.attrs['metastable'] = self.metastable
        g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        g.attrs['temperature'] = self.temperature

        # Write energy grid
        g.create_dataset('energy', data=self.energy)

        # Write reaction data
        rxs_group = g.create_group('reactions')
        for rx in self.reactions.values():
            rx_group = rxs_group.create_group('reaction_{:03}'.format(rx.mt))
            rx.to_hdf5(rx_group)

            # Write total nu data if available
            if len(rx.derived_products) > 0 and 'total_nu' not in g:
                tgroup = g.create_group('total_nu')
                rx.derived_products[0].to_hdf5(tgroup)

        # Write unresolved resonance probability tables
        if self.urr is not None:
            urr_group = g.create_group('urr')
            self.urr.to_hdf5(urr_group)

        # Write fission energy release data
        if self.fission_energy is not None:
            fer_group = g.create_group('fission_energy_release')
            self.fission_energy.to_hdf5(fer_group)

        f.close()

    @classmethod
    def from_hdf5(cls, group_or_filename):
        """Generate continuous-energy neutron interaction data from HDF5 group

        Parameters
        ----------
        group_or_filename : h5py.Group or str
            HDF5 group containing interaction data. If given as a string, it is
            assumed to be the filename for the HDF5 file, and the first group is
            used to read from.

        Returns
        -------
        openmc.data.IncidentNeutron
            Continuous-energy neutron interaction data

        """
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(group_or_filename, 'r')
            group = list(h5file.values())[0]

        name = group.name[1:]
        atomic_number = group.attrs['Z']
        mass_number = group.attrs['A']
        metastable = group.attrs['metastable']
        atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        temperature = group.attrs['temperature']

        data = cls(name, atomic_number, mass_number, metastable,
                   atomic_weight_ratio, temperature)

        # Read energy grid
        data.energy = group['energy'].value

        # Read reaction data
        rxs_group = group['reactions']
        for name, obj in sorted(rxs_group.items()):
            if name.startswith('reaction_'):
                rx = Reaction.from_hdf5(obj, data.energy)
                data.reactions[rx.mt] = rx

                # Read total nu data if available
                if rx.mt in (18, 19, 20, 21, 38) and 'total_nu' in group:
                    tgroup = group['total_nu']
                    rx.derived_products.append(Product.from_hdf5(tgroup))

        # Build summed reactions.  Start from the highest MT number because high
        # MTs never depend on lower MTs.
        for mt_sum in sorted(SUM_RULES, reverse=True):
            if mt_sum not in data:
                xs_components = [data[mt].xs for mt in SUM_RULES[mt_sum]
                                 if mt in data]
                if len(xs_components) > 0:
                    rxn = Reaction(mt_sum)
                    rxn.xs = Sum(xs_components)
                    data.summed_reactions[mt_sum] = rxn

        # Read unresolved resonance probability tables
        if 'urr' in group:
            urr_group = group['urr']
            data.urr = ProbabilityTables.from_hdf5(urr_group)

        # Read fission energy release data
        if 'fission_energy_release' in group:
            fer_group = group['fission_energy_release']
            data.fission_energy = FissionEnergyRelease.from_hdf5(fer_group)

        return data

    @classmethod
    def from_ace(cls, ace_or_filename, metastable_scheme='nndc'):
        """Generate incident neutron continuous-energy data from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table or str
            ACE table to read from. If given as a string, it is assumed to be
            the filename for the ACE file.
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
        openmc.data.IncidentNeutron
            Incident neutron continuous-energy data

        """
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # If mass number hasn't been specified, make an educated guess
        zaid, xs = ace.name.split('.')
        zaid = int(zaid)
        Z = zaid // 1000
        mass_number = zaid % 1000

        if metastable_scheme == 'mcnp':
            if zaid > 1000000:
                # New SZA format
                Z = Z % 1000
                if zaid == 1095242:
                    metastable = 0
                else:
                    metastable = zaid // 1000000
            else:
                if zaid == 95242:
                    metastable = 1
                elif zaid == 95642:
                    metastable = 0
                else:
                    metastable = 1 if mass_number > 300 else 0
        elif metastable_scheme == 'nndc':
            metastable = 1 if mass_number > 300 else 0

        while mass_number > 3*Z:
            mass_number -= 100

        # Determine name for group
        element = ATOMIC_SYMBOL[Z]
        if metastable > 0:
            name = '{}{}_m{}.{}'.format(element, mass_number, metastable, xs)
        else:
            name = '{}{}.{}'.format(element, mass_number, xs)

        data = cls(name, Z, mass_number, metastable,
                   ace.atomic_weight_ratio, ace.temperature)

        # Read energy grid
        n_energy = ace.nxs[3]
        energy = ace.xss[ace.jxs[1]:ace.jxs[1] + n_energy]
        data.energy = energy
        total_xs = ace.xss[ace.jxs[1] + n_energy:ace.jxs[1] + 2*n_energy]
        absorption_xs = ace.xss[ace.jxs[1] + 2*n_energy:ace.jxs[1] + 3*n_energy]

        # Create summed reactions (total and absorption)
        total = Reaction(1)
        total.xs = Tabulated1D(energy, total_xs)
        data.summed_reactions[1] = total
        absorption = Reaction(27)
        absorption.xs = Tabulated1D(energy, absorption_xs)
        data.summed_reactions[27] = absorption

        # Read each reaction
        n_reaction = ace.nxs[4] + 1
        for i in range(n_reaction):
            rx = Reaction.from_ace(ace, i)
            data.reactions[rx.mt] = rx

        # Some photon production reactions may be assigned to MTs that don't
        # exist, usually MT=4. In this case, we create a new reaction and add
        # them
        n_photon_reactions = ace.nxs[6]
        photon_mts = ace.xss[ace.jxs[13]:ace.jxs[13] +
                             n_photon_reactions].astype(int)

        for mt in np.unique(photon_mts // 1000):
            if mt not in data:
                if mt not in SUM_RULES:
                    warn('Photon production is present for MT={} but no '
                         'cross section is given.'.format(mt))
                    continue

                # Create summed reaction with appropriate cross section
                rx = Reaction(mt)
                mts = data.get_reaction_components(mt)
                if len(mts) == 0:
                    warn('Photon production is present for MT={} but no '
                         'reaction components exist.'.format(mt))
                    continue
                rx.xs = Sum([data.reactions[mt_i].xs for mt_i in mts])

                # Determine summed cross section
                rx.products += _get_photon_products(ace, rx)
                data.summed_reactions[mt] = rx

        # Read unresolved resonance probability tables
        data.urr = ProbabilityTables.from_ace(ace)

        return data
