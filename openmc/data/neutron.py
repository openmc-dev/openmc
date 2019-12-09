from collections import OrderedDict
from collections.abc import Mapping, MutableMapping
from io import StringIO
from math import log10
from numbers import Integral, Real
import os
import tempfile
from warnings import warn

import numpy as np
import h5py

from . import HDF5_VERSION, HDF5_VERSION_MAJOR
from .ace import Library, Table, get_table, get_metadata
from .data import ATOMIC_SYMBOL, K_BOLTZMANN, EV_PER_MEV
from .endf import (
    Evaluation, SUM_RULES, get_head_record, get_tab1_record, get_evaluations)
from .fission_energy import FissionEnergyRelease
from .function import Tabulated1D, Sum, ResonancesWithBackground
from .grid import linearize, thin
from .njoy import make_ace
from .product import Product
from .reaction import Reaction, _get_photon_products_ace
from . import resonance as res
from . import resonance_covariance as res_cov
from .urr import ProbabilityTables
import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin


# Fractions of resonance widths used for reconstructing resonances
_RESONANCE_ENERGY_GRID = np.logspace(-3, 3, 61)


class IncidentNeutron(EqualityMixin):
    """Continuous-energy neutron interaction data.

    This class stores data derived from an ENDF-6 format neutron interaction
    sublibrary. Instances of this class are not normally instantiated by the
    user but rather created using the factory methods
    :meth:`IncidentNeutron.from_hdf5`, :meth:`IncidentNeutron.from_ace`, and
    :meth:`IncidentNeutron.from_endf`.

    Parameters
    ----------
    name : str
        Name of the nuclide using the GND naming convention
    atomic_number : int
        Number of protons in the target nucleus
    mass_number : int
        Number of nucleons in the target nucleus
    metastable : int
        Metastable state of the target nucleus. A value of zero indicates ground
        state.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    kTs : Iterable of float
        List of temperatures of the target nuclide in the data set.
        The temperatures have units of eV.

    Attributes
    ----------
    atomic_number : int
        Number of protons in the target nucleus
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
    resonances : openmc.data.Resonances or None
        Resonance parameters
    resonance_covariance : openmc.data.ResonanceCovariance or None
        Covariance for resonance parameters
    temperatures : list of str
        List of string representations of the temperatures of the target nuclide
        in the data set.  The temperatures are strings of the temperature,
        rounded to the nearest integer; e.g., '294K'
    kTs : Iterable of float
        List of temperatures of the target nuclide in the data set.
        The temperatures have units of eV.
    urr : dict
        Dictionary whose keys are temperatures (e.g., '294K') and values are
        unresolved resonance region probability tables.

    """

    def __init__(self, name, atomic_number, mass_number, metastable,
                 atomic_weight_ratio, kTs):
        self.name = name
        self.atomic_number = atomic_number
        self.mass_number = mass_number
        self.metastable = metastable
        self.atomic_weight_ratio = atomic_weight_ratio
        self.kTs = kTs
        self.energy = {}
        self._fission_energy = None
        self.reactions = OrderedDict()
        self._urr = {}
        self._resonances = None

    def __contains__(self, mt):
        return mt in self.reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        else:
            # Try to create a redundant cross section
            mts = self.get_reaction_components(mt)
            if len(mts) > 0:
                return self._get_redundant_reaction(mt, mts)
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
    def fission_energy(self):
        return self._fission_energy

    @property
    def reactions(self):
        return self._reactions

    @property
    def resonances(self):
        return self._resonances

    @property
    def resonance_covariance(self):
        return self._resonance_covariance

    @property
    def urr(self):
        return self._urr

    @property
    def temperatures(self):
        return ["{}K".format(int(round(kT / K_BOLTZMANN))) for kT in self.kTs]

    @name.setter
    def name(self, name):
        cv.check_type('name', name, str)
        self._name = name

    @property
    def atomic_symbol(self):
        return ATOMIC_SYMBOL[self.atomic_number]

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type('atomic number', atomic_number, Integral)
        cv.check_greater_than('atomic number', atomic_number, 0, True)
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

    @fission_energy.setter
    def fission_energy(self, fission_energy):
        cv.check_type('fission energy release', fission_energy,
                      FissionEnergyRelease)
        self._fission_energy = fission_energy

    @reactions.setter
    def reactions(self, reactions):
        cv.check_type('reactions', reactions, Mapping)
        self._reactions = reactions

    @resonances.setter
    def resonances(self, resonances):
        cv.check_type('resonances', resonances, res.Resonances)
        self._resonances = resonances

    @resonance_covariance.setter
    def resonance_covariance(self, resonance_covariance):
        cv.check_type('resonance covariance', resonance_covariance,
                      res_cov.ResonanceCovariances)
        self._resonance_covariance = resonance_covariance

    @urr.setter
    def urr(self, urr):
        cv.check_type('probability table dictionary', urr, MutableMapping)
        for key, value in urr:
            cv.check_type('probability table temperature', key, str)
            cv.check_type('probability tables', value, ProbabilityTables)
        self._urr = urr

    def add_temperature_from_ace(self, ace_or_filename, metastable_scheme='nndc'):
        """Append data from an ACE file at a different temperature.

        Parameters
        ----------
        ace_or_filename : openmc.data.ace.Table or str
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

        """

        data = IncidentNeutron.from_ace(ace_or_filename, metastable_scheme)

        # Check if temprature already exists
        strT = data.temperatures[0]
        if strT in self.temperatures:
            warn('Cross sections at T={} already exist.'.format(strT))
            return

        # Check that name matches
        if data.name != self.name:
            raise ValueError('Data provided for an incorrect nuclide.')

        # Add temperature
        self.kTs += data.kTs

        # Add energy grid
        self.energy[strT] = data.energy[strT]

        # Add normal and redundant reactions
        for mt in data.reactions:
            if mt in self:
                self[mt].xs[strT] = data[mt].xs[strT]
            else:
                warn("Tried to add cross sections for MT={} at T={} but this "
                     "reaction doesn't exist.".format(mt, strT))

        # Add probability tables
        if strT in data.urr:
            self.urr[strT] = data.urr[strT]

    def add_elastic_0K_from_endf(self, filename, overwrite=False):
        """Append 0K elastic scattering cross section from an ENDF file.

        Parameters
        ----------
        filename : str
            Path to ENDF file
        overwrite : bool
            If existing 0 K data is present, this flag can be used to indicate
            that it should be overwritten. Otherwise, an exception will be
            thrown.

        Raises
        ------
        ValueError
            If 0 K data is already present and the `overwrite` parameter is
            False.

        """
        # Check for existing data
        if '0K' in self.energy and not overwrite:
            raise ValueError('0 K data already exists for this nuclide.')

        data = type(self).from_endf(filename)
        if data.resonances is not None:
            x = []
            y = []
            for rr in data.resonances:
                if isinstance(rr, res.RMatrixLimited):
                    raise TypeError('R-Matrix Limited not supported.')
                elif isinstance(rr, res.Unresolved):
                    continue

                # Get energies/widths for resonances
                e_peak = rr.parameters['energy'].values
                if isinstance(rr, res.MultiLevelBreitWigner):
                    gamma = rr.parameters['totalWidth'].values
                elif isinstance(rr, res.ReichMoore):
                    df = rr.parameters
                    gamma = (df['neutronWidth'] +
                             df['captureWidth'] +
                             abs(df['fissionWidthA']) +
                             abs(df['fissionWidthB'])).values

                # Determine peak energies and widths
                e_min, e_max = rr.energy_min, rr.energy_max
                in_range = (e_peak > e_min) & (e_peak < e_max)
                e_peak = e_peak[in_range]
                gamma = gamma[in_range]

                # Get midpoints between resonances (use min/max energy of
                # resolved region as absolute lower/upper bound)
                e_mid = np.concatenate(
                    ([e_min], (e_peak[1:] + e_peak[:-1])/2, [e_max]))

                # Add grid around each resonance that includes the peak +/- the
                # width times each value in _RESONANCE_ENERGY_GRID. Values are
                # constrained so that points around one resonance don't overlap
                # with points around another. This algorithm is from Fudge.
                energies = []
                for e, g, e_lower, e_upper in zip(e_peak, gamma, e_mid[:-1],
                                                  e_mid[1:]):
                    e_left = e - g*_RESONANCE_ENERGY_GRID
                    energies.append(e_left[e_left > e_lower][::-1])
                    e_right = e + g*_RESONANCE_ENERGY_GRID[1:]
                    energies.append(e_right[e_right < e_upper])

                # Concatenate all points
                energies = np.concatenate(energies)

                # Create 1000 equal log-spaced energies over RRR, combine with
                # resonance peaks and half-height energies
                e_log = np.logspace(log10(e_min), log10(e_max), 1000)
                energies = np.union1d(e_log, energies)

                # Linearize and thin cross section
                xi, yi = linearize(energies, data[2].xs['0K'])
                xi, yi = thin(xi, yi)

                # If there are multiple resolved resonance ranges (e.g. Pu239 in
                # ENDF/B-VII.1), combine them
                x = np.concatenate((x, xi))
                y = np.concatenate((y, yi))
        else:
            energies = data[2].xs['0K'].x
            x, y = linearize(energies, data[2].xs['0K'])
            x, y = thin(x, y)

        # Set 0K energy grid and elastic scattering cross section
        self.energy['0K'] = x
        self[2].xs['0K'] = Tabulated1D(x, y)

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

    def export_to_hdf5(self, path, mode='a', libver='earliest'):
        """Export incident neutron data to an HDF5 file.

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
        # If data come from ENDF, don't allow exporting to HDF5
        if hasattr(self, '_evaluation'):
            raise NotImplementedError('Cannot export incident neutron data that '
                                      'originated from an ENDF file.')

        # Open file and write version
        f = h5py.File(str(path), mode, libver=libver)
        f.attrs['filetype'] = np.string_('data_neutron')
        f.attrs['version'] = np.array(HDF5_VERSION)

        # Write basic data
        g = f.create_group(self.name)
        g.attrs['Z'] = self.atomic_number
        g.attrs['A'] = self.mass_number
        g.attrs['metastable'] = self.metastable
        g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        ktg = g.create_group('kTs')
        for i, temperature in enumerate(self.temperatures):
            ktg.create_dataset(temperature, data=self.kTs[i])

        # Write energy grid
        eg = g.create_group('energy')
        for temperature in self.temperatures:
            eg.create_dataset(temperature, data=self.energy[temperature])

        # Write 0K energy grid if needed
        if '0K' in self.energy and '0K' not in eg:
            eg.create_dataset('0K', data=self.energy['0K'])

        # Write reaction data
        rxs_group = g.create_group('reactions')
        for rx in self.reactions.values():
            # Skip writing redundant reaction if it doesn't have photon
            # production or is a summed transmutation reaction. MT=4 is also
            # sometimes needed for probability tables. Also write gas
            # production, heating, and damage energy production.
            if rx.redundant:
                photon_rx = any(p.particle == 'photon' for p in rx.products)
                keep_mts = (4, 16, 103, 104, 105, 106, 107,
                            203, 204, 205, 206, 207, 301, 444, 901)
                if not (photon_rx or rx.mt in keep_mts):
                    continue

            rx_group = rxs_group.create_group('reaction_{:03}'.format(rx.mt))
            rx.to_hdf5(rx_group)

            # Write total nu data if available
            if len(rx.derived_products) > 0 and 'total_nu' not in g:
                tgroup = g.create_group('total_nu')
                rx.derived_products[0].to_hdf5(tgroup)

        # Write unresolved resonance probability tables
        if self.urr:
            urr_group = g.create_group('urr')
            for temperature, urr in self.urr.items():
                tgroup = urr_group.create_group(temperature)
                urr.to_hdf5(tgroup)

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
            h5file = h5py.File(str(group_or_filename), 'r')

            # Make sure version matches
            if 'version' in h5file.attrs:
                major, minor = h5file.attrs['version']
                # For now all versions of HDF5 data can be read
            else:
                raise IOError(
                    'HDF5 data does not indicate a version. Your installation of '
                    'the OpenMC Python API expects version {}.x data.'
                    .format(HDF5_VERSION_MAJOR))

            group = list(h5file.values())[0]

        name = group.name[1:]
        atomic_number = group.attrs['Z']
        mass_number = group.attrs['A']
        metastable = group.attrs['metastable']
        atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        kTg = group['kTs']
        kTs = []
        for temp in kTg:
            kTs.append(kTg[temp][()])

        data = cls(name, atomic_number, mass_number, metastable,
                   atomic_weight_ratio, kTs)

        # Read energy grid
        e_group = group['energy']
        for temperature, dset in e_group.items():
            data.energy[temperature] = dset[()]

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

        # Read unresolved resonance probability tables
        if 'urr' in group:
            urr_group = group['urr']
            for temperature, tgroup in urr_group.items():
                data.urr[temperature] = ProbabilityTables.from_hdf5(tgroup)

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
        openmc.data.IncidentNeutron
            Incident neutron continuous-energy data

        """

        # First obtain the data for the first provided ACE table/file
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # If mass number hasn't been specified, make an educated guess
        zaid, xs = ace.name.split('.')
        if not xs.endswith('c'):
            raise TypeError(
                "{} is not a continuous-energy neutron ACE table.".format(ace))
        name, element, Z, mass_number, metastable = \
            get_metadata(int(zaid), metastable_scheme)

        # Assign temperature to the running list
        kTs = [ace.temperature*EV_PER_MEV]

        data = cls(name, Z, mass_number, metastable,
                   ace.atomic_weight_ratio, kTs)

        # Get string of temperature to use as a dictionary key
        strT = data.temperatures[0]

        # Read energy grid
        n_energy = ace.nxs[3]
        i = ace.jxs[1]
        energy = ace.xss[i : i + n_energy]*EV_PER_MEV
        data.energy[strT] = energy
        total_xs = ace.xss[i + n_energy : i + 2*n_energy]
        absorption_xs = ace.xss[i + 2*n_energy : i + 3*n_energy]
        heating_number = ace.xss[i + 4*n_energy : i + 5*n_energy]*EV_PER_MEV

        # Create redundant reactions (total, absorption, and heating)
        total = Reaction(1)
        total.xs[strT] = Tabulated1D(energy, total_xs)
        total.redundant = True
        data.reactions[1] = total

        if np.count_nonzero(absorption_xs) > 0:
            absorption = Reaction(101)
            absorption.xs[strT] = Tabulated1D(energy, absorption_xs)
            absorption.redundant = True
            data.reactions[101] = absorption

        heating = Reaction(301)
        heating.xs[strT] = Tabulated1D(energy, heating_number*total_xs)
        heating.redundant = True
        data.reactions[301] = heating

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

                # Create redundant reaction with appropriate cross section
                mts = data.get_reaction_components(mt)
                if len(mts) == 0:
                    warn('Photon production is present for MT={} but no '
                         'reaction components exist.'.format(mt))
                    continue

                # Determine redundant cross section
                rx = data._get_redundant_reaction(mt, mts)
                rx.products += _get_photon_products_ace(ace, rx)
                data.reactions[mt] = rx

        # For transmutation reactions, sometimes only individual levels are
        # present in an ACE file, e.g. MT=600-649 instead of the summation
        # MT=103. In this case, if a user wants to tally (n,p), OpenMC doesn't
        # know about the total cross section. Here, we explicitly create a
        # redundant reaction for this purpose.
        for mt in (16, 103, 104, 105, 106, 107):
            if mt not in data:
                # Determine if any individual levels are present
                mts = data.get_reaction_components(mt)
                if len(mts) == 0:
                    continue

                # Determine redundant cross section
                rx = data._get_redundant_reaction(mt, mts)
                data.reactions[mt] = rx

        # Make sure redundant cross sections that are present in an ACE file get
        # marked as such
        for rx in data:
            mts = data.get_reaction_components(rx.mt)
            if mts != [rx.mt]:
                rx.redundant = True
            if rx.mt in (203, 204, 205, 206, 207, 444):
                rx.redundant = True

        # Read unresolved resonance probability tables
        urr = ProbabilityTables.from_ace(ace)
        if urr is not None:
            data.urr[strT] = urr

        return data

    @classmethod
    def from_endf(cls, ev_or_filename, covariance=False):
        """Generate incident neutron continuous-energy data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : openmc.data.endf.Evaluation or str
            ENDF evaluation to read from. If given as a string, it is assumed to
            be the filename for the ENDF file.

        covariance : bool
            Flag to indicate whether or not covariance data from File 32 should be
            retrieved

        Returns
        -------
        openmc.data.IncidentNeutron
            Incident neutron continuous-energy data

        """
        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        atomic_number = ev.target['atomic_number']
        mass_number = ev.target['mass_number']
        metastable = ev.target['isomeric_state']
        atomic_weight_ratio = ev.target['mass']
        temperature = ev.target['temperature']

        # Determine name
        element = ATOMIC_SYMBOL[atomic_number]
        if metastable > 0:
            name = '{}{}_m{}'.format(element, mass_number, metastable)
        else:
            name = '{}{}'.format(element, mass_number)

        # Instantiate incident neutron data
        data = cls(name, atomic_number, mass_number, metastable,
                   atomic_weight_ratio, [temperature])

        if (2, 151) in ev.section:
            data.resonances = res.Resonances.from_endf(ev)

        if (32, 151) in ev.section and covariance:
            data.resonance_covariance = (
                res_cov.ResonanceCovariances.from_endf(ev, data.resonances)
            )

        # Read each reaction
        for mf, mt, nc, mod in ev.reaction_list:
            if mf == 3:
                data.reactions[mt] = Reaction.from_endf(ev, mt)

        # Replace cross sections for elastic, capture, fission
        try:
            if any(isinstance(r, res._RESOLVED) for r in data.resonances):
                for mt in (2, 102, 18):
                    if mt in data.reactions:
                        rx = data.reactions[mt]
                        rx.xs['0K'] = ResonancesWithBackground(
                            data.resonances, rx.xs['0K'], mt)
        except ValueError:
            # Thrown if multiple resolved ranges (e.g. Pu239 in ENDF/B-VII.1)
            pass

        # If first-chance, second-chance, etc. fission are present, check
        # whether energy distributions were specified in MF=5. If not, copy the
        # energy distribution from MT=18.
        for mt, rx in data.reactions.items():
            if mt in (19, 20, 21, 38):
                if (5, mt) not in ev.section:
                    if rx.products:
                        neutron = data.reactions[18].products[0]
                        rx.products[0].applicability = neutron.applicability
                        rx.products[0].distribution = neutron.distribution

        # Read fission energy release (requires that we already know nu for
        # fission)
        if (1, 458) in ev.section:
            data.fission_energy = FissionEnergyRelease.from_endf(ev, data)

        data._evaluation = ev
        return data

    @classmethod
    def from_njoy(cls, filename, temperatures=None, evaluation=None, **kwargs):
        """Generate incident neutron data by running NJOY.

        Parameters
        ----------
        filename : str
            Path to ENDF file
        temperatures : iterable of float
            Temperatures in Kelvin to produce data at. If omitted, data is
            produced at room temperature (293.6 K)
        evaluation : openmc.data.endf.Evaluation, optional
            If the ENDF file contains multiple material evaluations, this
            argument indicates which evaluation to use.
        **kwargs
            Keyword arguments passed to :func:`openmc.data.njoy.make_ace`

        Returns
        -------
        data : openmc.data.IncidentNeutron
            Incident neutron continuous-energy data

        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run NJOY to create an ACE library
            kwargs.setdefault("output_dir", tmpdir)
            for key in ("acer", "pendf", "heatr", "broadr", "gaspr", "purr"):
                kwargs.setdefault(key, os.path.join(kwargs["output_dir"], key))
            kwargs['evaluation'] = evaluation
            make_ace(filename, temperatures, **kwargs)

            # Create instance from ACE tables within library
            lib = Library(kwargs['acer'])
            data = cls.from_ace(lib.tables[0])
            for table in lib.tables[1:]:
                data.add_temperature_from_ace(table)

            # Add 0K elastic scattering cross section
            if '0K' not in data.energy:
                pendf = Evaluation(kwargs['pendf'])
                file_obj = StringIO(pendf.section[3, 2])
                get_head_record(file_obj)
                params, xs = get_tab1_record(file_obj)
                data.energy['0K'] = xs.x
                data[2].xs['0K'] = xs

            # Add fission energy release data
            ev = evaluation if evaluation is not None else Evaluation(filename)
            if (1, 458) in ev.section:
                data.fission_energy = f = FissionEnergyRelease.from_endf(ev, data)
            else:
                f = None

            # For energy deposition, we want to store two different KERMAs:
            # one calculated assuming outgoing photons deposit their energy
            # locally, and one calculated assuming they carry their energy
            # away. This requires two HEATR runs (which make_ace does by
            # default). Here, we just need to correct for the fact that NJOY
            # uses a fission heating number of h = EFR, whereas we want:
            #
            # 1) h = EFR + EGP + EGD + EB (for local case)
            # 2) h = EFR + EB (for non-local case)
            #
            # The best way to handle this is to subtract off the fission
            # KERMA that NJOY calculates and add back exactly what we want.

            # If NJOY is not run with HEATR at all, skip everything below
            if not kwargs["heatr"]:
                return data

            # Helper function to get a cross section from an ENDF file on a
            # given energy grid
            def get_file3_xs(ev, mt, E):
                file_obj = StringIO(ev.section[3, mt])
                get_head_record(file_obj)
                _, xs = get_tab1_record(file_obj)
                return xs(E)

            heating_local = Reaction(901)
            heating_local.redundant = True

            heatr_evals = get_evaluations(kwargs["heatr"])
            heatr_local_evals = get_evaluations(kwargs["heatr"] + "_local")
            for ev, ev_local in zip(heatr_evals, heatr_local_evals):
                temp = "{}K".format(round(ev.target["temperature"]))

                # Get total KERMA (originally from ACE file) and energy grid
                kerma = data.reactions[301].xs[temp]
                E = kerma.x

                if f is not None:
                    # Replace fission KERMA with (EFR + EB)*sigma_f
                    fission = data.reactions[18].xs[temp]
                    kerma_fission = get_file3_xs(ev, 318, E)
                    kerma.y = kerma.y - kerma_fission + (
                        f.fragments(E) + f.betas(E)) * fission(E)

                # For local KERMA, we first need to get the values from the
                # HEATR run with photon energy deposited locally and put
                # them on the same energy grid
                kerma_local = get_file3_xs(ev_local, 301, E)

                if f is not None:
                    # When photons deposit their energy locally, we replace the
                    # fission KERMA with (EFR + EGP + EGD + EB)*sigma_f
                    kerma_fission_local = get_file3_xs(ev_local, 318, E)
                    kerma_local = kerma_local - kerma_fission_local + (
                        f.fragments(E) + f.prompt_photons(E)
                        + f.delayed_photons(E) + f.betas(E))*fission(E)

                heating_local.xs[temp] = Tabulated1D(E, kerma_local)

            data.reactions[901] = heating_local

        return data

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

        rx = Reaction(mt)
        # Get energy grid
        for strT in self.temperatures:
            energy = self.energy[strT]
            xss = [self.reactions[mt_i].xs[strT] for mt_i in mts]
            idx = min([xs._threshold_idx if hasattr(xs, '_threshold_idx')
                       else 0 for xs in xss])
            rx.xs[strT] = Tabulated1D(energy[idx:], Sum(xss)(energy[idx:]))
            rx.xs[strT]._threshold_idx = idx

        rx.redundant = True

        return rx
