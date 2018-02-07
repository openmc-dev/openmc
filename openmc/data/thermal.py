from collections.abc import Iterable
from difflib import get_close_matches
from numbers import Real
import itertools
import os
import re
import shutil
import tempfile
from warnings import warn

import numpy as np
import h5py

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from . import HDF5_VERSION, HDF5_VERSION_MAJOR
from .data import K_BOLTZMANN, ATOMIC_SYMBOL, EV_PER_MEV, NATURAL_ABUNDANCE
from .ace import Table, get_table, Library
from .angle_energy import AngleEnergy
from .function import Tabulated1D
from .correlated import CorrelatedAngleEnergy
from .njoy import make_ace_thermal
from openmc.stats import Discrete, Tabular


_THERMAL_NAMES = {
    'c_Al27': ('al', 'al27'),
    'c_Be': ('be', 'be-metal'),
    'c_BeO': ('beo'),
    'c_Be_in_BeO': ('bebeo', 'be-o', 'be/o'),
    'c_C6H6': ('benz', 'c6h6'),
    'c_C_in_SiC': ('csic',),
    'c_Ca_in_CaH2': ('cah'),
    'c_D_in_D2O': ('dd2o', 'hwtr', 'hw'),
    'c_Fe56': ('fe', 'fe56'),
    'c_Graphite': ('graph', 'grph', 'gr'),
    'c_H_in_CaH2': ('hcah2'),
    'c_H_in_CH2': ('hch2', 'poly', 'pol'),
    'c_H_in_CH4_liquid': ('lch4', 'lmeth'),
    'c_H_in_CH4_solid': ('sch4', 'smeth'),
    'c_H_in_H2O': ('hh2o', 'lwtr', 'lw'),
    'c_H_in_H2O_solid': ('hice',),
    'c_H_in_C5O2H8': ('lucite', 'c5o2h8'),
    'c_H_in_YH2': ('hyh2'),
    'c_H_in_ZrH': ('hzrh', 'h-zr', 'h/zr', 'hzr'),
    'c_Mg24': ('mg', 'mg24'),
    'c_O_in_BeO': ('obeo', 'o-be', 'o/be'),
    'c_O_in_D2O': ('od2o'),
    'c_O_in_H2O_ice': ('oice'),
    'c_O_in_UO2': ('ouo2', 'o2-u', 'o2/u'),
    'c_ortho_D': ('orthod', 'dortho'),
    'c_ortho_H': ('orthoh', 'hortho'),
    'c_Si_in_SiC': ('sisic'),
    'c_SiO2_alpha': ('sio2', 'sio2a'),
    'c_SiO2_beta': ('sio2b'),
    'c_para_D': ('parad', 'dpara'),
    'c_para_H': ('parah', 'hpara'),
    'c_U_in_UO2': ('uuo2', 'u-o2', 'u/o2'),
    'c_Y_in_YH2': ('yyh2'),
    'c_Zr_in_ZrH': ('zrzrh', 'zr-h', 'zr/h')
}


def get_thermal_name(name):
    """Get proper S(a,b) table name, e.g. 'HH2O' -> 'c_H_in_H2O'

    Parameters
    ----------
    name : str
        Name of an ACE thermal scattering table

    Returns
    -------
    str
        GND-format thermal scattering name

    """
    if name in _THERMAL_NAMES:
        return name
    else:
        for proper_name, names in _THERMAL_NAMES.items():
            if name.lower() in names:
                return proper_name
        else:
            # Make an educated guess?? This actually works well for
            # JEFF-3.2 which stupidly uses names like lw00.32t,
            # lw01.32t, etc. for different temperatures

            # First, construct a list of all the values/keys in the names
            # dictionary
            all_names = itertools.chain(_THERMAL_NAMES.keys(),
                                        *_THERMAL_NAMES.values())

            matches = get_close_matches(name, all_names, cutoff=0.5)
            if len(matches) > 0:
                # Figure out the key for the corresponding match
                match = matches[0]
                if match not in _THERMAL_NAMES:
                    for key, value_list in _THERMAL_NAMES.items():
                        if match in value_list:
                            match = key
                            break

                warn('Thermal scattering material "{}" is not recognized. '
                     'Assigning a name of {}.'.format(name, match))
                return match
            else:
                # OK, we give up. Just use the ACE name.
                warn('Thermal scattering material "{0}" is not recognized. '
                     'Assigning a name of c_{0}.'.format(name))
                return 'c_' + name


class CoherentElastic(EqualityMixin):
    r"""Coherent elastic scattering data from a crystalline material

    Parameters
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in eV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} S_i`

    Attributes
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in eV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} S_i`

    """

    def __init__(self, bragg_edges, factors):
        self.bragg_edges = bragg_edges
        self.factors = factors

    def __call__(self, E):
        if isinstance(E, Iterable):
            E = np.asarray(E)
        idx = np.searchsorted(self.bragg_edges, E)
        return self.factors[idx] / E

    def __len__(self):
        return len(self.bragg_edges)

    @property
    def bragg_edges(self):
        return self._bragg_edges

    @property
    def factors(self):
        return self._factors

    @bragg_edges.setter
    def bragg_edges(self, bragg_edges):
        cv.check_type('Bragg edges', bragg_edges, Iterable, Real)
        self._bragg_edges = np.asarray(bragg_edges)

    @factors.setter
    def factors(self, factors):
        cv.check_type('structure factor cumulative sums', factors,
                      Iterable, Real)
        self._factors = np.asarray(factors)

    def to_hdf5(self, group, name):
        """Write coherent elastic scattering to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        name : str
            Name of the dataset to create

        """
        dataset = group.create_dataset(name, data=np.vstack(
            [self.bragg_edges, self.factors]))
        dataset.attrs['type'] = np.string_('bragg')

    @classmethod
    def from_hdf5(cls, dataset):
        """Read coherent elastic scattering from an HDF5 dataset

        Parameters
        ----------
        group : h5py.Dataset
            HDF5 group to write to

        Returns
        -------
        openmc.data.CoherentElastic
            Coherent elastic scattering cross section

        """
        bragg_edges = dataset.value[0, :]
        factors = dataset.value[1, :]
        return cls(bragg_edges, factors)


class ThermalScattering(EqualityMixin):
    """A ThermalScattering object contains thermal scattering data as represented by
    an S(alpha, beta) table.

    Parameters
    ----------
    name : str
        Name of the material using GND convention, e.g. c_H_in_H2O
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    kTs : Iterable of float
        List of temperatures of the target nuclide in the data set.
        The temperatures have units of eV.

    Attributes
    ----------
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    elastic_xs : openmc.data.Tabulated1D or openmc.data.CoherentElastic
        Elastic scattering cross section derived in the coherent or incoherent
        approximation
    inelastic_xs : openmc.data.Tabulated1D
        Inelastic scattering cross section derived in the incoherent
        approximation
    name : str
        Name of the material using GND convention, e.g. c_H_in_H2O
    temperatures : Iterable of str
        List of string representations the temperatures of the target nuclide
        in the data set.  The temperatures are strings of the temperature,
        rounded to the nearest integer; e.g., '294K'
    kTs : Iterable of float
        List of temperatures of the target nuclide in the data set.
        The temperatures have units of eV.
    nuclides : Iterable of str
        Nuclide names that the thermal scattering data applies to

    """

    def __init__(self, name, atomic_weight_ratio, kTs):
        self.name = name
        self.atomic_weight_ratio = atomic_weight_ratio
        self.kTs = kTs
        self.elastic_xs = {}
        self.elastic_mu_out = {}
        self.inelastic_xs = {}
        self.inelastic_e_out = {}
        self.inelastic_mu_out = {}
        self.inelastic_dist = {}
        self.secondary_mode = None
        self.nuclides = []

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<Thermal Scattering Data: {0}>".format(self.name)
        else:
            return "<Thermal Scattering Data>"

    @property
    def temperatures(self):
        return ["{}K".format(int(round(kT / K_BOLTZMANN))) for kT in self.kTs]

    def export_to_hdf5(self, path, mode='a', libver='earliest'):
        """Export table to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        mode : {'r', r+', 'w', 'x', 'a'}
            Mode that is used to open the HDF5 file. This is the second argument
            to the :class:`h5py.File` constructor.
        libver : {'earliest', 'latest'}
            Compatibility mode for the HDF5 file. 'latest' will produce files
            that are less backwards compatible but have performance benefits.

        """
        # Open file and write version
        f = h5py.File(path, mode, libver=libver)
        f.attrs['version'] = np.array(HDF5_VERSION)

        # Write basic data
        g = f.create_group(self.name)
        g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        g.attrs['nuclides'] = np.array(self.nuclides, dtype='S')
        g.attrs['secondary_mode'] = np.string_(self.secondary_mode)
        ktg = g.create_group('kTs')
        for i, temperature in enumerate(self.temperatures):
            ktg.create_dataset(temperature, data=self.kTs[i])

        for T in self.temperatures:
            Tg = g.create_group(T)
            # Write thermal elastic scattering
            if self.elastic_xs:
                elastic_group = Tg.create_group('elastic')

                self.elastic_xs[T].to_hdf5(elastic_group, 'xs')
                if self.elastic_mu_out:
                    elastic_group.create_dataset('mu_out',
                                                 data=self.elastic_mu_out[T])

            # Write thermal inelastic scattering
            if self.inelastic_xs:
                inelastic_group = Tg.create_group('inelastic')
                self.inelastic_xs[T].to_hdf5(inelastic_group, 'xs')
                if self.secondary_mode in ('equal', 'skewed'):
                    inelastic_group.create_dataset('energy_out',
                                                   data=self.inelastic_e_out[T])
                    inelastic_group.create_dataset('mu_out',
                                                   data=self.inelastic_mu_out[T])
                elif self.secondary_mode == 'continuous':
                    self.inelastic_dist[T].to_hdf5(inelastic_group)

        f.close()

    def add_temperature_from_ace(self, ace_or_filename, name=None):
        """Add data to the ThermalScattering object from an ACE file at a
        different temperature.

        Parameters
        ----------
        ace_or_filename : openmc.data.ace.Table or str
            ACE table to read from. If given as a string, it is assumed to be
            the filename for the ACE file.
        name : str
            GND-conforming name of the material, e.g. c_H_in_H2O. If none is
            passed, the appropriate name is guessed based on the name of the ACE
            table.

        Returns
        -------
        openmc.data.ThermalScattering
            Thermal scattering data

        """
        data = ThermalScattering.from_ace(ace_or_filename, name)

        # Check if temprature already exists
        strT = data.temperatures[0]
        if strT in self.temperatures:
            warn('S(a,b) data at T={} already exists.'.format(strT))
            return

        # Check that name matches
        if data.name != self.name:
            raise ValueError('Data provided for an incorrect material.')

        # Add temperature
        self.kTs += data.kTs

        # Add inelastic cross section and distributions
        if strT in data.inelastic_xs:
            self.inelastic_xs[strT] = data.inelastic_xs[strT]
        if strT in data.inelastic_e_out:
            self.inelastic_e_out[strT] = data.inelastic_e_out[strT]
        if strT in data.inelastic_mu_out:
            self.inelastic_mu_out[strT] = data.inelastic_mu_out[strT]
        if strT in data.inelastic_dist:
            self.inelastic_dist[strT] = data.inelastic_dist[strT]

        # Add elastic cross sectoin and angular distribution
        if strT in data.elastic_xs:
            self.elastic_xs[strT] = data.elastic_xs[strT]
        if strT in data.elastic_mu_out:
            self.elastic_mu_out[strT] = data.elastic_mu_out[strT]

    @classmethod
    def from_hdf5(cls, group_or_filename):
        """Generate thermal scattering data from HDF5 group

        Parameters
        ----------
        group_or_filename : h5py.Group or str
            HDF5 group containing interaction data. If given as a string, it is
            assumed to be the filename for the HDF5 file, and the first group
            is used to read from.

        Returns
        -------
        openmc.data.ThermalScattering
            Neutron thermal scattering data

        """
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(group_or_filename, 'r')

            # Make sure version matches
            if 'version' in h5file.attrs:
                major, minor = h5file.attrs['version']
                if major != HDF5_VERSION_MAJOR:
                    raise IOError(
                        'HDF5 data format uses version {}.{} whereas your '
                        'installation of the OpenMC Python API expects version '
                        '{}.x.'.format(major, minor, HDF5_VERSION_MAJOR))
            else:
                raise IOError(
                    'HDF5 data does not indicate a version. Your installation of '
                    'the OpenMC Python API expects version {}.x data.'
                    .format(HDF5_VERSION_MAJOR))

            group = list(h5file.values())[0]

        name = group.name[1:]
        atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        kTg = group['kTs']
        kTs = []
        for temp in kTg:
            kTs.append(kTg[temp].value)
        temperatures = [str(int(round(kT / K_BOLTZMANN))) + "K" for kT in kTs]

        table = cls(name, atomic_weight_ratio, kTs)
        table.nuclides = [nuc.decode() for nuc in group.attrs['nuclides']]
        table.secondary_mode = group.attrs['secondary_mode'].decode()

        # Read thermal elastic scattering
        for T in temperatures:
            Tgroup = group[T]
            if 'elastic' in Tgroup:
                elastic_group = Tgroup['elastic']

                # Cross section
                elastic_xs_type = elastic_group['xs'].attrs['type'].decode()
                if elastic_xs_type == 'Tabulated1D':
                    table.elastic_xs[T] = Tabulated1D.from_hdf5(
                        elastic_group['xs'])
                elif elastic_xs_type == 'bragg':
                    table.elastic_xs[T] = CoherentElastic.from_hdf5(
                        elastic_group['xs'])

                # Angular distribution
                if 'mu_out' in elastic_group:
                    table.elastic_mu_out[T] = elastic_group['mu_out'].value

            # Read thermal inelastic scattering
            if 'inelastic' in Tgroup:
                inelastic_group = Tgroup['inelastic']
                table.inelastic_xs[T] = Tabulated1D.from_hdf5(
                    inelastic_group['xs'])
                if table.secondary_mode in ('equal', 'skewed'):
                    table.inelastic_e_out[T] = inelastic_group['energy_out'].value
                    table.inelastic_mu_out[T] = inelastic_group['mu_out'].value
                elif table.secondary_mode == 'continuous':
                    table.inelastic_dist[T] = AngleEnergy.from_hdf5(
                        inelastic_group)

        return table

    @classmethod
    def from_ace(cls, ace_or_filename, name=None):
        """Generate thermal scattering data from an ACE table

        Parameters
        ----------
        ace_or_filename : openmc.data.ace.Table or str
            ACE table to read from. If given as a string, it is assumed to be
            the filename for the ACE file.
        name : str
            GND-conforming name of the material, e.g. c_H_in_H2O. If none is
            passed, the appropriate name is guessed based on the name of the ACE
            table.

        Returns
        -------
        openmc.data.ThermalScattering
            Thermal scattering data

        """
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # Get new name that is GND-consistent
        ace_name, xs = ace.name.split('.')
        name = get_thermal_name(ace_name)

        # Assign temperature to the running list
        kTs = [ace.temperature*EV_PER_MEV]
        temperatures = [str(int(round(ace.temperature*EV_PER_MEV
                                      / K_BOLTZMANN))) + "K"]

        table = cls(name, ace.atomic_weight_ratio, kTs)

        # Incoherent inelastic scattering cross section
        idx = ace.jxs[1]
        n_energy = int(ace.xss[idx])
        energy = ace.xss[idx+1 : idx+1+n_energy]*EV_PER_MEV
        xs = ace.xss[idx+1+n_energy : idx+1+2*n_energy]
        table.inelastic_xs[temperatures[0]] = Tabulated1D(energy, xs)

        if ace.nxs[7] == 0:
            table.secondary_mode = 'equal'
        elif ace.nxs[7] == 1:
            table.secondary_mode = 'skewed'
        elif ace.nxs[7] == 2:
            table.secondary_mode = 'continuous'

        n_energy_out = ace.nxs[4]
        if table.secondary_mode in ('equal', 'skewed'):
            n_mu = ace.nxs[3]
            idx = ace.jxs[3]
            table.inelastic_e_out[temperatures[0]] = \
                ace.xss[idx:idx + n_energy * n_energy_out * (n_mu + 2):
                        n_mu + 2]*EV_PER_MEV
            table.inelastic_e_out[temperatures[0]].shape = \
                (n_energy, n_energy_out)

            table.inelastic_mu_out[temperatures[0]] = \
                ace.xss[idx:idx + n_energy * n_energy_out * (n_mu + 2)]
            table.inelastic_mu_out[temperatures[0]].shape = \
                (n_energy, n_energy_out, n_mu+2)
            table.inelastic_mu_out[temperatures[0]] = \
                table.inelastic_mu_out[temperatures[0]][:, :, 1:]
        else:
            n_mu = ace.nxs[3] - 1
            idx = ace.jxs[3]
            locc = ace.xss[idx:idx + n_energy].astype(int)
            n_energy_out = \
                ace.xss[idx + n_energy:idx + 2 * n_energy].astype(int)
            energy_out = []
            mu_out = []
            for i in range(n_energy):
                idx = locc[i]

                # Outgoing energy distribution for incoming energy i
                e = ace.xss[idx + 1:idx + 1 + n_energy_out[i]*(n_mu + 3):
                            n_mu + 3]*EV_PER_MEV
                p = ace.xss[idx + 2:idx + 2 + n_energy_out[i]*(n_mu + 3):
                            n_mu + 3]/EV_PER_MEV
                c = ace.xss[idx + 3:idx + 3 + n_energy_out[i]*(n_mu + 3):
                            n_mu + 3]
                eout_i = Tabular(e, p, 'linear-linear', ignore_negative=True)
                eout_i.c = c

                # Outgoing angle distribution for each
                # (incoming, outgoing) energy pair
                mu_i = []
                for j in range(n_energy_out[i]):
                    mu = ace.xss[idx + 4:idx + 4 + n_mu]
                    p_mu = 1. / n_mu * np.ones(n_mu)
                    mu_ij = Discrete(mu, p_mu)
                    mu_ij.c = np.cumsum(p_mu)
                    mu_i.append(mu_ij)
                    idx += 3 + n_mu

                energy_out.append(eout_i)
                mu_out.append(mu_i)

            # Create correlated angle-energy distribution
            breakpoints = [n_energy]
            interpolation = [2]
            energy = table.inelastic_xs[temperatures[0]].x
            table.inelastic_dist[temperatures[0]] = CorrelatedAngleEnergy(
                breakpoints, interpolation, energy, energy_out, mu_out)

        # Incoherent/coherent elastic scattering cross section
        idx = ace.jxs[4]
        n_mu = ace.nxs[6] + 1
        if idx != 0:
            n_energy = int(ace.xss[idx])
            energy = ace.xss[idx + 1: idx + 1 + n_energy]*EV_PER_MEV
            P = ace.xss[idx + 1 + n_energy: idx + 1 + 2 * n_energy]

            if ace.nxs[5] == 4:
                # Coherent elastic
                table.elastic_xs[temperatures[0]] = CoherentElastic(
                    energy, P*EV_PER_MEV)

                # Coherent elastic shouldn't have angular distributions listed
                assert n_mu == 0
            else:
                # Incoherent elastic
                table.elastic_xs[temperatures[0]] = Tabulated1D(energy, P)

                # Angular distribution
                assert n_mu > 0
                idx = ace.jxs[6]
                table.elastic_mu_out[temperatures[0]] = \
                    ace.xss[idx:idx + n_energy * n_mu]
                table.elastic_mu_out[temperatures[0]].shape = \
                    (n_energy, n_mu)

        # Get relevant nuclides -- NJOY only allows one to specify three
        # nuclides that the S(a,b) table applies to. Thus, for all elements
        # other than H and Fe, we automatically add all the naturally-occurring
        # isotopes.
        for zaid, awr in ace.pairs:
            if zaid > 0:
                Z, A = divmod(zaid, 1000)
                element = ATOMIC_SYMBOL[Z]
                if element in ['H', 'Fe']:
                    table.nuclides.append(element + str(A))
                else:
                    if element + '0' not in table.nuclides:
                        table.nuclides.append(element + '0')
                    for isotope in sorted(NATURAL_ABUNDANCE):
                        if re.match(r'{}\d+'.format(element), isotope):
                            if isotope not in table.nuclides:
                                table.nuclides.append(isotope)

        return table

    @classmethod
    def from_njoy(cls, filename, filename_thermal, temperatures=None, **kwargs):
        """Generate incident neutron data by running NJOY.

        Parameters
        ----------
        filename : str
            Path to ENDF neutron sublibrary file
        filename_thermal : str
            Path to ENDF thermal scattering sublibrary file
        temperatures : iterable of float
            Temperatures in Kelvin to produce data at. If omitted, data is
            produced at all temperatures in the ENDF thermal scattering
            sublibrary.
        **kwargs
            Keyword arguments passed to :func:`openmc.data.njoy.make_ace_thermal`

        Returns
        -------
        data : openmc.data.ThermalScattering
            Thermal scattering data

        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run NJOY to create an ACE library
            ace_file = os.path.join(tmpdir, 'ace')
            xsdir_file = os.path.join(tmpdir, 'xsdir')
            make_ace_thermal(filename, filename_thermal, temperatures,
                             ace_file, xsdir_file, **kwargs)

            # Create instance from ACE tables within library
            lib = Library(ace_file)
            data = cls.from_ace(lib.tables[0])
            for table in lib.tables[1:]:
                data.add_temperature_from_ace(table)

        return data
