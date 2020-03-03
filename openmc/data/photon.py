from collections import OrderedDict
from collections.abc import Mapping, Callable
from copy import deepcopy
from io import StringIO
from numbers import Integral, Real
from math import pi, sqrt
import os

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.integrate import quad

from openmc.mixin import EqualityMixin
import openmc.checkvalue as cv
from . import HDF5_VERSION
from .ace import Table, get_metadata, get_table
from .data import ATOMIC_SYMBOL, EV_PER_MEV
from .endf import Evaluation, get_head_record, get_tab1_record, get_list_record
from .function import Tabulated1D


# Constants
MASS_ELECTRON_EV = 0.5109989461e6  # Electron mass energy
PLANCK_C = 1.2398419739062977e4  # Planck's constant times c in eV-Angstroms
FINE_STRUCTURE = 137.035999139  # Inverse fine structure constant
CM_PER_ANGSTROM = 1.0e-8
# classical electron radius in cm
R0 = CM_PER_ANGSTROM * PLANCK_C / (2.0 * pi * FINE_STRUCTURE * MASS_ELECTRON_EV)

# Electron subshell labels
_SUBSHELLS = [None, 'K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5',
              'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3',
              'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'P1', 'P2', 'P3', 'P4',
              'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'Q1', 'Q2', 'Q3']

_REACTION_NAME = {
    501: ('Total photon interaction', 'total'),
    502: ('Photon coherent scattering', 'coherent'),
    504: ('Photon incoherent scattering', 'incoherent'),
    515: ('Pair production, electron field', 'pair_production_electron'),
    516: ('Total pair production', 'pair_production_total'),
    517: ('Pair production, nuclear field', 'pair_production_nuclear'),
    522: ('Photoelectric absorption', 'photoelectric'),
    525: ('Heating', 'heating'),
    526: ('Electro-atomic scattering', 'electro_atomic_scat'),
    527: ('Electro-atomic bremsstrahlung', 'electro_atomic_brem'),
    528: ('Electro-atomic excitation', 'electro_atomic_excit'),
    534: ('K (1s1/2) subshell photoelectric', 'K'),
    535: ('L1 (2s1/2) subshell photoelectric', 'L1'),
    536: ('L2 (2p1/2) subshell photoelectric', 'L2'),
    537: ('L3 (2p3/2) subshell photoelectric', 'L3'),
    538: ('M1 (3s1/2) subshell photoelectric', 'M1'),
    539: ('M2 (3p1/2) subshell photoelectric', 'M2'),
    540: ('M3 (3p3/2) subshell photoelectric', 'M3'),
    541: ('M4 (3d3/2) subshell photoelectric', 'M4'),
    542: ('M5 (3d5/2) subshell photoelectric', 'M5'),
    543: ('N1 (4s1/2) subshell photoelectric', 'N1'),
    544: ('N2 (4p1/2) subshell photoelectric', 'N2'),
    545: ('N3 (4p3/2) subshell photoelectric', 'N3'),
    546: ('N4 (4d3/2) subshell photoelectric', 'N4'),
    547: ('N5 (4d5/2) subshell photoelectric', 'N5'),
    548: ('N6 (4f5/2) subshell photoelectric', 'N6'),
    549: ('N7 (4f7/2) subshell photoelectric', 'N7'),
    550: ('O1 (5s1/2) subshell photoelectric', 'O1'),
    551: ('O2 (5p1/2) subshell photoelectric', 'O2'),
    552: ('O3 (5p3/2) subshell photoelectric', 'O3'),
    553: ('O4 (5d3/2) subshell photoelectric', 'O4'),
    554: ('O5 (5d5/2) subshell photoelectric', 'O5'),
    555: ('O6 (5f5/2) subshell photoelectric', 'O6'),
    556: ('O7 (5f7/2) subshell photoelectric', 'O7'),
    557: ('O8 (5g7/2) subshell photoelectric', 'O8'),
    558: ('O9 (5g9/2) subshell photoelectric', 'O9'),
    559: ('P1 (6s1/2) subshell photoelectric', 'P1'),
    560: ('P2 (6p1/2) subshell photoelectric', 'P2'),
    561: ('P3 (6p3/2) subshell photoelectric', 'P3'),
    562: ('P4 (6d3/2) subshell photoelectric', 'P4'),
    563: ('P5 (6d5/2) subshell photoelectric', 'P5'),
    564: ('P6 (6f5/2) subshell photoelectric', 'P6'),
    565: ('P7 (6f7/2) subshell photoelectric', 'P7'),
    566: ('P8 (6g7/2) subshell photoelectric', 'P8'),
    567: ('P9 (6g9/2) subshell photoelectric', 'P9'),
    568: ('P10 (6h9/2) subshell photoelectric', 'P10'),
    569: ('P11 (6h11/2) subshell photoelectric', 'P11'),
    570: ('Q1 (7s1/2) subshell photoelectric', 'Q1'),
    571: ('Q2 (7p1/2) subshell photoelectric', 'Q2'),
    572: ('Q3 (7p3/2) subshell photoelectric', 'Q3')
}

# Compton profiles are read from a pre-generated HDF5 file when they are first
# needed. The dictionary stores an array of electron momentum values (at which
# the profiles are tabulated) with the key 'pz' and the profile for each element
# is a 2D array with shape (n_shells, n_momentum_values) stored on the key Z
_COMPTON_PROFILES = {}

# Scaled bremsstrahlung DCSs are read from a data file provided by Selzter and
# Berger when they are first needed. The dictionary stores an array of n
# incident electron kinetic energies with key 'electron_energies', an array of
# k reduced photon energies with key 'photon_energies', and the cross sections
# for each element are in a 2D array with shape (n, k) stored on the key 'Z'.
# It also stores data used for calculating the density effect correction and
# stopping power, namely, the mean excitation energy with the key 'I', number
# of electrons per subshell with the key 'num_electrons', and binding energies
# with the key 'ionization_energy'.
_BREMSSTRAHLUNG = {}


class AtomicRelaxation(EqualityMixin):
    """Atomic relaxation data.

    This class stores the binding energy, number of electrons, and electron
    transitions possible from ioniziation for each electron subshell of an
    atom. All of the data originates from an ENDF-6 atomic relaxation
    sub-library (NSUB=6). Instances of this class are not normally instantiated
    directly but rather created using the factory method
    :math:`AtomicRelaxation.from_endf`.

    Parameters
    ----------
    binding_energy : dict
        Dictionary indicating the binding energy in eV (values) for given
        subshells (keys). The subshells should be given as strings, e.g., 'K',
        'L1', 'L2', etc.
    num_electrons : dict
        Dictionary indicating the number of electrons in a subshell when neutral
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc.
    transitions : pandas.DataFrame
        Dictionary indicating allowed transitions and their probabilities
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc. The transitions are represented as
        a DataFrame with columns indicating the secondary and tertiary subshell,
        the energy of the transition in eV, and the fractional probability of
        the transition.

    Attributes
    ----------
    binding_energy : dict
        Dictionary indicating the binding energy in eV (values) for given
        subshells (keys). The subshells should be given as strings, e.g., 'K',
        'L1', 'L2', etc.
    num_electrons : dict
        Dictionary indicating the number of electrons in a subshell when neutral
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc.
    transitions : pandas.DataFrame
        Dictionary indicating allowed transitions and their probabilities
        (values) for given subshells (keys). The subshells should be given as
        strings, e.g., 'K', 'L1', 'L2', etc. The transitions are represented as
        a DataFrame with columns indicating the secondary and tertiary subshell,
        the energy of the transition in eV, and the fractional probability of
        the transition.

    See Also
    --------
    IncidentPhoton

    """
    def __init__(self, binding_energy, num_electrons, transitions):
        self.binding_energy = binding_energy
        self.num_electrons = num_electrons
        self.transitions = transitions
        self._e_fluorescence = {}

    @property
    def binding_energy(self):
        return self._binding_energy

    @property
    def num_electrons(self):
        return self._num_electrons

    @property
    def subshells(self):
        return list(sorted(self.binding_energy.keys()))

    @property
    def transitions(self):
        return self._transitions

    @binding_energy.setter
    def binding_energy(self, binding_energy):
        cv.check_type('binding energies', binding_energy, Mapping)
        for subshell, energy in binding_energy.items():
            cv.check_value('subshell', subshell, _SUBSHELLS)
            cv.check_type('binding energy', energy, Real)
            cv.check_greater_than('binding energy', energy, 0.0, True)
        self._binding_energy = binding_energy

    @num_electrons.setter
    def num_electrons(self, num_electrons):
        cv.check_type('number of electrons', num_electrons, Mapping)
        for subshell, num in num_electrons.items():
            cv.check_value('subshell', subshell, _SUBSHELLS)
            cv.check_type('number of electrons', num, Real)
            cv.check_greater_than('number of electrons', num, 0.0, True)
        self._num_electrons = num_electrons

    @transitions.setter
    def transitions(self, transitions):
        cv.check_type('transitions', transitions, Mapping)
        for subshell, df in transitions.items():
            cv.check_value('subshell', subshell, _SUBSHELLS)
            cv.check_type('transitions', df, pd.DataFrame)
        self._transitions = transitions

    @classmethod
    def from_ace(cls, ace):
        """Generate atomic relaxation data from an ACE file

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from

        Returns
        -------
        openmc.data.AtomicRelaxation
            Atomic relaxation data

        """
        # Create data dictionaries
        binding_energy = {}
        num_electrons = {}
        transitions = {}

        # Get shell designators
        n = ace.nxs[7]
        idx = ace.jxs[11]
        shells = [_SUBSHELLS[int(i)] for i in ace.xss[idx : idx+n]]

        # Get number of electrons for each shell
        idx = ace.jxs[12]
        for shell, num in zip(shells, ace.xss[idx : idx+n]):
            num_electrons[shell] = num

        # Get binding energy for each shell
        idx = ace.jxs[13]
        for shell, e in zip(shells, ace.xss[idx : idx+n]):
            binding_energy[shell] = e*EV_PER_MEV

        # Get transition table
        columns = ['secondary', 'tertiary', 'energy (eV)', 'probability']
        idx = ace.jxs[18]
        for i, subi in enumerate(shells):
            n_transitions = int(ace.xss[ace.jxs[15] + i])
            if n_transitions > 0:
                records = []
                for j in range(n_transitions):
                    subj = _SUBSHELLS[int(ace.xss[idx])]
                    subk = _SUBSHELLS[int(ace.xss[idx + 1])]
                    etr = ace.xss[idx + 2]*EV_PER_MEV
                    if j == 0:
                        ftr = ace.xss[idx + 3]
                    else:
                        ftr = ace.xss[idx + 3] - ace.xss[idx - 1]
                    records.append((subj, subk, etr, ftr))
                    idx += 4

                # Create dataframe for transitions
                transitions[subi] = pd.DataFrame.from_records(
                    records, columns=columns)

        return cls(binding_energy, num_electrons, transitions)

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate atomic relaxation data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : str or openmc.data.endf.Evaluation
            ENDF atomic relaxation evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.

        Returns
        -------
        openmc.data.AtomicRelaxation
            Atomic relaxation data

        """
        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        # Atomic relaxation data is always MF=28, MT=533
        if (28, 533) not in ev.section:
            raise IOError('{} does not appear to be an atomic relaxation '
                          'sublibrary.'.format(ev))

        # Determine number of subshells
        file_obj = StringIO(ev.section[28, 533])
        params = get_head_record(file_obj)
        n_subshells = params[4]

        # Create data dictionaries
        binding_energy = {}
        num_electrons = {}
        transitions = {}
        columns = ['secondary', 'tertiary', 'energy (eV)', 'probability']

        # Read data for each subshell
        for i in range(n_subshells):
            params, list_items = get_list_record(file_obj)
            subi = _SUBSHELLS[int(params[0])]
            n_transitions = int(params[5])
            binding_energy[subi] = list_items[0]
            num_electrons[subi] = list_items[1]

            if n_transitions > 0:
                # Read transition data
                records = []
                for j in range(n_transitions):
                    subj = _SUBSHELLS[int(list_items[6*(j+1)])]
                    subk = _SUBSHELLS[int(list_items[6*(j+1) + 1])]
                    etr = list_items[6*(j+1) + 2]
                    ftr = list_items[6*(j+1) + 3]
                    records.append((subj, subk, etr, ftr))

                # Create dataframe for transitions
                transitions[subi] = pd.DataFrame.from_records(
                    records, columns=columns)

        # Return instance of class
        return cls(binding_energy, num_electrons, transitions)

    @classmethod
    def from_hdf5(cls, group):
        """Generate atomic relaxation data from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.AtomicRelaxation
            Atomic relaxation data

        """
        # Create data dictionaries
        binding_energy = {}
        num_electrons = {}
        transitions = {}

        designators = [s.decode() for s in group.attrs['designators']]
        columns = ['secondary', 'tertiary', 'energy (eV)', 'probability']
        for shell in designators:
            # Shell group
            sub_group = group[shell]

            # Read subshell binding energy and number of electrons
            if 'binding_energy' in sub_group.attrs:
                binding_energy[shell] = sub_group.attrs['binding_energy']
            if 'num_electrons' in sub_group.attrs:
                num_electrons[shell] = sub_group.attrs['num_electrons']

            # Read transition data
            if 'transitions' in sub_group:
                df = pd.DataFrame(sub_group['transitions'][()],
                                  columns=columns)
                # Replace float indexes back to subshell strings
                df[columns[:2]] = df[columns[:2]].replace(
                              np.arange(float(len(_SUBSHELLS))), _SUBSHELLS)
                transitions[shell] = df

        return cls(binding_energy, num_electrons, transitions)

    def to_hdf5(self, group, shell):
        """Write atomic relaxation data to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        shell : str
            The subshell to write data for

        """

        # Write subshell binding energy and number of electrons
        group.attrs['binding_energy'] = self.binding_energy[shell]
        group.attrs['num_electrons'] = self.num_electrons[shell]

        # Write transition data with replacements
        if shell in self.transitions:
            df = self.transitions[shell].replace(
                 _SUBSHELLS, range(len(_SUBSHELLS)))
            group.create_dataset('transitions', data=df.values.astype(float))


class IncidentPhoton(EqualityMixin):
    r"""Photon interaction data.

    This class stores photo-atomic, photo-nuclear, atomic relaxation,
    Compton profile, stopping power, and bremsstrahlung data assembled from
    different sources. To create an instance, the factory method
    :meth:`IncidentPhoton.from_endf` can be used. To add atomic relaxation or
    Compton profile data, set the :attr:`IncidentPhoton.atomic_relaxation` and
    :attr:`IncidentPhoton.compton_profiles` attributes directly.

    Parameters
    ----------
    atomic_number : int
        Number of protons in the target nucleus

    Attributes
    ----------
    atomic_number : int
        Number of protons in the target nucleus
    atomic_relaxation : openmc.data.AtomicRelaxation or None
        Atomic relaxation data
    bremsstrahlung : dict
        Dictionary of bremsstrahlung data with keys 'I' (mean excitation energy
        in [eV]), 'num_electrons' (number of electrons in each subshell),
        'ionization_energy' (ionization potential of each subshell),
        'electron_energy' (incident electron kinetic energy values in [eV]),
        'photon_energy' (ratio of the energy of the emitted photon to the
        incident electron kinetic energy), and 'dcs' (cross section values in
        [b]). The cross sections are in scaled form: :math:`(\beta^2/Z^2) E_k
        (d\sigma/dE_k)`, where :math:`E_k` is the energy of the emitted photon.
        A negative number of electrons in a subshell indicates conduction
        electrons.
    compton_profiles : dict
        Dictionary of Compton profile data with keys 'num_electrons' (number of
        electrons in each subshell), 'binding_energy' (ionization potential of
        each subshell), and 'J' (Hartree-Fock Compton profile as a function of
        the projection of the electron momentum on the scattering vector,
        :math:`p_z` for each subshell). Note that subshell occupancies may not
        match the atomic relaxation data.
    reactions : collections.OrderedDict
        Contains the cross sections for each photon reaction. The keys are MT
        values and the values are instances of :class:`PhotonReaction`.

    """

    def __init__(self, atomic_number):
        self.atomic_number = atomic_number
        self._atomic_relaxation = None
        self.reactions = OrderedDict()
        self.compton_profiles = {}
        self.bremsstrahlung = {}

    def __contains__(self, mt):
        return mt in self.reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        else:
            raise KeyError('No reaction with MT={}.'.format(mt))

    def __repr__(self):
        return "<IncidentPhoton: {}>".format(self.name)

    def __iter__(self):
        return iter(self.reactions.values())

    @property
    def atomic_number(self):
        return self._atomic_number

    @property
    def atomic_relaxation(self):
        return self._atomic_relaxation

    @property
    def name(self):
        return ATOMIC_SYMBOL[self.atomic_number]

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type('atomic number', atomic_number, Integral)
        cv.check_greater_than('atomic number', atomic_number, 0, True)
        self._atomic_number = atomic_number

    @atomic_relaxation.setter
    def atomic_relaxation(self, atomic_relaxation):
        cv.check_type('atomic relaxation data', atomic_relaxation,
                      AtomicRelaxation)
        self._atomic_relaxation = atomic_relaxation

    @classmethod
    def from_ace(cls, ace_or_filename):
        """Generate incident photon data from an ACE table

        Parameters
        ----------
        ace_or_filename : str or openmc.data.ace.Table
            ACE table to read from. If given as a string, it is assumed to be
            the filename for the ACE file.

        Returns
        -------
        openmc.data.IncidentPhoton
            Photon interaction data

        """
        # First obtain the data for the first provided ACE table/file
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # Get atomic number based on name of ACE table
        zaid, xs = ace.name.split('.')
        if not xs.endswith('p'):
            raise TypeError("{} is not a photoatomic transport ACE table.".format(ace))
        Z = get_metadata(int(zaid))[2]

        # Read each reaction
        data = cls(Z)
        for mt in (502, 504, 515, 522, 525):
            data.reactions[mt] = PhotonReaction.from_ace(ace, mt)

        # Get heating cross sections [eV-barn] from factors [eV per collision]
        # by multiplying with total xs
        data.reactions[525].xs.y *= sum([data.reactions[mt].xs.y for mt in
                                         (502, 504, 515, 522)])

        # Compton profiles
        n_shell = ace.nxs[5]
        if n_shell != 0:
            # Get number of electrons in each shell
            idx = ace.jxs[6]
            data.compton_profiles['num_electrons'] = ace.xss[idx : idx+n_shell]

            # Get binding energy for each shell
            idx = ace.jxs[7]
            e = ace.xss[idx : idx+n_shell]*EV_PER_MEV
            data.compton_profiles['binding_energy'] = e

            # Create Compton profile for each electron shell
            profiles = []
            for k in range(n_shell):
                # Get number of momentum values and interpolation scheme
                loca = int(ace.xss[ace.jxs[9] + k])
                jj = int(ace.xss[ace.jxs[10] + loca - 1])
                m = int(ace.xss[ace.jxs[10] + loca])

                # Read momentum and PDF
                idx = ace.jxs[10] + loca + 1
                pz = ace.xss[idx : idx+m]
                pdf = ace.xss[idx+m : idx+2*m]

                # Create proflie function
                J_k = Tabulated1D(pz, pdf, [m], [jj])
                profiles.append(J_k)
            data.compton_profiles['J'] = profiles

        # Subshell photoelectric xs and atomic relaxation data
        if ace.nxs[7] > 0:
            data.atomic_relaxation = AtomicRelaxation.from_ace(ace)

            # Get subshell designators
            n_subshells = ace.nxs[7]
            idx = ace.jxs[11]
            designators = [int(i) for i in ace.xss[idx : idx+n_subshells]]

            # Get energy grid for subshell photoionization
            n_energy = ace.nxs[3]
            idx = ace.jxs[1]
            energy = np.exp(ace.xss[idx : idx+n_energy])*EV_PER_MEV

            # Get cross section for each subshell
            idx = ace.jxs[16]
            for d in designators:
                # Create photon reaction
                mt = 533 + d
                rx = PhotonReaction(mt)
                data.reactions[mt] = rx

                # Store cross section, determining threshold
                xs = ace.xss[idx : idx+n_energy].copy()
                nonzero = (xs != 0.0)
                xs[nonzero] = np.exp(xs[nonzero])
                threshold = np.where(xs > 0.0)[0][0]
                rx.xs = Tabulated1D(energy[threshold:], xs[threshold:],
                                    [n_energy - threshold], [5])
                idx += n_energy

                # Copy binding energy
                shell = _SUBSHELLS[d]
                e = data.atomic_relaxation.binding_energy[shell]
                rx.subshell_binding_energy = e
        else:
            raise ValueError("ACE table {} does not have subshell data. Only "
                             "newer ACE photoatomic libraries are supported "
                             "(e.g., eprdata14).".format(ace.name))

        # Add bremsstrahlung DCS data
        data._add_bremsstrahlung()

        return data

    @classmethod
    def from_endf(cls, photoatomic, relaxation=None):
        """Generate incident photon data from an ENDF evaluation

        Parameters
        ----------
        photoatomic : str or openmc.data.endf.Evaluation
            ENDF photoatomic data evaluation to read from. If given as a string,
            it is assumed to be the filename for the ENDF file.
        relaxation : str or openmc.data.endf.Evaluation, optional
            ENDF atomic relaxation data evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.

        Returns
        -------
        openmc.data.IncidentPhoton
            Photon interaction data

        """
        if isinstance(photoatomic, Evaluation):
            ev = photoatomic
        else:
            ev = Evaluation(photoatomic)

        Z = ev.target['atomic_number']
        data = cls(Z)

        # Read each reaction
        for mf, mt, nc, mod in ev.reaction_list:
            if mf == 23:
                data.reactions[mt] = PhotonReaction.from_endf(ev, mt)

        # Add atomic relaxation data if it hasn't been added already
        if relaxation is not None:
            data.atomic_relaxation = AtomicRelaxation.from_endf(relaxation)

        # If Compton profile data hasn't been loaded, do so
        if not _COMPTON_PROFILES:
            filename = os.path.join(os.path.dirname(__file__), 'compton_profiles.h5')
            with h5py.File(filename, 'r') as f:
                _COMPTON_PROFILES['pz'] = f['pz'][()]
                for i in range(1, 101):
                    group = f['{:03}'.format(i)]
                    num_electrons = group['num_electrons'][()]
                    binding_energy = group['binding_energy'][()]*EV_PER_MEV
                    J = group['J'][()]
                    _COMPTON_PROFILES[i] = {'num_electrons': num_electrons,
                                            'binding_energy': binding_energy,
                                            'J': J}

        # Add Compton profile data
        pz = _COMPTON_PROFILES['pz']
        profile = _COMPTON_PROFILES[Z]
        data.compton_profiles['num_electrons'] = profile['num_electrons']
        data.compton_profiles['binding_energy'] = profile['binding_energy']
        data.compton_profiles['J'] = [Tabulated1D(pz, J_k) for J_k in profile['J']]

        # Add bremsstrahlung DCS data
        data._add_bremsstrahlung()

        return data

    @classmethod
    def from_hdf5(cls, group_or_filename):
        """Generate photon reaction from an HDF5 group

        Parameters
        ----------
        group_or_filename : h5py.Group or str
            HDF5 group containing interaction data. If given as a string, it is
            assumed to be the filename for the HDF5 file, and the first group is
            used to read from.

        Returns
        -------
        openmc.data.IncidentPhoton
            Photon interaction data

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
                    'HDF5 data does not indicate a version. Your installation '
                    'of the OpenMC Python API expects version {}.x data.'
                    .format(HDF5_VERSION_MAJOR))

            group = list(h5file.values())[0]

        Z = group.attrs['Z']
        data = cls(Z)

        # Read energy grid
        energy = group['energy'][()]

        # Read cross section data
        for mt, (name, key) in _REACTION_NAME.items():
            if key in group:
                rgroup = group[key]
            elif key in group['subshells']:
                rgroup = group['subshells'][key]
            else:
                continue

            data.reactions[mt] = PhotonReaction.from_hdf5(rgroup, mt, energy)

        # Check for necessary reactions
        for mt in (502, 504, 522):
            assert mt in data, "Reaction {} not found".format(mt)

        # Read atomic relaxation
        data.atomic_relaxation = AtomicRelaxation.from_hdf5(group['subshells'])

        # Read Compton profiles
        if 'compton_profiles' in group:
            rgroup = group['compton_profiles']
            profile = data.compton_profiles
            profile['num_electrons'] = rgroup['num_electrons'][()]
            profile['binding_energy'] = rgroup['binding_energy'][()]

            # Get electron momentum values
            pz = rgroup['pz'][()]
            J = rgroup['J'][()]
            if pz.size != J.shape[1]:
                raise ValueError("'J' array shape is not consistent with the "
                                 "'pz' array shape")
            profile['J'] = [Tabulated1D(pz, Jk) for Jk in J]

        # Read bremsstrahlung
        if 'bremsstrahlung' in group:
            rgroup = group['bremsstrahlung']
            data.bremsstrahlung['I'] = rgroup.attrs['I']
            for key in ('dcs', 'electron_energy', 'ionization_energy',
                        'num_electrons', 'photon_energy'):
                data.bremsstrahlung[key] = rgroup[key][()]

        return data

    def export_to_hdf5(self, path, mode='a', libver='earliest'):
        """Export incident photon data to an HDF5 file.

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
        f.attrs['filetype'] = np.string_('data_photon')
        if 'version' not in f.attrs:
            f.attrs['version'] = np.array(HDF5_VERSION)

        group = f.create_group(self.name)
        group.attrs['Z'] = Z = self.atomic_number

        # Determine union energy grid
        union_grid = np.array([])
        for rx in self:
            union_grid = np.union1d(union_grid, rx.xs.x)
        group.create_dataset('energy', data=union_grid)

        # Write cross sections
        shell_group = group.create_group('subshells')
        designators = []
        for mt, rx in self.reactions.items():
            name, key = _REACTION_NAME[mt]
            if mt in (502, 504, 515, 517, 522, 525):
                sub_group = group.create_group(key)
            elif mt >= 534 and mt <= 572:
                # Subshell
                designators.append(key)
                sub_group = shell_group.create_group(key)

                # Write atomic relaxation
                if self.atomic_relaxation is not None:
                    if key in self.atomic_relaxation.subshells:
                        self.atomic_relaxation.to_hdf5(sub_group, key)
            else:
                continue

            rx.to_hdf5(sub_group, union_grid, Z)

        shell_group.attrs['designators'] = np.array(designators, dtype='S')

        # Write Compton profiles
        if self.compton_profiles:
            compton_group = group.create_group('compton_profiles')

            profile = self.compton_profiles
            compton_group.create_dataset('num_electrons',
                                         data=profile['num_electrons'])
            compton_group.create_dataset('binding_energy',
                                         data=profile['binding_energy'])

            # Get electron momentum values
            compton_group.create_dataset('pz', data=profile['J'][0].x)

            # Create/write 2D array of profiles
            J = np.array([Jk.y for Jk in profile['J']])
            compton_group.create_dataset('J', data=J)

        # Write bremsstrahlung
        if self.bremsstrahlung:
            brem_group = group.create_group('bremsstrahlung')
            for key, value in self.bremsstrahlung.items():
                if key == 'I':
                    brem_group.attrs[key] = value
                else:
                    brem_group.create_dataset(key, data=value)

    def _add_bremsstrahlung(self):
        """Add the data used in the thick-target bremsstrahlung approximation

        """
        # Load bremsstrahlung data if it has not yet been loaded
        if not _BREMSSTRAHLUNG:
            # Add data used for density effect correction
            filename = os.path.join(os.path.dirname(__file__), 'density_effect.h5')
            with h5py.File(filename, 'r') as f:
                for i in range(1, 101):
                    group = f['{:03}'.format(i)]
                    _BREMSSTRAHLUNG[i] = {
                        'I': group.attrs['I'],
                        'num_electrons': group['num_electrons'][()],
                        'ionization_energy': group['ionization_energy'][()]
                    }

            filename = os.path.join(os.path.dirname(__file__), 'BREMX.DAT')
            brem = open(filename, 'r').read().split()

            # Incident electron kinetic energy grid in eV
            _BREMSSTRAHLUNG['electron_energy'] = np.logspace(3, 9, 200)
            log_energy = np.log(_BREMSSTRAHLUNG['electron_energy'])

            # Get number of tabulated electron and photon energy values
            n = int(brem[37])
            k = int(brem[38])

            # Index in data
            p = 39

            # Get log of incident electron kinetic energy values, used for
            # cubic spline interpolation in log energy. Units are in MeV, so
            # convert to eV.
            logx = np.log(np.fromiter(brem[p:p+n], float, n)*EV_PER_MEV)
            p += n

            # Get reduced photon energy values
            _BREMSSTRAHLUNG['photon_energy'] = np.fromiter(brem[p:p+k], float, k)
            p += k

            for i in range(1, 101):
                dcs = np.empty([len(log_energy), k])

                # Get the scaled cross section values for each electron energy
                # and reduced photon energy for this Z. Units are in mb, so
                # convert to b.
                y = np.reshape(np.fromiter(brem[p:p+n*k], float, n*k), (n, k))*1.0e-3
                p += k*n

                for j in range(k):
                    # Cubic spline interpolation in log energy and linear DCS
                    cs = CubicSpline(logx, y[:, j])

                    # Get scaled DCS values (millibarns) on new energy grid
                    dcs[:, j] = cs(log_energy)

                _BREMSSTRAHLUNG[i]['dcs'] = dcs

        # Add bremsstrahlung DCS data
        self.bremsstrahlung['electron_energy'] = _BREMSSTRAHLUNG['electron_energy']
        self.bremsstrahlung['photon_energy'] = _BREMSSTRAHLUNG['photon_energy']
        self.bremsstrahlung.update(_BREMSSTRAHLUNG[self.atomic_number])


class PhotonReaction(EqualityMixin):
    """Photon-induced reaction

    Parameters
    ----------
    mt : int
        The ENDF MT number for this reaction.

    Attributes
    ----------
    anomalous_real : openmc.data.Tabulated1D
        Real part of the anomalous scattering factor
    anomlaous_imag : openmc.data.Tabulated1D
        Imaginary part of the anomalous scatttering factor
    mt : int
        The ENDF MT number for this reaction.
    scattering_factor : openmc.data.Tabulated1D
        Coherent or incoherent form factor.
    xs : Callable
        Cross section as a function of incident photon energy

    """

    def __init__(self, mt):
        self.mt = mt
        self._xs = None
        self._scattering_factor = None
        self._anomalous_real = None
        self._anomalous_imag = None

    def __repr__(self):
        if self.mt in _REACTION_NAME:
            return "<Photon Reaction: MT={} {}>".format(
                self.mt, _REACTION_NAME[self.mt][0])
        else:
            return "<Photon Reaction: MT={}>".format(self.mt)

    @property
    def anomalous_real(self):
        return self._anomalous_real

    @property
    def anomalous_imag(self):
        return self._anomalous_imag

    @property
    def scattering_factor(self):
        return self._scattering_factor

    @property
    def xs(self):
        return self._xs

    @anomalous_real.setter
    def anomalous_real(self, anomalous_real):
        cv.check_type('real part of anomalous scattering factor',
                      anomalous_real, Callable)
        self._anomalous_real = anomalous_real

    @anomalous_imag.setter
    def anomalous_imag(self, anomalous_imag):
        cv.check_type('imaginary part of anomalous scattering factor',
                      anomalous_imag, Callable)
        self._anomalous_imag = anomalous_imag

    @scattering_factor.setter
    def scattering_factor(self, scattering_factor):
        cv.check_type('scattering factor', scattering_factor, Callable)
        self._scattering_factor = scattering_factor

    @xs.setter
    def xs(self, xs):
        cv.check_type('reaction cross section', xs, Callable)
        self._xs = xs

    @classmethod
    def from_ace(cls, ace, mt):
        """Generate photon reaction from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        openmc.data.PhotonReaction
            Photon reaction data

        """
        # Create instance
        rx = cls(mt)

        # Get energy grid (stored as logarithms)
        n = ace.nxs[3]
        idx = ace.jxs[1]
        energy = np.exp(ace.xss[idx : idx+n])*EV_PER_MEV

        # Get index for appropriate reaction
        if mt == 502:
            # Coherent scattering
            idx = ace.jxs[1] + 2*n
        elif mt == 504:
            # Incoherent scattering
            idx = ace.jxs[1] + n
        elif mt == 515:
            # Pair production
            idx = ace.jxs[1] + 4*n
        elif mt == 522:
            # Photoelectric
            idx = ace.jxs[1] + 3*n
        elif mt == 525:
            # Heating
            idx = ace.jxs[5]
        else:
            raise ValueError('ACE photoatomic cross sections do not have '
                             'data for MT={}.'.format(mt))

        # Store cross section
        xs = ace.xss[idx : idx+n].copy()
        if mt == 525:
            # Get heating factors in [eV per collision]
            xs *= EV_PER_MEV
        else:
            nonzero = (xs != 0.0)
            xs[nonzero] = np.exp(xs[nonzero])
        rx.xs = Tabulated1D(energy, xs, [n], [5])

        # Get form factors for incoherent/coherent scattering
        new_format = (ace.nxs[6] > 0)
        if mt == 502:
            idx = ace.jxs[3]
            if new_format:
                n = (ace.jxs[4] - ace.jxs[3]) // 3
                x = ace.xss[idx : idx+n]
                idx += n
            else:
                x = np.array([
                    0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12,
                    0.15, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                    0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                    1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                    3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6,
                    5.8, 6.0])
                n = x.size
            ff = ace.xss[idx+n : idx+2*n]
            rx.scattering_factor = Tabulated1D(x, ff)

        elif mt == 504:
            idx = ace.jxs[2]
            if new_format:
                n = (ace.jxs[3] - ace.jxs[2]) // 2
                x = ace.xss[idx : idx+n]
                idx += n
            else:
                x = np.array([
                    0.0, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6,
                    0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 8.0
                ])
                n = x.size
            ff = ace.xss[idx : idx+n]
            rx.scattering_factor = Tabulated1D(x, ff)

        return rx

    @classmethod
    def from_endf(cls, ev, mt):
        """Generate photon reaction from an ENDF evaluation

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF photo-atomic interaction data evaluation
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        openmc.data.PhotonReaction
            Photon reaction data

        """
        rx = cls(mt)

        # Read photon cross section
        if (23, mt) in ev.section:
            file_obj = StringIO(ev.section[23, mt])
            get_head_record(file_obj)
            params, rx.xs = get_tab1_record(file_obj)

            # Set subshell binding energy and/or fluorescence yield
            if mt >= 534 and mt <= 599:
                rx.subshell_binding_energy = params[0]
            if mt >= 534 and mt <= 572:
                rx.fluorescence_yield = params[1]

        # Read form factors / scattering functions
        if (27, mt) in ev.section:
            file_obj = StringIO(ev.section[27, mt])
            get_head_record(file_obj)
            params, rx.scattering_factor = get_tab1_record(file_obj)

        # Check for anomalous scattering factor
        if mt == 502:
            if (27, 506) in ev.section:
                file_obj = StringIO(ev.section[27, 506])
                get_head_record(file_obj)
                params, rx.anomalous_real = get_tab1_record(file_obj)

            if (27, 505) in ev.section:
                file_obj = StringIO(ev.section[27, 505])
                get_head_record(file_obj)
                params, rx.anomalous_imag = get_tab1_record(file_obj)

        return rx

    @classmethod
    def from_hdf5(cls, group, mt, energy):
        """Generate photon reaction from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from
        mt : int
            The MT value of the reaction to get data for
        energy : Iterable of float
            arrays of energies at which cross sections are tabulated at

        Returns
        -------
        openmc.data.PhotonReaction
            Photon reaction data

        """
        # Create instance
        rx = cls(mt)

        # Cross sections
        xs = group['xs'][()]
        # Replace zero elements to small non-zero to enable log-log
        xs[xs == 0.0] = np.exp(-500.0)

        # Threshold
        threshold_idx = 0
        if 'threshold_idx' in group['xs'].attrs:
            threshold_idx = group['xs'].attrs['threshold_idx']

        # Store cross section
        rx.xs = Tabulated1D(energy[threshold_idx:], xs, [len(xs)], [5])

        # Check for anomalous scattering factor
        if 'anomalous_real' in group:
            rx.anomalous_real = Tabulated1D.from_hdf5(group['anomalous_real'])
        if 'anomalous_imag' in group:
            rx.anomalous_imag = Tabulated1D.from_hdf5(group['anomalous_imag'])

        # Check for factors / scattering functions
        if 'scattering_factor' in group:
            rx.scattering_factor = Tabulated1D.from_hdf5(group['scattering_factor'])

        return rx

    def to_hdf5(self, group, energy, Z):
        """Write photon reaction to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        energy : Iterable of float
            arrays of energies at which cross sections are tabulated at
        Z : int
            atomic number

        """

        # Write cross sections
        if self.mt >= 534 and self.mt <= 572:
            # Determine threshold
            threshold = self.xs.x[0]
            idx = np.searchsorted(energy, threshold, side='right') - 1

            # Interpolate cross section onto union grid and write
            photoionization = self.xs(energy[idx:])
            group.create_dataset('xs', data=photoionization)
            assert len(energy) == len(photoionization) + idx
            group['xs'].attrs['threshold_idx'] = idx
        else:
            group.create_dataset('xs', data=self.xs(energy))

        # Write scattering factor
        if self.scattering_factor is not None:
            if self.mt == 502:
                # Create integrated form factor
                ff = deepcopy(self.scattering_factor)
                ff.x *= ff.x
                ff.y *= ff.y/Z**2
                int_ff = Tabulated1D(ff.x, ff.integral())
                int_ff.to_hdf5(group, 'integrated_scattering_factor')
            self.scattering_factor.to_hdf5(group, 'scattering_factor')
        if self.anomalous_real is not None:
            self.anomalous_real.to_hdf5(group, 'anomalous_real')
        if self.anomalous_imag is not None:
            self.anomalous_imag.to_hdf5(group, 'anomalous_imag')
