from collections import namedtuple
from collections.abc import Iterable
from io import StringIO
from math import log
from numbers import Real
import re
from warnings import warn

import numpy as np
from uncertainties import ufloat, unumpy, UFloat

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from .data import ATOMIC_SYMBOL, ATOMIC_NUMBER
from .endf import Evaluation, get_head_record, get_list_record, get_tab1_record


# Gives name and (change in A, change in Z) resulting from decay
_DECAY_MODES = {
    0: ('gamma', (0, 0)),
    1: ('beta-', (0, 1)),
    2: ('ec/beta+', (0, -1)),
    3: ('IT', (0, 0)),
    4: ('alpha', (-4, -2)),
    5: ('n', (-1, 0)),
    6: ('sf', None),
    7: ('p', (-1, -1)),
    8: ('e-', (0, 0)),
    9: ('xray', (0, 0)),
    10: ('unknown', None)
}

_RADIATION_TYPES = {
    0: 'gamma',
    1: 'beta-',
    2: 'ec/beta+',
    4: 'alpha',
    5: 'n',
    6: 'sf',
    7: 'p',
    8: 'e-',
    9: 'xray',
    10: 'anti-neutrino',
    11: 'neutrino'
}


def get_decay_modes(value):
    """Return sequence of decay modes given an ENDF RTYP value.

    Parameters
    ----------
    value : float
        ENDF definition of sequence of decay modes

    Returns
    -------
    list of str
        List of successive decays, e.g. ('beta-', 'neutron')

    """
    return [_DECAY_MODES[int(x)][0] for x in
            str(value).strip('0').replace('.', '')]


class FissionProductYields(EqualityMixin):
    """Independent and cumulative fission product yields.

    Parameters
    ----------
    ev_or_filename : str of openmc.data.endf.Evaluation
        ENDF fission product yield evaluation to read from. If given as a
        string, it is assumed to be the filename for the ENDF file.

    Attributes
    ----------
    cumulative : list of dict
        Cumulative yields for each tabulated energy. Each item in the list is a
        dictionary whose keys are nuclide names and values are cumulative
        yields. The i-th dictionary corresponds to the i-th incident neutron
        energy.
    energies : Iterable of float or None
        Energies at which fission product yields are tabulated.
    independent : list of dict
        Independent yields for each tabulated energy. Each item in the list is a
        dictionary whose keys are nuclide names and values are independent
        yields. The i-th dictionary corresponds to the i-th incident neutron
        energy.
    nuclide : dict
        Properties of the fissioning nuclide.

    Notes
    -----
    Neutron fission yields are typically not measured with a monoenergetic
    source of neutrons. As such, if the fission yields are given at, e.g.,
    0.0253 eV, one should interpret this as meaning that they are derived from a
    typical thermal reactor flux spectrum as opposed to a monoenergetic source
    at 0.0253 eV.

    """
    def __init__(self, ev_or_filename):
        # Define function that can be used to read both independent and
        # cumulative yields
        def get_yields(file_obj):
            # Determine number of energies
            n_energy = get_head_record(file_obj)[2]
            energies = np.zeros(n_energy)

            data = []
            for i in range(n_energy):
                # Determine i-th energy and number of products
                items, values = get_list_record(file_obj)
                energies[i] = items[0]
                n_products = items[5]

                # Get yields for i-th energy
                yields = {}
                for j in range(n_products):
                    Z, A = divmod(int(values[4*j]), 1000)
                    isomeric_state = int(values[4*j + 1])
                    name = ATOMIC_SYMBOL[Z] + str(A)
                    if isomeric_state > 0:
                        name += '_m{}'.format(isomeric_state)
                    yield_j = ufloat(values[4*j + 2], values[4*j + 3])
                    yields[name] = yield_j

                data.append(yields)

            return energies, data

        # Get evaluation if str is passed
        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        # Assign basic nuclide properties
        self.nuclide = {
            'name': ev.gnd_name,
            'atomic_number': ev.target['atomic_number'],
            'mass_number': ev.target['mass_number'],
            'isomeric_state': ev.target['isomeric_state']
        }

        # Read independent yields
        if (8, 454) in ev.section:
            file_obj = StringIO(ev.section[8, 454])
            self.energies, self.independent = get_yields(file_obj)

        # Read cumulative yields
        if (8, 459) in ev.section:
            file_obj = StringIO(ev.section[8, 459])
            energies, self.cumulative = get_yields(file_obj)
            assert np.all(energies == self.energies)

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate fission product yield data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : str or openmc.data.endf.Evaluation
            ENDF fission product yield evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.

        Returns
        -------
        openmc.data.FissionProductYields
            Fission product yield data

        """
        return cls(ev_or_filename)


class DecayMode(EqualityMixin):
    """Radioactive decay mode.

    Parameters
    ----------
    parent : str
        Parent decaying nuclide
    modes : list of str
        Successive decay modes
    daughter_state : int
        Metastable state of the daughter nuclide
    energy : uncertainties.UFloat
        Total decay energy in eV available in the decay process.
    branching_ratio : uncertainties.UFloat
        Fraction of the decay of the parent nuclide which proceeds by this mode.

    Attributes
    ----------
    branching_ratio : uncertainties.UFloat
        Fraction of the decay of the parent nuclide which proceeds by this mode.
    daughter : str
        Name of daughter nuclide produced from decay
    energy : uncertainties.UFloat
        Total decay energy in eV available in the decay process.
    modes : list of str
        Successive decay modes
    parent : str
        Parent decaying nuclide

    """

    def __init__(self, parent, modes, daughter_state, energy,
                 branching_ratio):
        self._daughter_state = daughter_state
        self.parent = parent
        self.modes = modes
        self.energy = energy
        self.branching_ratio = branching_ratio

    def __repr__(self):
        return ('<DecayMode: ({}), {} -> {}, {}>'.format(
            ','.join(self.modes), self.parent, self.daughter,
            self.branching_ratio))

    @property
    def branching_ratio(self):
        return self._branching_ratio

    @property
    def daughter(self):
        # Determine atomic number and mass number of parent
        symbol, A = re.match(r'([A-Zn][a-z]*)(\d+)', self.parent).groups()
        A = int(A)
        Z = ATOMIC_NUMBER[symbol]

        # Process changes
        for mode in self.modes:
            for name, changes in _DECAY_MODES.values():
                if name == mode:
                    if changes is not None:
                        delta_A, delta_Z = changes
                        A += delta_A
                        Z += delta_Z

        if self._daughter_state > 0:
            return '{}{}_m{}'.format(ATOMIC_SYMBOL[Z], A, self._daughter_state)
        else:
            return '{}{}'.format(ATOMIC_SYMBOL[Z], A)

    @property
    def energy(self):
        return self._energy

    @property
    def modes(self):
        return self._modes

    @property
    def parent(self):
        return self._parent

    @branching_ratio.setter
    def branching_ratio(self, branching_ratio):
        cv.check_type('branching ratio', branching_ratio, UFloat)
        cv.check_greater_than('branching ratio',
                              branching_ratio.nominal_value, 0.0, True)
        if branching_ratio.nominal_value == 0.0:
            warn('Decay mode {} of parent {} has a zero branching ratio.'
                 .format(self.modes, self.parent))
        cv.check_greater_than('branching ratio uncertainty',
                              branching_ratio.std_dev, 0.0, True)
        self._branching_ratio = branching_ratio

    @energy.setter
    def energy(self, energy):
        cv.check_type('decay energy', energy, UFloat)
        cv.check_greater_than('decay energy', energy.nominal_value, 0.0, True)
        cv.check_greater_than('decay energy uncertainty',
                              energy.std_dev, 0.0, True)
        self._energy = energy

    @modes.setter
    def modes(self, modes):
        cv.check_type('decay modes', modes, Iterable, str)
        self._modes = modes

    @parent.setter
    def parent(self, parent):
        cv.check_type('parent nuclide', parent, str)
        self._parent = parent


class Decay(EqualityMixin):
    """Radioactive decay data.

    Parameters
    ----------
    ev_or_filename : str of openmc.data.endf.Evaluation
        ENDF radioactive decay data evaluation to read from. If given as a
        string, it is assumed to be the filename for the ENDF file.

    Attributes
    ----------
    average_energies : dict
        Average decay energies in eV of each type of radiation for decay heat
        applications.
    decay_constant : uncertainties.UFloat
        Decay constant in inverse seconds.
    half_life : uncertainties.UFloat
        Half-life of the decay in seconds.
    modes : list
        Decay mode information for each mode of decay.
    nuclide : dict
        Dictionary describing decaying nuclide with keys 'name',
        'excited_state', 'mass', 'stable', 'spin', and 'parity'.
    spectra : dict
        Resulting radiation spectra for each radiation type.

    """
    def __init__(self, ev_or_filename):
        # Get evaluation if str is passed
        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        file_obj = StringIO(ev.section[8, 457])

        self.nuclide = {}
        self.modes = []
        self.spectra = {}
        self.average_energies = {}

        # Get head record
        items = get_head_record(file_obj)
        Z, A = divmod(items[0], 1000)
        metastable = items[3]
        self.nuclide['atomic_number'] = Z
        self.nuclide['mass_number'] = A
        self.nuclide['isomeric_state'] = metastable
        if metastable > 0:
            self.nuclide['name'] = '{}{}_m{}'.format(ATOMIC_SYMBOL[Z], A,
                                                     metastable)
        else:
            self.nuclide['name'] = '{}{}'.format(ATOMIC_SYMBOL[Z], A)
        self.nuclide['mass'] = items[1]  # AWR
        self.nuclide['excited_state'] = items[2]  # State of the original nuclide
        self.nuclide['stable'] = (items[4] == 1)  # Nucleus stability flag

        # Determine if radioactive/stable
        if not self.nuclide['stable']:
            NSP = items[5]  # Number of radiation types

            # Half-life and decay energies
            items, values = get_list_record(file_obj)
            self.half_life = ufloat(items[0], items[1])
            NC = items[4]//2
            pairs = [x for x in zip(values[::2], values[1::2])]
            ex = self.average_energies
            ex['light'] = ufloat(*pairs[0])
            ex['electromagnetic'] = ufloat(*pairs[1])
            ex['heavy'] = ufloat(*pairs[2])
            if NC == 17:
                ex['beta-'] = ufloat(*pairs[3])
                ex['beta+'] = ufloat(*pairs[4])
                ex['auger'] = ufloat(*pairs[5])
                ex['conversion'] = ufloat(*pairs[6])
                ex['gamma'] = ufloat(*pairs[7])
                ex['xray'] = ufloat(*pairs[8])
                ex['Bremsstrahlung'] = ufloat(*pairs[9])
                ex['annihilation'] = ufloat(*pairs[10])
                ex['alpha'] = ufloat(*pairs[11])
                ex['recoil'] = ufloat(*pairs[12])
                ex['SF'] = ufloat(*pairs[13])
                ex['neutron'] = ufloat(*pairs[14])
                ex['proton'] = ufloat(*pairs[15])
                ex['neutrino'] = ufloat(*pairs[16])

            items, values = get_list_record(file_obj)
            spin = items[0]
            if spin == -77.777:
                self.nuclide['spin'] = None
            else:
                self.nuclide['spin'] = spin
            self.nuclide['parity'] = items[1]  # Parity of the nuclide

            # Decay mode information
            n_modes = items[5]  # Number of decay modes
            for i in range(n_modes):
                decay_type = get_decay_modes(values[6*i])
                isomeric_state = int(values[6*i + 1])
                energy = ufloat(*values[6*i + 2:6*i + 4])
                branching_ratio = ufloat(*values[6*i + 4:6*(i + 1)])

                mode = DecayMode(self.nuclide['name'], decay_type, isomeric_state,
                                 energy, branching_ratio)
                self.modes.append(mode)

            discrete_type = {0.0: None, 1.0: 'allowed', 2.0: 'first-forbidden',
                             3.0: 'second-forbidden', 4.0: 'third-forbidden',
                             5.0: 'fourth-forbidden', 6.0: 'fifth-forbidden'}

            # Read spectra
            for i in range(NSP):
                spectrum = {}

                items, values = get_list_record(file_obj)
                # Decay radiation type
                spectrum['type'] = _RADIATION_TYPES[items[1]]
                # Continuous spectrum flag
                spectrum['continuous_flag'] = {0: 'discrete', 1: 'continuous',
                                               2: 'both'}[items[2]]
                spectrum['discrete_normalization'] = ufloat(*values[0:2])
                spectrum['energy_average'] = ufloat(*values[2:4])
                spectrum['continuous_normalization'] = ufloat(*values[4:6])

                NER = items[5]  # Number of tabulated discrete energies

                if not spectrum['continuous_flag'] == 'continuous':
                    # Information about discrete spectrum
                    spectrum['discrete'] = []
                    for j in range(NER):
                        items, values = get_list_record(file_obj)
                        di = {}
                        di['energy'] = ufloat(*items[0:2])
                        di['from_mode'] = get_decay_modes(values[0])
                        di['type'] = discrete_type[values[1]]
                        di['intensity'] = ufloat(*values[2:4])
                        if spectrum['type'] == 'ec/beta+':
                            di['positron_intensity'] = ufloat(*values[4:6])
                        elif spectrum['type'] == 'gamma':
                            if len(values) >= 6:
                                di['internal_pair'] = ufloat(*values[4:6])
                            if len(values) >= 8:
                                di['total_internal_conversion'] = ufloat(*values[6:8])
                            if len(values) == 12:
                                di['k_shell_conversion'] = ufloat(*values[8:10])
                                di['l_shell_conversion'] = ufloat(*values[10:12])
                        spectrum['discrete'].append(di)

                if not spectrum['continuous_flag'] == 'discrete':
                    # Read continuous spectrum
                    ci = {}
                    params, ci['probability'] = get_tab1_record(file_obj)
                    ci['type'] = get_decay_modes(params[0])

                    # Read covariance (Ek, Fk) table
                    LCOV = params[3]
                    if LCOV != 0:
                        items, values = get_list_record(file_obj)
                        ci['covariance_lb'] = items[3]
                        ci['covariance'] = zip(values[0::2], values[1::2])

                    spectrum['continuous'] = ci

                # Add spectrum to dictionary
                self.spectra[spectrum['type']] = spectrum

        else:
            items, values = get_list_record(file_obj)
            items, values = get_list_record(file_obj)
            self.nuclide['spin'] = items[0]
            self.nuclide['parity'] = items[1]
            self.half_life = ufloat(float('inf'), float('inf'))

    @property
    def decay_constant(self):
        if hasattr(self.half_life, 'n'):
            return log(2.)/self.half_life
        else:
            mu, sigma = self.half_life
            return ufloat(log(2.)/mu, log(2.)/mu**2*sigma)

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate radioactive decay data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : str or openmc.data.endf.Evaluation
            ENDF radioactive decay data evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.

        Returns
        -------
        openmc.data.Decay
            Radioactive decay data

        """
        return cls(ev_or_filename)
