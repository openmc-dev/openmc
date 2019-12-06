from collections.abc import Callable
from copy import deepcopy
from io import StringIO
import sys

import h5py
import numpy as np

from .data import EV_PER_MEV
from .endf import get_cont_record, get_list_record, get_tab1_record, Evaluation
from .function import Function1D, Tabulated1D, Polynomial, sum_functions
import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin


_NAMES = (
    'fragments', 'prompt_neutrons', 'delayed_neutrons',
    'prompt_photons', 'delayed_photons', 'betas',
    'neutrinos', 'recoverable', 'total'
)


class FissionEnergyRelease(EqualityMixin):
    """Energy relased by fission reactions.

    Energy is carried away from fission reactions by many different particles.
    The attributes of this class specify how much energy is released in the form
    of fission fragments, neutrons, photons, etc.  Each component is also (in
    general) a function of the incident neutron energy.

    Following a fission reaction, most of the energy release is carried by the
    daughter nuclei fragments.  These fragments accelerate apart from the
    Coulomb force on the time scale of ~10^-20 s [1].  Those fragments emit
    prompt neutrons between ~10^-18 and ~10^-13 s after scission (although some
    prompt neutrons may come directly from the scission point) [1].  Prompt
    photons follow with a time scale of ~10^-14 to ~10^-7 s [1].  The fission
    products then emit delayed neutrons with half lives between 0.1 and 100 s.
    The remaining fission energy comes from beta decays of the fission products
    which release beta particles, photons, and neutrinos (that escape the
    reactor and do not produce usable heat).

    Use the class methods to instantiate this class from an HDF5 or ENDF
    dataset.  The :meth:`FissionEnergyRelease.from_hdf5` method builds this
    class from the usual OpenMC HDF5 data files.
    :meth:`FissionEnergyRelease.from_endf` uses ENDF-formatted data.

    References
    ----------
    [1] D. G. Madland, "Total prompt energy release in the neutron-induced
    fission of ^235U, ^238U, and ^239Pu", Nuclear Physics A 772:113--137 (2006).
    <http://dx.doi.org/10.1016/j.nuclphysa.2006.03.013>

    Attributes
    ----------
    fragments : Callable
        Function that accepts incident neutron energy value(s) and returns the
        kinetic energy of the fission daughter nuclides (after prompt neutron
        emission).
    prompt_neutrons : Callable
        Function of energy that returns the kinetic energy of prompt fission
        neutrons.
    delayed_neutrons : Callable
        Function of energy that returns the kinetic energy of delayed neutrons
        emitted from fission products.
    prompt_photons : Callable
        Function of energy that returns the kinetic energy of prompt fission
        photons.
    delayed_photons : Callable
        Function of energy that returns the kinetic energy of delayed photons.
    betas : Callable
        Function of energy that returns the kinetic energy of delayed beta
        particles.
    neutrinos : Callable
        Function of energy that returns the kinetic energy of neutrinos.
    recoverable : Callable
        Function of energy that returns the kinetic energy of all products that
        can be absorbed in the reactor (all of the energy except for the
        neutrinos).
    total : Callable
        Function of energy that returns the kinetic energy of all products.
    q_prompt : Callable
        Function of energy that returns the prompt fission Q-value (fragments +
        prompt neutrons + prompt photons - incident neutron energy).
    q_recoverable : Callable
        Function of energy that returns the recoverable fission Q-value
        (total release - neutrinos - incident neutron energy).  This value is
        sometimes referred to as the pseudo-Q-value.
    q_total : Callable
        Function of energy that returns the total fission Q-value (total release
        - incident neutron energy).

    """
    def __init__(self, fragments, prompt_neutrons, delayed_neutrons,
                 prompt_photons, delayed_photons, betas, neutrinos):
        self.fragments = fragments
        self.prompt_neutrons = prompt_neutrons
        self.delayed_neutrons = delayed_neutrons
        self.prompt_photons = prompt_photons
        self.delayed_photons = delayed_photons
        self.betas = betas
        self.neutrinos = neutrinos

    @property
    def fragments(self):
        return self._fragments

    @property
    def prompt_neutrons(self):
        return self._prompt_neutrons

    @property
    def delayed_neutrons(self):
        return self._delayed_neutrons

    @property
    def prompt_photons(self):
        return self._prompt_photons

    @property
    def delayed_photons(self):
        return self._delayed_photons

    @property
    def betas(self):
        return self._betas

    @property
    def neutrinos(self):
        return self._neutrinos

    @property
    def recoverable(self):
        components = ['fragments', 'prompt_neutrons', 'delayed_neutrons',
                      'prompt_photons', 'delayed_photons', 'betas']
        return sum_functions(getattr(self, c) for c in components)

    @property
    def total(self):
        components = ['fragments', 'prompt_neutrons', 'delayed_neutrons',
                      'prompt_photons', 'delayed_photons', 'betas',
                      'neutrinos']
        return sum_functions(getattr(self, c) for c in components)

    @property
    def q_prompt(self):
        # Use a polynomial to subtract incident energy.
        funcs = [self.fragments, self.prompt_neutrons, self.prompt_photons,
                 Polynomial((0.0, -1.0))]
        return sum_functions(funcs)

    @property
    def q_recoverable(self):
        # Use a polynomial to subtract incident energy.
        return sum_functions([self.recoverable, Polynomial((0.0, -1.0))])

    @property
    def q_total(self):
        # Use a polynomial to subtract incident energy.
        return sum_functions([self.total, Polynomial((0.0, -1.0))])

    @fragments.setter
    def fragments(self, energy_release):
        cv.check_type('fragments', energy_release, Callable)
        self._fragments = energy_release

    @prompt_neutrons.setter
    def prompt_neutrons(self, energy_release):
        cv.check_type('prompt_neutrons', energy_release, Callable)
        self._prompt_neutrons = energy_release

    @delayed_neutrons.setter
    def delayed_neutrons(self, energy_release):
        cv.check_type('delayed_neutrons', energy_release, Callable)
        self._delayed_neutrons = energy_release

    @prompt_photons.setter
    def prompt_photons(self, energy_release):
        cv.check_type('prompt_photons', energy_release, Callable)
        self._prompt_photons = energy_release

    @delayed_photons.setter
    def delayed_photons(self, energy_release):
        cv.check_type('delayed_photons', energy_release, Callable)
        self._delayed_photons = energy_release

    @betas.setter
    def betas(self, energy_release):
        cv.check_type('betas', energy_release, Callable)
        self._betas = energy_release

    @neutrinos.setter
    def neutrinos(self, energy_release):
        cv.check_type('neutrinos', energy_release, Callable)
        self._neutrinos = energy_release

    @classmethod
    def from_endf(cls, ev, incident_neutron):
        """Generate fission energy release data from an ENDF file.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        incident_neutron : openmc.data.IncidentNeutron
            Corresponding incident neutron dataset

        Returns
        -------
        openmc.data.FissionEnergyRelease
            Fission energy release data

        """
        cv.check_type('evaluation', ev, Evaluation)

        # Check to make sure this ENDF file matches the expected isomer.
        if ev.target['atomic_number'] != incident_neutron.atomic_number:
            raise ValueError('The atomic number of the ENDF evaluation does '
                             'not match the given IncidentNeutron.')
        if ev.target['mass_number'] != incident_neutron.mass_number:
            raise ValueError('The atomic mass of the ENDF evaluation does '
                             'not match the given IncidentNeutron.')
        if ev.target['isomeric_state'] != incident_neutron.metastable:
            raise ValueError('The metastable state of the ENDF evaluation '
                             'does not match the given IncidentNeutron.')
        if not ev.target['fissionable']:
            raise ValueError('The ENDF evaluation is not fissionable.')

        if (1, 458) not in ev.section:
            raise ValueError('ENDF evaluation does not have MF=1, MT=458.')

        file_obj = StringIO(ev.section[1, 458])

        # Read first record and check whether any components appear as
        # tabulated functions
        items = get_cont_record(file_obj)
        lfc = items[3]
        nfc = items[5]

        # Parse the ENDF LIST into an array.
        items, data = get_list_record(file_obj)
        npoly = items[3]

        # Associate each set of values and uncertainties with its label.
        functions = {}
        for i, name in enumerate(_NAMES):
            coeffs = data[2*i::18]

            # Ignore recoverable and total since we recalculate those directly
            if name in ('recoverable', 'total'):
                continue

            # In ENDF/B-VII.1, data for 2nd-order coefficients were mistakenly
            # not converted from MeV to eV.  Check for this error and fix it if
            # present.
            if npoly == 2:  # Only check 2nd-order data.
                # If a 5 MeV neutron causes a change of more than 100 MeV, we
                # know something is wrong.
                second_order = coeffs[2]
                if abs(second_order) * (5e6)**2 > 1e8:
                    # If we found the error, reduce 2nd-order coeff by 10**6.
                    coeffs[2] /= EV_PER_MEV

            # If multiple coefficients were given, we can create the polynomial
            # and move on to the next component
            if npoly > 0:
                functions[name] = Polynomial(coeffs)
                continue

            # If a single coefficient was given, we need to use the Sher-Beck
            # formula for energy dependence
            zeroth_order = coeffs[0]
            if name in ('delayed_photons', 'betas'):
                func = Polynomial((zeroth_order, -0.075))
            elif name == 'neutrinos':
                func = Polynomial((zeroth_order, -0.105))
            elif name == 'prompt_neutrons':
                # Prompt neutrons require nu-data.  It is not clear from
                # ENDF-102 whether prompt or total nu value should be used, but
                # the delayed neutron fraction is so small that the difference
                # is negligible. MT=18 (n, fission) might not be available so
                # try MT=19 (n, f) as well.
                if 18 in incident_neutron and not incident_neutron[18].redundant:
                    nu = [p.yield_ for p in incident_neutron[18].products
                          if p.particle == 'neutron'
                          and p.emission_mode in ('prompt', 'total')]
                elif 19 in incident_neutron:
                    nu = [p.yield_ for p in incident_neutron[19].products
                          if p.particle == 'neutron'
                          and p.emission_mode in ('prompt', 'total')]
                else:
                    raise ValueError('IncidentNeutron data has no fission '
                                     'reaction.')
                if len(nu) == 0:
                    raise ValueError(
                        'Nu data is needed to compute fission energy '
                        'release with the Sher-Beck format.'
                    )
                if len(nu) > 1:
                    raise ValueError('Ambiguous prompt/total nu value.')

                nu = nu[0]
                if isinstance(nu, Tabulated1D):
                    # Evaluate Sher-Beck polynomial form at each tabulated value
                    func = deepcopy(nu)
                    func.y = (zeroth_order + 1.307*nu.x - 8.07e6*(nu.y - nu.y[0]))
                elif isinstance(nu, Polynomial):
                    # Combine polynomials
                    if len(nu) == 1:
                        func = Polynomial([zeroth_order, 1.307])
                    else:
                        func = Polynomial(
                            [zeroth_order, 1.307 - 8.07e6*nu.coef[1]]
                            + [-8.07e6*c for c in nu.coef[2:]])
            else:
                func = Polynomial(coeffs)

            functions[name] = func

        # Check for tabulated data
        if lfc == 1:
            for _ in range(nfc):
                # Get tabulated function
                items, eifc = get_tab1_record(file_obj)

                # Determine which component it is
                ifc = items[3]
                name = _NAMES[ifc - 1]

                # Replace value in dictionary
                functions[name] = eifc

        # Build the object
        return cls(**functions)

    @classmethod
    def from_hdf5(cls, group):
        """Generate fission energy release data from an HDF5 group.

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.FissionEnergyRelease
            Fission energy release data

        """

        fragments = Function1D.from_hdf5(group['fragments'])
        prompt_neutrons = Function1D.from_hdf5(group['prompt_neutrons'])
        delayed_neutrons = Function1D.from_hdf5(group['delayed_neutrons'])
        prompt_photons = Function1D.from_hdf5(group['prompt_photons'])
        delayed_photons = Function1D.from_hdf5(group['delayed_photons'])
        betas = Function1D.from_hdf5(group['betas'])
        neutrinos = Function1D.from_hdf5(group['neutrinos'])

        return cls(fragments, prompt_neutrons, delayed_neutrons, prompt_photons,
                   delayed_photons, betas, neutrinos)

    def to_hdf5(self, group):
        """Write energy release data to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        self.fragments.to_hdf5(group, 'fragments')
        self.prompt_neutrons.to_hdf5(group, 'prompt_neutrons')
        self.delayed_neutrons.to_hdf5(group, 'delayed_neutrons')
        self.prompt_photons.to_hdf5(group, 'prompt_photons')
        self.delayed_photons.to_hdf5(group, 'delayed_photons')
        self.betas.to_hdf5(group, 'betas')
        self.neutrinos.to_hdf5(group, 'neutrinos')
        self.q_prompt.to_hdf5(group, 'q_prompt')
        self.q_recoverable.to_hdf5(group, 'q_recoverable')
