from collections import Callable
import sys
#from warnings import warn

import numpy as np
from numpy.polynomial.polynomial import Polynomial

from .function import Tabulated1D, Sum
from .endf_utils import read_float, read_CONT_line, identify_nuclide
import openmc.checkvalue as cv

if sys.version_info[0] >= 3:
    basestring = str


class FissionEnergyRelease(object):
    def __init__(self):
        self._fragments = None
        self._prompt_neutrons = None
        self._delayed_neutrons = None
        self._prompt_photons = None
        self._delayed_photons = None
        self._betas = None
        self._neutrinos = None
        self._form = None

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
        return Sum([self.fragments, self.prompt_neutrons, self.delayed_neutrons,
                    self.prompt_photons, self.delayed_photons, self.betas])

    @property
    def total(self):
        return Sum([self.fragments, self.prompt_neutrons, self.delayed_neutrons,
                    self.prompt_photons, self.delayed_photons, self.betas,
                    self.neutrinos])

    @property
    def prompt_q(self):
        return Sum([self.fragments, self.prompt_neutrons, self.prompt_photons,
                    lambda E: -E])

    @property
    def recoverable_q(self):
        return Sum([self.recoverable, lambda E: -E])

    @property
    def total_q(self):
        return Sum([self.total, lambda E: -E])
    
    @property
    def form(self):
        return self._form

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

    @form.setter
    def form(self, form):
        cv.check_value('format', form, ('Madland', 'Sher-Beck'))
        self._form = form

    @classmethod
    def from_endf(cls, filename, incident_neutron):
        """Generate fission energy release data from an ENDF file.

        Parameters
        ----------
        filename : str
            Name of the ENDF file containing fission energy release data

        incident_neutron : openmc.data.IncidentNeutron
            Corresponding incident neutron dataset

        Returns
        -------
        openmc.data.FissionEnergyRelease
            Fission energy release data

        """

        # Check to make sure this ENDF file matches the expected isomer.
        ident = identify_nuclide(filename)
        if ident['Z'] != incident_neutron.atomic_number:
            pass
        if ident['A'] != incident_neutron.mass_number:
            pass
        if ident['LISO'] != incident_neutron.metastable:
            pass

        # Extract the MF=1, MT=458 section.
        lines = []
        with open(filename, 'r') as fh:
            line = fh.readline()
            while line != '':
                if line[70:75] == ' 1458':
                    lines.append(line)
                line = fh.readline()

        # Read the number of coefficients in this LIST record.
        NPL = read_CONT_line(lines[1])[4]

        # Parse the ENDF LIST into an array.
        data = []
        for i in range(NPL):
            row, column = divmod(i, 6)
            data.append(read_float(lines[2 + row][11*column:11*(column+1)]))

        # Declare the coefficient names and the order they are given in.  The
        # LIST contains a value followed immediately by an uncertainty for each
        # of these components, times the polynomial order + 1.  If we only find
        # one value for each of these components, then we need to use the
        # Sher-Beck formula for energy dependence.  Otherwise, it is a
        # polynomial.
        labels = ('EFR', 'ENP', 'END', 'EGP', 'EGD', 'EB', 'ENU', 'ER', 'ET')

        # Associate each set of values and uncertainties with its label.
        value = dict()
        uncertainty = dict()
        for i in range(len(labels)):
            value[labels[i]] = data[2*i::18]
            uncertainty[labels[i]] = data[2*i + 1::18]

        # In ENDF/B-7.1, data for 2nd-order coefficients were mistakenly not
        # converted from MeV to eV.  Check for this error and fix it if present.
        n_coeffs = len(value['EFR'])
        if n_coeffs == 3:  # Only check 2nd-order data.
            # Check each energy component for the error.  If a 1 MeV neutron
            # causes a change of more than 100 MeV, we know something is wrong.
            error_present = False
            for coeffs in value.values():
                second_order = coeffs[2]
                if abs(second_order) * 1e12 > 1e8:
                    error_present = True
                    break

            # If we found the error, reduce all 2nd-order coeffs by 10**6.
            if error_present:
                for coeffs in value.values(): coeffs[2] *= 1e-6
                for coeffs in uncertainty.values(): coeffs[2] *= 1e-6

            # Perform the sanity check again... just in case.
            for coeffs in value.values():
                second_order = coeffs[2]
                if abs(second_order) * 1e12 > 1e8:
                    raise ValueError("Encountered a ludicrously large second-"
                                     "order polynomial coefficient.")

        # Convert eV to MeV.
        for coeffs in value.values():
            for i in range(len(coeffs)):
                coeffs[i] *= 10**(-6 + 6*i)
        for coeffs in uncertainty.values():
            for i in range(len(coeffs)):
                coeffs[i] *= 10**(-6 + 6*i)

        out = cls()
        if n_coeffs > 1:
            out.form = 'Madland'
            out.fragments = Polynomial(value['EFR'])
            out.prompt_neutrons = Polynomial(value['ENP'])
            out.delayed_neutrons = Polynomial(value['END'])
            out.prompt_photons = Polynomial(value['EGP'])
            out.delayed_photons = Polynomial(value['EGD'])
            out.betas = Polynomial(value['EB'])
            out.neutrinos = Polynomial(value['ENU'])
        else:
            out.form = 'Sher-Beck'
            raise NotImplemented

        return out

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

        obj = cls()
        if group.attrs['format'] == 'Madland':
            obj.fragments = Polynomial(group['fragments'].value)
            obj.prompt_neutrons = Polynomial(group['prompt_neutrons'].value)
            obj.delayed_neutrons = Polynomial(group['delayed_neutrons'].value)
            obj.prompt_photons = Polynomial(group['prompt_photons'].value)
            obj.delayed_photons = Polynomial(group['delayed_photons'].value)
            obj.betas = Polynomial(group['betas'].value)
            obj.neutrinos = Polynomial(group['neutrinos'].value)
        elif group.attrs['format'] == 'Sher-Beck':
            raise NotImplemented
        else:
            raise ValueError('Unrecognized energy release format')

        return obj

    def to_hdf5(self, group):
        """Write energy release data to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        if self.form == 'Madland':
            group.attrs['format'] = np.string_('Madland')
            group.create_dataset('fragments', data=self.fragments.coef)
            group.create_dataset('prompt_neutrons',
                                 data=self.prompt_neutrons.coef)
            group.create_dataset('delayed_neutrons',
                                 data=self.delayed_neutrons.coef)
            group.create_dataset('prompt_photons',
                                 data=self.prompt_photons.coef)
            group.create_dataset('delayed_photons',
                                 data=self.delayed_photons.coef)
            group.create_dataset('betas', data=self.betas.coef)
            group.create_dataset('neutrinos', data=self.neutrinos.coef)

            q_prompt = (self.fragments + self.prompt_neutrons +
                        self.prompt_photons + Polynomial((-1.0, 0.0)))
            group.create_dataset('q_prompt', data=q_prompt.coef)
            q_recoverable = (self.fragments + self.prompt_neutrons +
                             self.delayed_neutrons + self.prompt_photons +
                             self.delayed_photons + self.betas +
                             Polynomial((-1.0, 0.0)))
            group.create_dataset('q_recoverable', data=q_recoverable.coef)
            q_total = (self.fragments + self.prompt_neutrons +
                       self.delayed_neutrons + self.prompt_photons +
                       self.delayed_photons + self.betas + self.neutrinos +
                       Polynomial((-1.0, 0.0)))
            group.create_dataset('q_total', data=q_total.coef)
        elif self.form == 'Sher-Beck':
            group.attrs['format'] = np.string_('Sher-Beck')
            self.fragments.to_hdf5(group, 'fragments')
            self.prompt_neutrons.to_hdf5(group, 'prompt_neutrons')
            self.delayed_neutrons.to_hdf5(group, 'delayed_neutrons')
            self.prompt_photons.to_hdf5(group, 'prompt_photons')
            self.delayed_photons.to_hdf5(group, 'delayed_photons')
            self.betas.to_hdf5(group, 'betas')
            self.neutrinos.to_hdf5(group, 'neutrinos')
        else:
            raise ValueError('Unrecognized energy release format')
