from collections import Callable
from copy import deepcopy
import sys

import h5py
import numpy as np
from numpy.polynomial.polynomial import Polynomial

from .function import Tabulated1D, Sum
from .endf_utils import read_float, read_CONT_line, identify_nuclide
import openmc.checkvalue as cv

if sys.version_info[0] >= 3:
    basestring = str


def _extract_458_data(filename):
    """Read an ENDF file and extract the MF=1, MT=458 values.

    Parameters
    ----------
    filename : str
        Path to and ENDF file

    Returns
    -------
    value : dict of str to list of float
        Dictionary that gives lists of coefficients for each energy component.
        The keys are the 2-3 letter strings used in ENDF-102, e.g. 'EFR' and
        'ET'.  The list will have a length of 1 for Sher-Beck data, more for
        polynomial data.
    uncertainty : dict of str to list of float
        A dictionary with the same format as above.  This is probably a
        one-standard deviation value, but that is not specified explicitly in
        ENDF-102.  Also, some evaluations will give zero uncertainty.  Use with
        caution.

    """
    ident = identify_nuclide(filename)

    if not ident['LFI']:
        # This nuclide isn't fissionable.
        return None

    # Extract the MF=1, MT=458 section.
    lines = []
    with open(filename, 'r') as fh:
        line = fh.readline()
        while line != '':
            if line[70:75] == ' 1458':
                lines.append(line)
            line = fh.readline()

    if len(lines) == 0:
        # No 458 data here.
        return None

    # Read the number of coefficients in this LIST record.
    NPL = read_CONT_line(lines[1])[4]

    # Parse the ENDF LIST into an array.
    data = []
    for i in range(NPL):
        row, column = divmod(i, 6)
        data.append(read_float(lines[2 + row][11*column:11*(column+1)]))

    # Declare the coefficient names and the order they are given in.  The LIST
    # contains a value followed immediately by an uncertainty for each of these
    # components, times the polynomial order + 1.
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

    return value, uncertainty


def write_compact_458_library(endf_files, output_name=None, comment=None,
                              verbose=False):
    """Read ENDF files, strip the MF=1 MT=458 data and write to small HDF5.

    Parameters
    ----------
    endf_files : Collection of str
        Strings giving the paths to the ENDF files that will be parsed for data.
    output_name : str
        Name of the output HDF5 file.  Default is 'fission_Q_data.h5'.
    comment : str
        Comment to write in the output HDF5 file.  Defaults to no comment.
    verbose : bool
        If True, print the name of each isomer as it is read.  Defaults to
        False.

    """
    # Open the output file.
    if output_name is None: output_name = 'fission_Q_data.h5'
    out = h5py.File(output_name, 'w', libver='latest')

    # Write comments, if given.  This commented out comment is the one used for
    # the library distributed with OpenMC.
    #comment = ('This data is extracted from ENDF/B-VII.1 library.  Thanks '
    #           'evaluators, for all your hard work :)  Citation:  '
    #           'M. B. Chadwick, M. Herman, P. Oblozinsky, '
    #           'M. E. Dunn, Y. Danon, A. C. Kahler, D. L. Smith, '
    #           'B. Pritychenko, G. Arbanas, R. Arcilla, R. Brewer, '
    #           'D. A. Brown, R. Capote, A. D. Carlson, Y. S. Cho, H. Derrien, '
    #           'K. Guber, G. M. Hale, S. Hoblit, S. Holloway, T. D. Johnson, '
    #           'T. Kawano, B. C. Kiedrowski, H. Kim, S. Kunieda, '
    #           'N. M. Larson, L. Leal, J. P. Lestone, R. C. Little, '
    #           'E. A. McCutchan, R. E. MacFarlane, M. MacInnes, '
    #           'C. M. Mattoon, R. D. McKnight, S. F. Mughabghab, '
    #           'G. P. A. Nobre, G. Palmiotti, A. Palumbo, M. T. Pigni, '
    #           'V. G. Pronyaev, R. O. Sayer, A. A. Sonzogni, N. C. Summers, '
    #           'P. Talou, I. J. Thompson, A. Trkov, R. L. Vogt, '
    #           'S. C. van der Marck, A. Wallner, M. C. White, D. Wiarda, '
    #           'and P. G. Young. ENDF/B-VII.1 nuclear data for science and '
    #           'technology: Cross sections, covariances, fission product '
    #           'yields and decay data", Nuclear Data Sheets, '
    #           '112(12):2887-2996 (2011).')
    if comment is not None:
        out.attrs['comment'] = np.string_(comment)

    # Declare the order of the components.  Use fixed-length numpy strings
    # because they work well with h5py.
    labels = np.array(('EFR', 'ENP', 'END', 'EGP', 'EGD', 'EB', 'ENU', 'ER',
                       'ET'), dtype='S3')
    out.attrs['component order'] = labels

    # Iterate over the given files.
    if verbose: print('Reading ENDF files:')
    for fname in endf_files:
        if verbose: print(fname)

        ident = identify_nuclide(fname)

        # Skip non-fissionable nuclides.
        if not ident['LFI']: continue

        # Get the important bits.
        data = _extract_458_data(fname)
        if data is None: continue
        value, uncertainty = data

        # Make a group for this isomer.
        name = str(ident['Z']) + str(ident['A'])
        if ident['LISO'] != 0:
            name += '_m' + str(ident['LISO'])
        nuclide_group = out.create_group(name)

        # Write all the coefficients into one array.  The first dimension gives
        # the component (e.g. fragments or prompt neutrons); the second switches
        # between value and uncertainty; the third gives the polynomial order.
        n_coeffs = len(value['EFR'])
        data_out = np.zeros((len(labels), 2, n_coeffs))
        for i, label in enumerate(labels):
            data_out[i, 0, :] = value[label.decode()]
            data_out[i, 1, :] = uncertainty[label.decode()]
        nuclide_group.create_dataset('data', data=data_out)

    out.close()


class FissionEnergyRelease(object):
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
    :meth:`FissionEnergyRelease.from_compact_hdf5` uses a different HDF5 format
    that is meant to be compact and store the exact same data as the ENDF
    format.  Files with this format can be generated with the
    :func:`openmc.data.write_compact_458_library` function.

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
    form : str
        Format used to compute the energy-dependence of the data.  Either
        'Sher-Beck' or 'Madland'.

    """
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
    def q_prompt(self):
        return Sum([self.fragments, self.prompt_neutrons, self.prompt_photons,
                    lambda E: -E])

    @property
    def q_recoverable(self):
        return Sum([self.recoverable, lambda E: -E])

    @property
    def q_total(self):
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
    def _from_dictionary(cls, energy_release, incident_neutron):
        """Generate fission energy release data from a dictionary.

        Parameters
        ----------
        energy_release : dict of str to list of float
            Dictionary that gives lists of coefficients for each energy
            component.  The keys are the 2-3 letter strings used in ENDF-102,
            e.g. 'EFR' and 'ET'.  The list will have a length of 1 for Sher-Beck
            data, more for polynomial data.

        incident_neutron : openmc.data.IncidentNeutron
            Corresponding incident neutron dataset

        Returns
        -------
        openmc.data.FissionEnergyRelease
            Fission energy release data

        """
        out = cls()

        # How many coefficients are given for each component?  If we only find
        # one value for each, then we need to use the Sher-Beck formula for
        # energy dependence.  Otherwise, it is a polynomial.
        n_coeffs = len(energy_release['EFR'])
        if n_coeffs > 1:
            out.form = 'Madland'
            out.fragments = Polynomial(energy_release['EFR'])
            out.prompt_neutrons = Polynomial(energy_release['ENP'])
            out.delayed_neutrons = Polynomial(energy_release['END'])
            out.prompt_photons = Polynomial(energy_release['EGP'])
            out.delayed_photons = Polynomial(energy_release['EGD'])
            out.betas = Polynomial(energy_release['EB'])
            out.neutrinos = Polynomial(energy_release['ENU'])
        else:
            out.form = 'Sher-Beck'

            # EFR and ENP are energy independent.  Polynomial is used because it
            # has a __call__ attribute that handles Iterable inputs.  The
            # energy-dependence of END is unspecified in ENDF-102 so assume it
            # is independent.
            out.fragments = Polynomial((energy_release['EFR'][0]))
            out.prompt_photons = Polynomial((energy_release['EGP'][0]))
            out.delayed_neutrons = Polynomial((energy_release['END'][0]))

            # EDP, EB, and ENU are linear.
            out.delayed_photons = Polynomial((energy_release['EGD'][0], -0.075))
            out.betas = Polynomial((energy_release['EB'][0], -0.075))
            out.neutrinos = Polynomial((energy_release['ENU'][0], -0.105))

            # Prompt neutrons require nu-data.  It is not clear from ENDF-102
            # whether prompt or total nu value should be used, but the delayed
            # neutron fraction is so small that the difference is negligible.
            # MT=18 (n, fission) might not be available so try MT=19 (n, f) as
            # well.
            if 18 in incident_neutron.reactions:
                nu_prompt = [p for p in incident_neutron[18].products
                             if p.particle == 'neutron'
                             and p.emission_mode == 'prompt']
            elif 19 in incident_neutron.reactions:
                nu_prompt = [p for p in incident_neutron[19].products
                             if p.particle == 'neutron'
                             and p.emission_mode == 'prompt']
            else:
                raise ValueError('IncidentNeutron data has no fission '
                                  'reaction.')
            if len(nu_prompt) == 0:
                raise ValueError('Nu data is needed to compute fission energy '
                                 'release with the Sher-Beck format.')
            if len(nu_prompt) > 1:
                raise ValueError('Ambiguous prompt value.')
            if not isinstance(nu_prompt[0].yield_, Tabulated1D):
                raise TypeError('Sher-Beck fission energy release currently '
                                'only supports Tabulated1D nu data.')
            ENP = deepcopy(nu_prompt[0].yield_)
            ENP.y = (energy_release['ENP'] + 1.307 * ENP.x
                     - 8.07 * (ENP.y - ENP.y[0]))
            out.prompt_neutrons = ENP

        return out

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
            raise ValueError('The atomic number of the ENDF evaluation does '
                             'not match the given IncidentNeutron.')
        if ident['A'] != incident_neutron.mass_number:
            raise ValueError('The atomic mass of the ENDF evaluation does '
                             'not match the given IncidentNeutron.')
        if ident['LISO'] != incident_neutron.metastable:
            raise ValueError('The metastable state of the ENDF evaluation does '
                             'not match the given IncidentNeutron.')
        if not ident['LFI']:
            raise ValueError('The ENDF evaluation is not fissionable.')

        # Read the 458 data from the ENDF file.
        value, uncertainty = _extract_458_data(filename)

        # Build the object.
        return cls._from_dictionary(value, incident_neutron)

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

        obj.fragments = Polynomial(group['fragments'].value)
        obj.delayed_neutrons = Polynomial(group['delayed_neutrons'].value)
        obj.prompt_photons = Polynomial(group['prompt_photons'].value)
        obj.delayed_photons = Polynomial(group['delayed_photons'].value)
        obj.betas = Polynomial(group['betas'].value)
        obj.neutrinos = Polynomial(group['neutrinos'].value)

        if group.attrs['format'].decode() == 'Madland':
            obj.form = 'Madland'
            obj.prompt_neutrons = Polynomial(group['prompt_neutrons'].value)
        elif group.attrs['format'].decode() == 'Sher-Beck':
            obj.form = 'Sher-Beck'
            obj.prompt_neutrons = Tabulated1D.from_hdf5(
                                                       group['prompt_neutrons'])
        else:
            raise ValueError('Unrecognized energy release format')

        return obj

    @classmethod
    def from_compact_hdf5(cls, fname, incident_neutron):
        """Generate fission energy release data from a small HDF5 library.

        Parameters
        ----------
        fname : str
            Path to an HDF5 file containing fission energy release data.  This
            file should have been generated form the
            :func:`openmc.data.write_compact_458_library` function.

        incident_neutron : openmc.data.IncidentNeutron
            Corresponding incident neutron dataset

        Returns
        -------
        openmc.data.FissionEnergyRelease or None
            Fission energy release data for the given nuclide if it is present
            in the data file

        """

        fin = h5py.File(fname, 'r')

        components = [s.decode() for s in fin.attrs['component order']]

        nuclide_name = str(incident_neutron.atomic_number)
        nuclide_name += str(incident_neutron.mass_number)
        if incident_neutron.metastable != 0:
            nuclide_name += '_m' + str(incident_neutron.metastable)

        if nuclide_name not in fin: return None

        data = {c : fin[nuclide_name + '/data'][i, 0, :]
                for i, c in enumerate(components)}

        return cls._from_dictionary(data, incident_neutron)

    def to_hdf5(self, group):
        """Write energy release data to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.create_dataset('fragments', data=self.fragments.coef)
        group.create_dataset('delayed_neutrons',
                             data=self.delayed_neutrons.coef)
        group.create_dataset('prompt_photons',
                             data=self.prompt_photons.coef)
        group.create_dataset('delayed_photons',
                             data=self.delayed_photons.coef)
        group.create_dataset('betas', data=self.betas.coef)
        group.create_dataset('neutrinos', data=self.neutrinos.coef)

        if self.form == 'Madland':
            group.attrs['format'] = np.string_('Madland')
            group.create_dataset('prompt_neutrons',
                                 data=self.prompt_neutrons.coef)

            q_prompt = (self.fragments + self.prompt_neutrons +
                        self.prompt_photons + Polynomial((0.0, -1.0)))
            group.create_dataset('q_prompt', data=q_prompt.coef)
            q_recoverable = (self.fragments + self.prompt_neutrons +
                             self.delayed_neutrons + self.prompt_photons +
                             self.delayed_photons + self.betas +
                             Polynomial((0.0, -1.0)))
            group.create_dataset('q_recoverable', data=q_recoverable.coef)
        elif self.form == 'Sher-Beck':
            group.attrs['format'] = np.string_('Sher-Beck')
            self.prompt_neutrons.to_hdf5(group, 'prompt_neutrons')

            q_prompt = deepcopy(self.prompt_neutrons)
            q_prompt.y += self.fragments(q_prompt.x)
            q_prompt.y += self.prompt_photons(q_prompt.x)
            q_prompt.to_hdf5(group, 'q_prompt')
            q_recoverable = q_prompt
            q_recoverable.y += self.delayed_neutrons(q_recoverable.x)
            q_recoverable.y += self.delayed_photons(q_recoverable.x)
            q_recoverable.y += self.betas(q_recoverable.x)
            q_recoverable.to_hdf5(group, 'q_recoverable')
        else:
            raise ValueError('Unrecognized energy release format')
