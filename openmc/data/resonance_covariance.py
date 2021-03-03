from collections.abc import MutableSequence
import warnings
import io
import copy

import numpy as np
import pandas as pd

from . import endf
import openmc.checkvalue as cv
from .resonance import Resonances


def _add_file2_contributions(file32params, file2params):
    """Function for aiding in adding resonance parameters from File 2 that are
    not always present in File 32. Uses already imported resonance data.

    Paramaters
    ----------
    file32params : pandas.Dataframe
        Incomplete set of resonance parameters contained in File 32.
    file2params : pandas.Dataframe
        Resonance parameters from File 2. Ordered by energy.

    Returns
    -------
    parameters : pandas.Dataframe
        Complete set of parameters ordered by L-values and then energy

    """
    # Use l-values and competitiveWidth from File 2 data
    # Re-sort File 2 by energy to match File 32
    file2params = file2params.sort_values(by=['energy'])
    file2params.reset_index(drop=True, inplace=True)
    # Sort File 32 parameters by energy as well (maintaining index)
    file32params.sort_values(by=['energy'], inplace=True)
    # Add in values (.values converts to array first to ignore index)
    file32params['L'] = file2params['L'].values
    if 'competitiveWidth' in file2params.columns:
        file32params['competitiveWidth'] = file2params['competitiveWidth'].values
    # Resort to File 32 order (by L then by E) for use with covariance
    file32params.sort_index(inplace=True)
    return file32params


class ResonanceCovariances(Resonances):
    """Resolved resonance covariance data

    Parameters
    ----------
    ranges : list of openmc.data.ResonanceCovarianceRange
        Distinct energy ranges for resonance data

    Attributes
    ----------
    ranges : list of openmc.data.ResonanceCovarianceRange
        Distinct energy ranges for resonance data

    """

    @property
    def ranges(self):
        return self._ranges

    @ranges.setter
    def ranges(self, ranges):
        cv.check_type('resonance ranges', ranges, MutableSequence)
        self._ranges = cv.CheckedList(ResonanceCovarianceRange,
                                      'resonance range', ranges)

    @classmethod
    def from_endf(cls, ev, resonances):
        """Generate resonance covariance data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        resonances : openmc.data.Resonance object
            openmc.data.Resonanance object generated from the same evaluation
            used to import values not contained in File 32

        Returns
        -------
        openmc.data.ResonanceCovariances
            Resonance covariance data

        """
        file_obj = io.StringIO(ev.section[32, 151])

        # Determine whether discrete or continuous representation
        items = endf.get_head_record(file_obj)
        n_isotope = items[4]  # Number of isotopes

        ranges = []
        for iso in range(n_isotope):
            items = endf.get_cont_record(file_obj)
            abundance = items[1]
            fission_widths = (items[3] == 1)  # Flag for fission widths
            n_ranges = items[4]  # Number of resonance energy ranges

            for j in range(n_ranges):
                items = endf.get_cont_record(file_obj)
                # Unresolved flags - 0: only scattering radius given
                #                    1: resolved parameters given
                #                    2: unresolved parameters given
                unresolved_flag = items[2]
                formalism = items[3]  # resonance formalism

                # Throw error for unsupported formalisms
                if formalism in [0, 7]:
                    error = 'LRF='+str(formalism)+' covariance not supported '\
                            'for this formalism'
                    raise NotImplementedError(error)

                if unresolved_flag in (0, 1):
                    # Resolved resonance region
                    resonance = resonances.ranges[j]
                    erange = _FORMALISMS[formalism].from_endf(ev, file_obj,
                                                              items, resonance)
                    ranges.append(erange)

                elif unresolved_flag == 2:
                    warn = 'Unresolved resonance not supported. Covariance '\
                           'values for the unresolved region not imported.'
                    warnings.warn(warn)

        return cls(ranges)


class ResonanceCovarianceRange:
    """Resonace covariance range. Base class for different formalisms.

    Parameters
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV

    Attributes
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    parameters : pandas.DataFrame
        Resonance parameters
    covariance : numpy.array
        The covariance matrix contained within the ENDF evaluation
    lcomp : int
        Flag indicating format of the covariance matrix within the ENDF file
    file2res : openmc.data.ResonanceRange object
        Corresponding resonance range with File 2 data.
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    formalism : str
        String descriptor of formalism
    """
    def __init__(self, energy_min, energy_max):
        self.energy_min = energy_min
        self.energy_max = energy_max

    def subset(self, parameter_str, bounds):
        """Produce a subset of resonance parameters and the corresponding
        covariance matrix to an IncidentNeutron object.

        Parameters
        ----------
        parameter_str : str
            parameter to be discriminated
            (i.e. 'energy', 'captureWidth', 'fissionWidthA'...)
        bounds : np.array
            [low numerical bound, high numerical bound]

        Returns
        -------
        res_cov_range : openmc.data.ResonanceCovarianceRange
            ResonanceCovarianceRange object that contains a subset of the
            covariance matrix (upper triangular) as well as a subset parameters
            within self.file2params

        """
        # Copy range and prevent change of original
        res_cov_range = copy.deepcopy(self)

        parameters = self.file2res.parameters
        cov = res_cov_range.covariance
        mpar = res_cov_range.mpar
        # Create mask
        mask1 = parameters[parameter_str] >= bounds[0]
        mask2 = parameters[parameter_str] <= bounds[1]
        mask = mask1 & mask2
        res_cov_range.parameters = parameters[mask]
        indices = res_cov_range.parameters.index.values
        # Build subset of covariance
        sub_cov_dim = len(indices)*mpar
        cov_subset_vals = []
        for index1 in indices:
            for i in range(mpar):
                for index2 in indices:
                    for j in range(mpar):
                        if index2*mpar+j >= index1*mpar+i:
                            cov_subset_vals.append(cov[index1*mpar+i,
                                                   index2*mpar+j])

        cov_subset = np.zeros([sub_cov_dim, sub_cov_dim])
        tri_indices = np.triu_indices(sub_cov_dim)
        cov_subset[tri_indices] = cov_subset_vals

        res_cov_range.file2res.parameters = parameters[mask]
        res_cov_range.covariance = cov_subset
        return res_cov_range

    def sample(self, n_samples):
        """Sample resonance parameters based on the covariances provided
        within an ENDF evaluation.

        Parameters
        ----------
        n_samples : int
            The number of samples to produce

        Returns
        -------
        samples : list of openmc.data.ResonanceCovarianceRange objects
            List of samples size `n_samples`

        """
        warn_str = 'Sampling routine does not guarantee positive values for '\
                   'parameters. This can lead to undefined behavior in the '\
                   'reconstruction routine.'
        warnings.warn(warn_str)
        parameters = self.parameters
        cov = self.covariance

        # Symmetrizing covariance matrix
        cov = cov + cov.T - np.diag(cov.diagonal())
        formalism = self.formalism
        mpar = self.mpar
        samples = []

        # Handling MLBW/SLBW sampling
        if formalism == 'mlbw' or formalism == 'slbw':
            params = ['energy', 'neutronWidth', 'captureWidth', 'fissionWidth',
                      'competitiveWidth']
            param_list = params[:mpar]
            mean_array = parameters[param_list].values
            mean = mean_array.flatten()
            par_samples = np.random.multivariate_normal(mean, cov,
                                                        size=n_samples)
            spin = parameters['J'].values
            l_value = parameters['L'].values
            for sample in par_samples:
                energy = sample[0::mpar]
                gn = sample[1::mpar]
                gg = sample[2::mpar]
                gf = sample[3::mpar] if mpar > 3 else parameters['fissionWidth'].values
                gx = sample[4::mpar] if mpar > 4 else parameters['competitiveWidth'].values
                gt = gn + gg + gf + gx

                records = []
                for j, E in enumerate(energy):
                    records.append([energy[j], l_value[j], spin[j], gt[j],
                                    gn[j], gg[j], gf[j], gx[j]])
                columns = ['energy', 'L', 'J', 'totalWidth', 'neutronWidth',
                           'captureWidth', 'fissionWidth', 'competitiveWidth']
                sample_params = pd.DataFrame.from_records(records,
                                                          columns=columns)
                # Copy ResonanceRange object
                res_range = copy.copy(self.file2res)
                res_range.parameters = sample_params
                samples.append(res_range)

        # Handling RM sampling
        elif formalism == 'rm':
            params = ['energy', 'neutronWidth', 'captureWidth',
                      'fissionWidthA', 'fissionWidthB']
            param_list = params[:mpar]
            mean_array = parameters[param_list].values
            mean = mean_array.flatten()
            par_samples = np.random.multivariate_normal(mean, cov,
                                                        size=n_samples)
            spin = parameters['J'].values
            l_value = parameters['L'].values
            for sample in par_samples:
                energy = sample[0::mpar]
                gn = sample[1::mpar]
                gg = sample[2::mpar]
                gfa = sample[3::mpar] if mpar > 3 else parameters['fissionWidthA'].values
                gfb = sample[4::mpar] if mpar > 3 else parameters['fissionWidthB'].values

                records = []
                for j, E in enumerate(energy):
                    records.append([energy[j], l_value[j], spin[j], gn[j],
                                    gg[j], gfa[j], gfb[j]])
                columns = ['energy', 'L', 'J', 'neutronWidth',
                           'captureWidth', 'fissionWidthA', 'fissionWidthB']
                sample_params = pd.DataFrame.from_records(records,
                                                          columns=columns)
                # Copy ResonanceRange object
                res_range = copy.copy(self.file2res)
                res_range.parameters = sample_params
                samples.append(res_range)

        return samples


class MultiLevelBreitWignerCovariance(ResonanceCovarianceRange):
    """Multi-level Breit-Wigner resolved resonance formalism covariance data.
    Parameters
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV

    Attributes
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    parameters : pandas.DataFrame
        Resonance parameters
    covariance : numpy.array
        The covariance matrix contained within the ENDF evaluation
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    lcomp : int
        Flag indicating format of the covariance matrix within the ENDF file
    file2res : openmc.data.ResonanceRange object
        Corresponding resonance range with File 2 data.
    formalism : str
        String descriptor of formalism

    """

    def __init__(self, energy_min, energy_max, parameters, covariance, mpar,
                 lcomp, file2res):
        super().__init__(energy_min, energy_max)
        self.parameters = parameters
        self.covariance = covariance
        self.mpar = mpar
        self.lcomp = lcomp
        self.file2res = copy.copy(file2res)
        self.formalism = 'mlbw'

    @classmethod
    def from_endf(cls, ev, file_obj, items, resonance):
        """Create MLBW covariance data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        file_obj : file-like object
            ENDF file positioned at the second record of a resonance range
            subsection in MF=32, MT=151
        items : list
            Items from the CONT record at the start of the resonance range
            subsection
        resonance : openmc.data.ResonanceRange object
            Corresponding resonance range with File 2 data.

        Returns
        -------
        openmc.data.MultiLevelBreitWignerCovariance
            Multi-level Breit-Wigner resonance covariance parameters

        """

        # Read energy-dependent scattering radius if present
        energy_min, energy_max = items[0:2]
        nro, naps = items[4:6]
        if nro != 0:
            params, ape = endf.get_tab1_record(file_obj)

        # Other scatter radius parameters
        items = endf.get_cont_record(file_obj)
        target_spin = items[0]
        lcomp = items[3]  # Flag for compatibility 0, 1, 2 - 2 is compact form
        nls = items[4]  # number of l-values

        # Build covariance matrix for General Resolved Resonance Formats
        if lcomp == 1:
            items = endf.get_cont_record(file_obj)
            # Number of short range type resonance covariances
            num_short_range = items[4]
            # Number of long range type resonance covariances
            num_long_range = items[5]

            # Read resonance widths, J values, etc
            records = []
            for i in range(num_short_range):
                items, values = endf.get_list_record(file_obj)
                mpar = items[2]
                num_res = items[5]
                num_par_vals = num_res*6
                res_values = values[:num_par_vals]
                cov_values = values[num_par_vals:]

                energy = res_values[0::6]
                spin = res_values[1::6]
                gt = res_values[2::6]
                gn = res_values[3::6]
                gg = res_values[4::6]
                gf = res_values[5::6]

                for i, E in enumerate(energy):
                    records.append([energy[i], spin[i], gt[i], gn[i],
                                    gg[i], gf[i]])

                # Build the upper-triangular covariance matrix
                cov_dim = mpar*num_res
                cov = np.zeros([cov_dim, cov_dim])
                indices = np.triu_indices(cov_dim)
                cov[indices] = cov_values

        # Compact format - Resonances and individual uncertainties followed by
        # compact correlations
        elif lcomp == 2:
            items, values = endf.get_list_record(file_obj)
            mean = items
            num_res = items[5]
            energy = values[0::12]
            spin = values[1::12]
            gt = values[2::12]
            gn = values[3::12]
            gg = values[4::12]
            gf = values[5::12]
            par_unc = []
            for i in range(num_res):
                res_unc = values[i*12+6 : i*12+12]
                # Delete 0 values (not provided, no fission width)
                # DAJ/DGT always zero, DGF sometimes nonzero [1, 2, 5]
                res_unc_nonzero = []
                for j in range(6):
                    if j in [1, 2, 5] and res_unc[j] != 0.0:
                        res_unc_nonzero.append(res_unc[j])
                    elif j in [0, 3, 4]:
                        res_unc_nonzero.append(res_unc[j])
                par_unc.extend(res_unc_nonzero)

            records = []
            for i, E in enumerate(energy):
                records.append([energy[i], spin[i], gt[i], gn[i],
                                gg[i], gf[i]])

            corr = endf.get_intg_record(file_obj)
            cov = np.diag(par_unc).dot(corr).dot(np.diag(par_unc))

        # Compatible resolved resonance format
        elif lcomp == 0:
            cov = np.zeros([4, 4])
            records = []
            cov_index = 0
            for i in range(nls):
                items, values = endf.get_list_record(file_obj)
                num_res = items[5]
                for j in range(num_res):
                    one_res = values[18*j:18*(j+1)]
                    res_values = one_res[:6]
                    cov_values = one_res[6:]
                    records.append(list(res_values))

                    # Populate the coviariance matrix for this resonance
                    # There are no covariances between resonances in lcomp=0
                    cov[cov_index, cov_index] = cov_values[0]
                    cov[cov_index+1, cov_index+1 : cov_index+2] = cov_values[1:2]
                    cov[cov_index+1, cov_index+3] = cov_values[4]
                    cov[cov_index+2, cov_index+2] = cov_values[3]
                    cov[cov_index+2, cov_index+3] = cov_values[5]
                    cov[cov_index+3, cov_index+3] = cov_values[6]

                    cov_index += 4
                    if j < num_res-1:  # Pad matrix for additional values
                        cov = np.pad(cov, ((0, 4), (0, 4)), 'constant',
                                     constant_values=0)

        # Create pandas DataFrame with resonance data, currently
        # redundant with data.IncidentNeutron.resonance
        columns = ['energy', 'J', 'totalWidth', 'neutronWidth',
                   'captureWidth', 'fissionWidth']
        parameters = pd.DataFrame.from_records(records, columns=columns)
        # Determine mpar (number of parameters for each resonance in
        # covariance matrix)
        nparams, params = parameters.shape
        covsize = cov.shape[0]
        mpar = int(covsize/nparams)
        # Add parameters from File 2
        parameters = _add_file2_contributions(parameters,
                                              resonance.parameters)
        # Create instance of class
        mlbw = cls(energy_min, energy_max, parameters, cov, mpar, lcomp,
                   resonance)
        return mlbw


class SingleLevelBreitWignerCovariance(MultiLevelBreitWignerCovariance):
    """Single-level Breit-Wigner resolved resonance formalism covariance data.
    Single-level Breit-Wigner resolved resonance data is is identified by LRF=1
    in the ENDF-6 format.

    Parameters
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV

    Attributes
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    parameters : pandas.DataFrame
        Resonance parameters
    covariance : numpy.array
        The covariance matrix contained within the ENDF evaluation
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    formalism : str
        String descriptor of formalism
    lcomp : int
        Flag indicating format of the covariance matrix within the ENDF file
    file2res : openmc.data.ResonanceRange object
        Corresponding resonance range with File 2 data.
    """

    def __init__(self, energy_min, energy_max, parameters, covariance, mpar,
                 lcomp, file2res):
        super().__init__(energy_min, energy_max, parameters, covariance, mpar,
                         lcomp, file2res)
        self.formalism = 'slbw'


class ReichMooreCovariance(ResonanceCovarianceRange):
    """Reich-Moore resolved resonance formalism covariance data.

    Reich-Moore resolved resonance data is identified by LRF=3 in the ENDF-6
    format.

    Parameters
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV

    Attributes
    ----------
    energy_min : float
        Minimum energy of the resolved resonance range in eV
    energy_max : float
        Maximum energy of the resolved resonance range in eV
    parameters : pandas.DataFrame
        Resonance parameters
    covariance : numpy.array
        The covariance matrix contained within the ENDF evaluation
    lcomp : int
        Flag indicating format of the covariance matrix within the ENDF file
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    file2res : openmc.data.ResonanceRange object
        Corresponding resonance range with File 2 data.
    formalism : str
        String descriptor of formalism
    """

    def __init__(self, energy_min, energy_max, parameters, covariance, mpar,
                 lcomp, file2res):
        super().__init__(energy_min, energy_max)
        self.parameters = parameters
        self.covariance = covariance
        self.mpar = mpar
        self.lcomp = lcomp
        self.file2res = copy.copy(file2res)
        self.formalism = 'rm'

    @classmethod
    def from_endf(cls, ev, file_obj, items, resonance):
        """Create Reich-Moore resonance covariance data from an ENDF
        evaluation. Includes the resonance parameters contained separately in
        File 32.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        file_obj : file-like object
            ENDF file positioned at the second record of a resonance range
            subsection in MF=2, MT=151
        items : list
            Items from the CONT record at the start of the resonance range
            subsection
        resonance : openmc.data.Resonance object
            openmc.data.Resonanance object generated from the same evaluation
            used to import values not contained in File 32

        Returns
        -------
        openmc.data.ReichMooreCovariance
            Reich-Moore resonance covariance parameters

        """
        # Read energy-dependent scattering radius if present
        energy_min, energy_max = items[0:2]
        nro, naps = items[4:6]
        if nro != 0:
            params, ape = endf.get_tab1_record(file_obj)

        # Other scatter radius parameters
        items = endf.get_cont_record(file_obj)
        target_spin = items[0]
        lcomp = items[3]  # Flag for compatibility 0, 1, 2 - 2 is compact form
        nls = items[4]  # Number of l-values

        # Build covariance matrix for General Resolved Resonance Formats
        if lcomp == 1:
            items = endf.get_cont_record(file_obj)
            # Number of short range type resonance covariances
            num_short_range = items[4]
            # Number of long range type resonance covariances
            num_long_range = items[5]
            # Read resonance widths, J values, etc
            channel_radius = {}
            scattering_radius = {}
            records = []
            for i in range(num_short_range):
                items, values = endf.get_list_record(file_obj)
                mpar = items[2]
                num_res = items[5]
                num_par_vals = num_res*6
                res_values = values[:num_par_vals]
                cov_values = values[num_par_vals:]

                energy = res_values[0::6]
                spin = res_values[1::6]
                gn = res_values[2::6]
                gg = res_values[3::6]
                gfa = res_values[4::6]
                gfb = res_values[5::6]

                for i, E in enumerate(energy):
                    records.append([energy[i], spin[i], gn[i], gg[i],
                                    gfa[i], gfb[i]])

                # Build the upper-triangular covariance matrix
                cov_dim = mpar*num_res
                cov = np.zeros([cov_dim, cov_dim])
                indices = np.triu_indices(cov_dim)
                cov[indices] = cov_values

        # Compact format - Resonances and individual uncertainties followed by
        # compact correlations
        elif lcomp == 2:
            items, values = endf.get_list_record(file_obj)
            num_res = items[5]
            energy = values[0::12]
            spin = values[1::12]
            gn = values[2::12]
            gg = values[3::12]
            gfa = values[4::12]
            gfb = values[5::12]
            par_unc = []
            for i in range(num_res):
                res_unc = values[i*12+6 : i*12+12]
                # Delete 0 values (not provided in evaluation)
                res_unc = [x for x in res_unc if x != 0.0]
                par_unc.extend(res_unc)

            records = []
            for i, E in enumerate(energy):
                records.append([energy[i], spin[i], gn[i], gg[i],
                                gfa[i], gfb[i]])

            corr = endf.get_intg_record(file_obj)
            cov = np.diag(par_unc).dot(corr).dot(np.diag(par_unc))

        # Create pandas DataFrame with resonacne data
        columns = ['energy', 'J', 'neutronWidth', 'captureWidth',
                   'fissionWidthA', 'fissionWidthB']
        parameters = pd.DataFrame.from_records(records, columns=columns)

        # Determine mpar (number of parameters for each resonance in
        # covariance matrix)
        nparams, params = parameters.shape
        covsize = cov.shape[0]
        mpar = int(covsize/nparams)

        # Add parameters from File 2
        parameters = _add_file2_contributions(parameters,
                                              resonance.parameters)
        # Create instance of ReichMooreCovariance
        rmc = cls(energy_min, energy_max, parameters, cov, mpar, lcomp,
                  resonance)
        return rmc


_FORMALISMS = {
    0: ResonanceCovarianceRange,
    1: SingleLevelBreitWignerCovariance,
    2: MultiLevelBreitWignerCovariance,
    3: ReichMooreCovariance
    # 7: RMatrixLimitedCovariance
}
