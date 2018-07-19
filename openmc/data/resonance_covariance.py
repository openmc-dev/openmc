from collections import defaultdict, MutableSequence, Iterable
import warnings
import io

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

    def __init__(self, ranges):
        self.ranges = ranges

    def __iter__(self):
        for r in self.ranges:
            yield r

    @property
    def ranges(self):
        return self._ranges

    @ranges.setter
    def ranges(self, ranges):
        cv.check_type('resonance ranges', ranges, MutableSequence)
        self._ranges = cv.CheckedList(ResonanceCovarianceRange, 'resonance range',
                                      ranges)

    @classmethod
    def from_endf(cls, ev, resonances):
        """Generate resonance covariance data from an ENDF evaluation.

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        resonances : openmc.data.Resonance object
            Resonanance object generated from the same evaluation

        Returns
        -------
        openmc.data.ResonanceCovariances
            Resonance covariance data

        """
        file_obj = io.StringIO(ev.section[32, 151])

        # Determine whether discrete or continuous representation
        items = endf.get_head_record(file_obj)
        n_isotope = items[4] # Number of isotopes

        ranges = []
        for iso in range(n_isotope):
            items = endf.get_cont_record(file_obj)
            abundance = items[1]
            fission_widths = (items[3] == 1) # Flag for fission widths
            n_ranges = items[4] # number of resonance energy ranges

            for j in range(n_ranges):
                items = endf.get_cont_record(file_obj)
                unresolved_flag = items[2]  # 0: only scattering radius given
                                            # 1: resolved parameters given
                                            # 2: unresolved parameters given
                formalism = items[3]  # resonance formalism

                # Throw error for unsupported formalisms
                if formalism in [0, 7]:
                    raise NotImplementedError('LRF= ', formalism,
                                    'covariance not supported for this formalism')

                if unresolved_flag in (0, 1):
                    # resolved resonance region
                    file2params = resonances.ranges[j].parameters
                    erange = _FORMALISMS[formalism].from_endf(ev, file_obj,
                                                              items, file2params)
                elif unresolved_flag == 2:
                    warn_str = 'Unresolved resonance not supported.'\
                               'Covariance values for the unresolved region not imported.'
                    warnings.warn(warn_str)
                                 
                ranges.append(erange)

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
        Flag indicating the format of the covariance matrix within the ENDF file
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    formalism : str
        String descriptor of formalism
    """
    def __init__(self, energy_min, energy_max):
        self.energy_min = energy_min
        self.energy_max = energy_max
    
    def res_subset(self, parameter_str, bounds):
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
        parameters_subset : pandas.Dataframe
            Subset of parameters (maintains indexing of original)
        cov_subset : np.array
            Subset of covariance matrix (upper triangular)
        
        """
        parameters = self.parameters
        cov = self.covariance
        mpar = self.mpar
        mask1 = parameters[parameter_str] >= bounds[0]
        mask2 = parameters[parameter_str] <= bounds[1]
        mask = mask1 & mask2
        parameters_subset = parameters[mask] 
        indices = parameters_subset.index.values
        sub_cov_dim = len(indices)*mpar
        cov_subset_vals = []
        for index1 in indices:
            for i in range(mpar):
                for index2 in indices:
                    for j in range(mpar):
                        if index2*mpar+j >= index1*mpar+i:
                            cov_subset_vals.append(cov[index1*mpar+i, index2*mpar+j])
    
        cov_subset = np.zeros([sub_cov_dim, sub_cov_dim])
        tri_indices = np.triu_indices(sub_cov_dim)
        cov_subset[tri_indices] = cov_subset_vals
    
        self.parameters_subset = parameters_subset
        self.cov_subset = cov_subset

    def sample_resonance_parameters(self, n_samples, use_subset=False):
        """Return a IncidentNeutron object with n_samples of xs
    
        Parameters
        ----------
        n_samples : int
            The number of samples to produce
        use_subset : bool, optional
            Flag on whether to sample from an already produced subset
    
        Returns
        -------
        samples : list of openmc.data.ResonanceCovarianceRange objects 
            List of samples size [n_samples]
    
        """
        if not use_subset:
            parameters = self.parameters
            cov = self.covariance
        else:
            if self.parameters_subset is None:
                raise ValueError('No subset of resonances defined')
            parameters = self.parameters_subset
            cov = self.cov_subset

        nparams, params = parameters.shape
        cov = cov + cov.T - np.diag(cov.diagonal()) # symmetrizing covariance matrix
        covsize = cov.shape[0]
        formalism = self.formalism
        mpar = self.mpar
        samples = []
    
    
        # Handling MLBW sampling
        if formalism == 'mlbw' or formalism == 'slbw':
            if mpar == 3:
                param_list = ['energy', 'neutronWidth', 'captureWidth']
                mean_array = pd.DataFrame.as_matrix(parameters[param_list])
                spin = pd.DataFrame.as_matrix(parameters['J'])
                l_value = pd.DataFrame.as_matrix(parameters['L'])
                gf = pd.DataFrame.as_matrix(parameters['fissionWidth'])
                gx = pd.DataFrame.as_matrix(parameters['competitiveWidth'])
                mean = mean_array.flatten()
                for i in range(n_samples):
                    sample = np.random.multivariate_normal(mean, cov)
                    energy = sample[0::3]
                    gn = sample[1::3]
                    gg = sample[2::3]
                    gt = gn + gg + gf
                    records = []
                    for j, E in enumerate(energy):
                        records.append([energy[j], l_value[j], spin[j], gt[j], gn[j],
                                        gg[j], gf[j], gx[j]])
                    columns = ['energy', 'L', 'J', 'totalWidth', 'neutronWidth',
                           'captureWidth', 'fissionWidth', 'competitiveWidth']
                    sample_params = pd.DataFrame.from_records(records, columns=columns)
                    samples.append(sample_params)
    
            elif mpar == 4:
                param_list = ['energy', 'neutronWidth', 'captureWidth', 'fissionWidth']
                mean_array = pd.DataFrame.as_matrix(parameters[param_list])
                spin = pd.DataFrame.as_matrix(parameters['J'])
                l_value = pd.DataFrame.as_matrix(parameters['L'])
                gx = pd.DataFrame.as_matrix(parameters['competitiveWidth'])
                mean = mean_array.flatten()
                for i in range(n_samples):
                    sample = np.random.multivariate_normal(mean, cov)
                    energy = sample[0::4]
                    gn = sample[1::4]
                    gg = sample[2::4]
                    gf = sample[3::4]
                    gt = gn + gg + gf
                    records = []
                    for j, E in enumerate(energy):
                        records.append([energy[j], l_value[j], spin[j], gt[j], gn[j],
                                        gg[j], gf[j], gx[j]])
                    columns = ['energy', 'L', 'J', 'totalWidth', 'neutronWidth',
                           'captureWidth', 'fissionWidth', 'competitiveWidth']
                    sample_params = pd.DataFrame.from_records(records, columns=columns)
                    samples.append(sample_params)
    
            elif mpar == 5:
                param_list = ['energy', 'neutronWidth', 'captureWidth',
                              'fissionWidth', 'competitiveWidth']
                mean_array = pd.DataFrame.as_matrix(parameters[param_list])
                spin = pd.DataFrame.as_matrix(parameters['J'])
                l_value = pd.DataFrame.as_matrix(parameters['L'])
                mean = mean_array.flatten()
                for i in range(n_samples):
                    sample = np.random.multivariate_normal(mean, cov)
                    energy = sample[0::5]
                    gn = sample[1::5]
                    gg = sample[2::5]
                    gf = sample[3::5]
                    gx = sample[4::5]
                    gt = gn + gg + gf
                    records = []
                    for j, E in enumerate(energy):
                        records.append([energy[j], l_value[j], spin[j], gt[j], gn[j],
                                        gg[j], gf[j], gx[j]])
                    columns = ['energy', 'L', 'J', 'totalWidth', 'neutronWidth',
                           'captureWidth', 'fissionWidth', 'competitveWidth']
                    sample_params = pd.DataFrame.from_records(records, columns=columns)
                    samples.append(sample_params)

        # Handling RM Sampling
        if formalism == 'rm':
            if mpar == 3:
                param_list = ['energy', 'neutronWidth', 'captureWidth']
                mean_array = pd.DataFrame.as_matrix(parameters[param_list])
                spin = pd.DataFrame.as_matrix(parameters['J'])
                l_value = pd.DataFrame.as_matrix(parameters['L'])
                gfa = pd.DataFrame.as_matrix(parameters['fissionWidthA'])
                gfb = pd.DataFrame.as_matrix(parameters['fissionWidthB'])
                mean = mean_array.flatten()
                for i in range(n_samples):
                    sample = np.random.multivariate_normal(mean, cov)
                    energy = sample[0::3]
                    gn = sample[1::3]
                    gg = sample[2::3]
                    records = []
                    for j, E in enumerate(energy):
                        records.append([energy[j], l_value[j], spin[j], gn[j],
                                        gg[j], gfa[j], gfb[j]])
                    columns = ['energy', 'L', 'J', 'neutronWidth',
                               'captureWidth', 'fissionWidthA', 'fissionWidthB']
                    sample_params = pd.DataFrame.from_records(records, columns=columns)
                    samples.append(sample_params)
    
            elif mpar == 5:
                param_list = ['energy', 'neutronWidth', 'captureWidth',
                              'fissionWidthA', 'fissionWidthB']
                mean_array = pd.DataFrame.as_matrix(parameters[param_list])
                spin = pd.DataFrame.as_matrix(parameters['J'])
                l_value = pd.DataFrame.as_matrix(parameters['L'])
                mean = mean_array.flatten()
                for i in range(n_samples):
                    sample = np.random.multivariate_normal(mean, cov)
                    energy = sample[0::5]
                    gn = sample[1::5]
                    gg = sample[2::5]
                    gfa = sample[3::5]
                    gfb = sample[4::5]
                    records = []
                    for j, E in enumerate(energy):
                        records.append([energy[j], l_value[j], spin[j], gn[j],
                                        gg[j], gfa[j], gfb[j]])
                    columns = ['energy', 'L', 'J', 'neutronWidth',
                               'captureWidth', 'fissionWidthA', 'fissionWidthB']
                    sample_params = pd.DataFrame.from_records(records, columns=columns)
                    samples.append(sample_params)
    
        self.samples = samples

    def reconstruct(self, energies, resonances, sampleN):
        """Evaluate the cross section at specified energies for an already
        sampled set of resonance parameters. 
    
        Parameters
        ----------
        energies : float or Iterable of float
            Energies at which the cross section should be evaluated
        resonances : openmc.data.Resonance object
            Corresponding resonance range with File 2 data. Used for
            reconstruction method
        sampleN : int
            Index of sample of resonance parameters to be used
    
        Returns
        -------
        3-tuple of float or numpy.ndarray
            Elastic, capture, and fission cross sections at the specified
            energies
    
        """
        if self.samples[sampleN] is None:
            raise ValueError("Sample of resonance parameters has not been set.")
        sample_parameters = self.samples[sampleN]
        xs_array = resonances.reconstruct(energies, use_sample = True,
                                          sample_parameters = sample_parameters)
        return xs_array


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
    lcomp : int
        Flag indicating the format of the covariance matrix within the ENDF file
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    formalism : str
        String descriptor of formalism
    """

    def __init__(self, energy_min, energy_max):
        super().__init__(energy_min, energy_max)
        self.parameters = None
        self.covariance = None
        self.mpar = None
        self.lcomp = None
        self.formalism = 'mlbw'
            
    @classmethod
    def from_endf(cls, ev, file_obj, items, file2params):
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
        resonances : openmc.data.IncidentNeutron.Resonance object

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
        LCOMP = items[3] # Flag for compatibility 0, 1, 2 - 2 is compact form
        NLS = items[4]  # number of l-values

        # Build covariance matrix for General Resolved Resonance Formats
        if LCOMP == 1:
            items = endf.get_cont_record(file_obj)
            num_short_range = items[4] # Number of short range type resonance 
                                       # covariances
            num_long_range = items[5] # Number of long range type resonance 
                                      # covariances

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

            # Create pandas DataFrame with resonance data, currently 
            # redundant with data.IncidentNeutron.resonance
            columns = ['energy', 'J', 'totalWidth', 'neutronWidth',
                       'captureWidth', 'fissionWidth']
            parameters = pd.DataFrame.from_records(records, columns=columns)

            # Add parameters from File 2
            parameters = _add_file2_contributions(parameters, file2params)

            # Create instance of class
            mlbw = cls(energy_min, energy_max)
            mlbw.parameters = parameters
            mlbw.covariance = cov
            mlbw.mpar  = mpar
            mlbw.lcomp = LCOMP

            return mlbw

        elif LCOMP == 2: # Compact format - Resonances and individual
                         # uncertainties followed by compact correlations
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
                res_unc = values[i*12+6:i*12+12]
                # Delete 0 values (not provided, no fission width)
                # DAJ/DGT always zero, DGF sometimes none zero [1, 2, 5]
                res_unc_nonzero = []
                for j in range(6):
                    if j in [1, 2, 5] and res_unc[j] != 0.0 :
                        res_unc_nonzero.append(res_unc[j])
                    elif j in [0,3,4]:
                        res_unc_nonzero.append(res_unc[j])
                par_unc.extend(res_unc_nonzero)

            records = []
            for i, E in enumerate(energy):
                records.append([energy[i], spin[i], gt[i], gn[i],
                                gg[i], gf[i]])

            corr = endf.get_intg_record(file_obj)
            cov = np.diag(par_unc).dot(corr).dot(np.diag(par_unc))

            # Create pandas DataFrame with resonacne data
            columns = ['energy', 'J', 'totalWidth', 'neutronWidth',
                       'captureWidth', 'fissionWidth']
            parameters = pd.DataFrame.from_records(records, columns=columns)

            # Determine mpar (number of parameters for each resonance in
            # covariance matrix)
            nparams, params = parameters.shape
            covsize = cov.shape[0]
            mpar = int(covsize/nparams)
            
            # Add parameters from File 2
            parameters = _add_file2_contributions(parameters, file2params)

            # Create instance of MultiLevelBreitWignerCovariance
            mlbw = cls(energy_min, energy_max)
            mlbw.parameters = parameters
            mlbw.covariance = cov
            mlbw.mpar  = mpar
            mlbw.lcomp = LCOMP

            return mlbw

        elif LCOMP == 0 :
            cov = np.zeros([4, 4])
            records = []
            cov_index = 0
            for i in range(NLS):
                items, values = endf.get_list_record(file_obj)
                num_res = items[5]
                for j in range(num_res):
                    one_res = values[18*j:18*(j+1)]
                    res_values = one_res[:6]
                    cov_values = one_res[6:]

                    energy = res_values[0]
                    spin = res_values[1]
                    gt = res_values[2]
                    gn = res_values[3]
                    gg = res_values[4]
                    gf = res_values[5]
                    records.append([energy, spin, gt, gn, gg, gf])

                    # Populate the coviariance matrix for this resonance
                    # There are no covariances between resonances in LCOMP=0
                    cov[cov_index, cov_index] = cov_values[0]
                    cov[cov_index+1, cov_index+1 : cov_index+2] = cov_values[1:2]
                    cov[cov_index+1, cov_index+3] = cov_values[4]
                    cov[cov_index+2, cov_index+2] = cov_values[3]
                    cov[cov_index+2, cov_index+3] = cov_values[5]
                    cov[cov_index+3, cov_index+3] = cov_values[6]

                    cov_index += 4
                    if j < num_res-1: # Pad matrix for additional values
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
            parameters = _add_file2_contributions(parameters, file2params)

            # Create instance of class
            mlbw = cls(energy_min, energy_max)
            mlbw.parameters = parameters
            mlbw.covariance = cov
            mlbw.mpar  = mpar
            mlbw.lcomp = LCOMP

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
    lcomp : int
        Flag indicating the format of the covariance matrix within the ENDF file
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    formalism : str
        String descriptor of formalism
    """

    def __init__(self, energy_min, energy_max):
        super().__init__(energy_min, energy_max)
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
        Flag indicating the format of the covariance matrix within the ENDF file
    mpar : int
        Number of parameters in covariance matrix for each individual resonance
    formalism : str
        String descriptor of formalism
    """

    def __init__(self, energy_min, energy_max):
        super().__init__(energy_min, energy_max)
        self.parameters = None
        self.covariance = None
        self.formalism = 'rm'

    @classmethod
    def from_endf(cls, ev, file_obj, items, file2params):
        """Create Reich-Moore resonance covariance data from an ENDF evaluation.
        Includes the resonance parameters contained separately in File 32.

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
        resonances : Resonance object

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
        LCOMP = items[3]  # Flag for compatibility 0, 1, 2 - 2 is compact form
        NLS = items[4]  # Number of l-values

        
        # Build covariance matrix for General Resolved Resonance Formats
        if LCOMP == 1:
            items = endf.get_cont_record(file_obj)
            num_short_range = items[4] # Number of short range type resonance 
                                       # covariances
            num_long_range = items[5] # Number of long range type resonance
                                      # covariances
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
                cov = np.zeros([cov_dim,cov_dim])
                indices = np.triu_indices(cov_dim)
                cov[indices] = cov_values

            # Create pandas DataFrame with resonance data
            columns = ['energy', 'J', 'neutronWidth', 'captureWidth',
                       'fissionWidthA', 'fissionWidthB']
            parameters = pd.DataFrame.from_records(records, columns=columns)

            # Add parameters from File 2
            parameters = _add_file2_contributions(parameters, file2params)

            # Create instance of ReichMooreCovariance
            rmc = cls(energy_min, energy_max)
            rmc.parameters = parameters
            rmc.covariance = cov
            rmc.mpar  = mpar
            rmc.lcomp = LCOMP

            return rmc

        elif LCOMP == 2: # Compact format - Resonances and individual
                         # uncertainties followed by compact correlations
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
                res_unc = values[i*12+6:i*12+12]
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
            nparams,params = parameters.shape
            covsize = cov.shape[0]
            mpar = int(covsize/nparams)

            # Add parameters from File 2
            parameters = _add_file2_contributions(parameters, file2params)

            # Create instance of ReichMooreCovariance
            rmc = cls(energy_min, energy_max)
            rmc.parameters = parameters
            rmc.covariance = cov
            rmc.mpar  = mpar
            rmc.lcomp = LCOMP

            return rmc

_FORMALISMS = {
               0: ResonanceCovarianceRange,
               1: SingleLevelBreitWignerCovariance,
               2: MultiLevelBreitWignerCovariance,
               3: ReichMooreCovariance
               # 7: RMatrixLimitedCovariance
               }


