from abc import ABC
from collections import OrderedDict
import copy
import sys

import h5py
import numpy as np
import scipy.sparse as sps

import openmc
import openmc.kinetics
import openmc.mgxs


class State(ABC):
    """State to store all the variables that describe a specific state of the system.

    Attributes
    ----------
    amplitude_mesh : openmc.RegularMesh
        Mesh by which shape is computed on.
    shape_mesh : openmc.RegularMesh
        Mesh by which power is reconstructed on.
    unity_mesh : openmc.RegularMesh
        Mesh with one cell convering the entire geometry..
    tally_mesh : openmc.RegularMesh
        Mesh to tally currents
    one_group : openmc.mgxs.groups.EnergyGroups
        EnergyGroups which specifies the a one-energy-group structure.
    energy_groups : openmc.mgxs.groups.EnergyGroups
        EnergyGroups which specifies the energy groups structure.
    fine_groups : openmc.mgxs.groups.EnergyGroups
        EnergyGroups used to tally the transport cross section that will be
        condensed to get the diffusion coefficients in the coarse group
        structure.
    tally_groups : openmc.mgxs.groups.EnergyGroups
        EnergyGroups used in all tallies except the diffusion coefficients
        (see fine_groups).
    shape : np.ndarray
        Numpy array used to store the shape function.
    amplitude : np.ndarray
        Numpy array used to store the amplitude.
    adjoint_flux : np.ndarray
        Numpy array used to store the adjoint flux.
    precursors : np.ndarray
        Numpy array used to store the precursor concentrations.
    mgxs_lib : OrderedDict of OrderedDict of openmc.tallies
        Dict of Dict of tallies. The first Dict is indexed by time point
        and the second Dict is indexed by rxn type.
    k_crit : float
        The initial eigenvalue.
    chi_delayed_by_delayed_group : bool
        Whether to use delayed groups in representing chi-delayed.
    chi_delayed_by_mesh : bool
        Whether to use a mesh in representing chi-delayed.
    use_agd : bool
        Whether to use artificial grid diffusion
    use_pcmfd : bool
        Whether to use p-CMFD
    num_delayed_groups : int
        The number of delayed neutron precursor groups.
    time_point : str
        The time point of this state.
    clock : openmc.kinetics.Clock
        A clock object to indicate the current simulation time.
    core_volume : float
        The core volume used to normalize the initial power.
    log_file : str
        Log file name (including directory prefix).
    multi_group : bool
        Whether the OpenMC run is multi-group or continuous-energy.

    """

    def __init__(self, states):

        # Initialize Solver class attributes
        self._amplitude_mesh = None
        self._shape_mesh = None
        self._unity_mesh = None
        self._tally_mesh = None

        self._one_group = None
        self._energy_groups = None
        self._fine_groups = None
        self._tally_groups = None
        self._num_delayed_groups = 6

        self._precursors = None
        self._adjoint_flux = None
        self._amplitude = None
        self._k_crit = 1.0

        self._time_point = None
        self._clock = None
        self._core_volume = 1.

        self._log_file = None
        self._multi_group = True
        self._use_agd = False
        self._use_pcmfd = False
        self.states = states

    @property
    def states(self):
        return self._states

    @property
    def amplitude_mesh(self):
        return self._amplitude_mesh

    @property
    def shape_mesh(self):
        return self._shape_mesh

    @property
    def unity_mesh(self):
        return self._unity_mesh

    @property
    def tally_mesh(self):
        return self._tally_mesh

    @property
    def one_group(self):
        return self._one_group

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def fine_groups(self):
        return self._fine_groups

    @property
    def tally_groups(self):
        return self._tally_groups

    @property
    def method(self):
        return self._method

    @property
    def amplitude(self):
        self._amplitude.shape = (self.amplitude_nxyz, self.ng)
        return self._amplitude

    @property
    def adjoint_flux(self):
        self._adjoint_flux.shape = (self.amplitude_nxyz, self.ng)
        return self._adjoint_flux

    @property
    def precursors(self):
        self._precursors.shape = (self.shape_nxyz, self.nd)
        return self._precursors

    @property
    def k_crit(self):
        return self._k_crit

    @property
    def num_delayed_groups(self):
        return self._num_delayed_groups

    @property
    def time_point(self):
        return self._time_point

    @property
    def clock(self):
        return self._clock

    @property
    def core_volume(self):
        return self._core_volume

    @property
    def log_file(self):
        return self._log_file

    @property
    def multi_group(self):
        return self._multi_group

    @states.setter
    def states(self, states):
        self._states = states

    @shape_mesh.setter
    def shape_mesh(self, mesh):
        self._shape_mesh = mesh

    @amplitude_mesh.setter
    def amplitude_mesh(self, mesh):
        self._amplitude_mesh = mesh

    @unity_mesh.setter
    def unity_mesh(self, mesh):
        self._unity_mesh = mesh

    @tally_mesh.setter
    def tally_mesh(self, mesh):
        self._tally_mesh = mesh

    @one_group.setter
    def one_group(self, one_group):
        self._one_group = one_group

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        self._energy_groups = energy_groups

    @fine_groups.setter
    def fine_groups(self, fine_groups):
        self._fine_groups = fine_groups

    @tally_groups.setter
    def tally_groups(self, tally_groups):
        self._tally_groups = tally_groups

    @method.setter
    def method(self, method):
        self._method = method

    @amplitude.setter
    def amplitude(self, amp):
        self._amplitude = copy.deepcopy(amp)
        self._amplitude.shape = (self.amplitude_nxyz, self.ng)

    @adjoint_flux.setter
    def adjoint_flux(self, adjoint_flux):
        self._adjoint_flux = copy.deepcopy(adjoint_flux)
        self._adjoint_flux.shape = (self.amplitude_nxyz, self.ng)

    @precursors.setter
    def precursors(self, precursors):
        self._precursors = copy.deepcopy(precursors)
        self._precursors.shape = (self.shape_nxyz, self.nd)

    @k_crit.setter
    def k_crit(self, k_crit):
        self._k_crit = k_crit

    @num_delayed_groups.setter
    def num_delayed_groups(self, num_delayed_groups):
        self._num_delayed_groups = num_delayed_groups

    @time_point.setter
    def time_point(self, time_point):
        self._time_point = time_point

    @clock.setter
    def clock(self, clock):
        self._clock = clock

    @core_volume.setter
    def core_volume(self, core_volume):
        self._core_volume = core_volume

    @log_file.setter
    def log_file(self, log_file):
        self._log_file = log_file

    @multi_group.setter
    def multi_group(self, multi_group):
        self._multi_group = multi_group

    @property
    def shape_dimension(self):
        return tuple(self.shape_mesh.dimension[::-1])

    @property
    def tally_dimension(self):
        return tuple(self.tally_mesh.dimension[::-1])

    @property
    def amplitude_dimension(self):
        return tuple(self.amplitude_mesh.dimension[::-1])

    @property
    def shape_zyxg(self):
        return self.shape_dimension + (self.ng,)

    @property
    def tally_zyxg(self):
        return self.tally_dimension + (self.ng,)

    @property
    def amplitude_zyxg(self):
        return self.amplitude_dimension + (self.ng,)

    @property
    def shape_zyxgg(self):
        return self.shape_dimension + (self.ng, self.ng)

    @property
    def amplitude_zyxgg(self):
        return self.amplitude_dimension + (self.ng, self.ng)

    @property
    def ng(self):
        return self.energy_groups.num_groups

    @property
    def nd(self):
        return self.num_delayed_groups

    @property
    def dt_inner(self):
        return self.clock.dt_inner

    @property
    def dt_outer(self):
        return self.clock.dt_outer

    @property
    def shape_nxyz(self):
        return np.prod(self.shape_dimension)

    @property
    def tally_nxyz(self):
        return np.prod(self.tally_dimension)

    @property
    def amplitude_nxyz(self):
        return np.prod(self.amplitude_dimension)

    @property
    def shape_dxyz(self):
        if self.shape_mesh.width:
            width = list(self.shape_mesh.width)
        else:
            width = list([(i-j)/k for i,j,k in 
                zip(self.shape_mesh.upper_right, 
                    self.shape_mesh.lower_left, 
                    self.shape_mesh.dimension)])
        dim = list(self.shape_dimension)[::-1]
        width = [i/j for i,j in zip(width,dim)]
        return np.prod(width)

    @property
    def amplitude_dxyz(self):
        if self.amplitude_mesh.width:
            width = list(self.amplitude_mesh.width)
        else:
            width = list([(i-j)/k for i,j,k in 
                zip(self.amplitude_mesh.upper_right, 
                    self.amplitude_mesh.lower_left, 
                    self.amplitude_mesh.dimension)])
        dim = list(self.amplitude_dimension)[::-1]
        width = [i/j for i,j in zip(width,dim)]
        return np.prod(width)

    @property
    def power(self):
        return (self.kappa_fission * self.flux).sum(axis=1)

    @property
    def core_power_density(self):
        return self.power.sum() / self.core_volume

    @property
    def flux(self):
        amp = self.amplitude
        amp = openmc.kinetics.map_array(amp, self.amplitude_zyxg, self.shape_zyxg, True)
        amp.shape = (self.shape_nxyz, self.ng)
        return amp * self.shape

    @property
    def delayed_production(self):
        chi_delayed = np.repeat(self.chi_delayed, self.ng)
        chi_delayed.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)
        delayed_nu_fission = np.tile(self.delayed_nu_fission, self.ng)
        delayed_nu_fission.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)
        return (chi_delayed * delayed_nu_fission)

    @property
    def prompt_production(self):
        chi_prompt = np.repeat(self.chi_prompt, self.ng)
        chi_prompt.shape = (self.shape_nxyz, self.ng, self.ng)
        prompt_nu_fission = np.tile(self.prompt_nu_fission, self.ng)
        prompt_nu_fission.shape = (self.shape_nxyz, self.ng, self.ng)
        return (chi_prompt * prompt_nu_fission)

    @property
    def delayed_production_matrix(self):
        delayed_production = self.delayed_production.sum(axis=1)
        shape = np.tile(self.shape, self.ng).reshape((self.shape_nxyz, self.ng, self.ng))
        delayed_production *= shape / self.k_crit

        # Condense down to the coarse mesh
        delayed_production = openmc.kinetics.map_array(delayed_production, self.shape_zyxgg, self.amplitude_zyxgg, False)
        delayed_production.shape = (self.amplitude_nxyz, self.ng, self.ng)

        return openmc.kinetics.block_diag(delayed_production)

    @property
    def production_matrix(self):
        return self.prompt_production_matrix + self.delayed_production_matrix

    @property
    def prompt_production_matrix(self):
        prompt_production = self.prompt_production / self.k_crit

        shape = np.tile(self.shape, self.ng).reshape((self.shape_nxyz, self.ng, self.ng))
        prompt_production = prompt_production * shape
        prompt_production = openmc.kinetics.map_array(prompt_production, self.shape_zyxgg, self.amplitude_zyxgg, False)
        prompt_production.shape = (self.amplitude_nxyz, self.ng, self.ng)

        return openmc.kinetics.block_diag(prompt_production)

    @property
    def destruction_matrix(self):

        linear, non_linear = self.coupling_matrix
        shape              = np.tile(self.shape, self.ng).reshape((self.shape_nxyz, self.ng, self.ng))
        inscatter          = self.inscatter * shape
        absorb_outscat     = self.outscatter + self.absorption
        absorb_outscat     = absorb_outscat * self.shape

        # Condense down to the coarse mesh
        inscatter       = openmc.kinetics.map_array(inscatter, self.shape_zyxgg, self.amplitude_zyxgg, False)
        absorb_outscat  = openmc.kinetics.map_array(absorb_outscat, self.shape_zyxg, self.amplitude_zyxg, False)
        inscatter.shape = (self.amplitude_nxyz, self.ng, self.ng)

        total = sps.diags(absorb_outscat.flatten()) - openmc.kinetics.block_diag(inscatter)

        return total + linear + non_linear

    @property
    def coupling_matrix(self):

        diags, dc_linear_data, dc_nonlinear_data = self.coupling_terms

        # Form a matrix of the surface diffusion coefficients corrections
        dc_linear_matrix    = sps.diags(dc_linear_data   , diags)
        dc_nonlinear_matrix = sps.diags(dc_nonlinear_data, diags)

        return dc_linear_matrix, dc_nonlinear_matrix

    @property
    def pnl(self):
        inv_vel = self.inverse_velocity * self.flux
        inv_vel = openmc.kinetics.map_array(inv_vel, self.shape_zyxg, self.amplitude_zyxg, False)
        inv_vel.shape = (self.amplitude_nxyz, self.ng)
        inv_vel *= self.adjoint_flux
        production = self.production_matrix * self.amplitude.flatten()
        production = production.flatten() * self.adjoint_flux.flatten()
        return inv_vel.sum() / production.sum()

    @property
    def reactivity(self):
        production = self.production_matrix * self.amplitude.flatten()
        destruction = self.destruction_matrix * self.amplitude.flatten()
        balance = production - destruction

        balance = balance * self.adjoint_flux.flatten()
        production = production * self.adjoint_flux.flatten()
        return balance.sum() / production.sum()

    @property
    def beta_eff(self):
        flux = np.tile(self.flux, self.nd).reshape((self.shape_nxyz, self.nd, self.ng))
        adjoint_flux = np.tile(self.adjoint_flux, self.nd).reshape((self.amplitude_nxyz, self.nd, self.ng))

        delayed_production = self.delayed_production
        delayed_production.shape = (self.shape_nxyz * self.nd, self.ng, self.ng)
        delayed_production = openmc.kinetics.block_diag(delayed_production)
        delayed_production /= self.k_crit

        delayed_production *= flux.flatten()
        old_shape = (self.shape_nxyz, self.nd, self.ng)
        new_shape = (self.amplitude_nxyz, self.nd, self.ng)
        delayed_production = openmc.kinetics.map_array(delayed_production, old_shape, new_shape, False)
        delayed_production = delayed_production.flatten() * adjoint_flux.flatten()
        delayed_production.shape = (self.amplitude_nxyz, self.nd, self.ng)
        delayed_production = delayed_production.sum(axis=(0,2))

        production = self.production_matrix * self.amplitude.flatten()
        production = production.flatten() * self.adjoint_flux.flatten()
        production.shape = (self.amplitude_nxyz, self.ng)
        production = production.sum(axis=(0,1))
        production = np.repeat(production, self.nd)

        return (delayed_production / production)

    def beta(self, integrated=False):
        flux = np.tile(self.flux, self.nd).reshape((self.shape_nxyz, self.nd, self.ng))

        delayed_production = self.delayed_production
        delayed_production.shape = (self.shape_nxyz * self.nd, self.ng, self.ng)
        delayed_production = openmc.kinetics.block_diag(delayed_production)
        delayed_production /= self.k_crit

        delayed_production *= flux.flatten()
        old_shape = (self.shape_nxyz, self.nd, self.ng)
        new_shape = (self.amplitude_nxyz, self.nd, self.ng)
        delayed_production = openmc.kinetics.map_array(delayed_production, old_shape, new_shape, False)
        delayed_production.shape = (self.amplitude_nxyz, self.nd, self.ng)
        if integrated:
            delayed_production = delayed_production.sum(axis=(0,2))
        else:
            delayed_production = delayed_production.sum(axis=(2,))

        production = self.production_matrix * self.amplitude.flatten()
        production.shape = (self.amplitude_nxyz, self.ng)
        if integrated:
            production = production.sum(axis=(0,1))
            production = np.repeat(production, self.nd)
        else:
            production = production.sum(axis=(1,))
            production = np.repeat(production, self.nd)
            production.shape = (self.amplitude_nxyz, self.nd)

        return delayed_production / production

    def beta_shape(self, integrated=False):
        flux = np.tile(self.flux, self.nd).reshape((self.shape_nxyz, self.nd, self.ng))

        delayed_production = self.delayed_production
        delayed_production.shape = (self.shape_nxyz * self.nd, self.ng, self.ng)
        delayed_production = openmc.kinetics.block_diag(delayed_production)
        delayed_production /= self.k_crit

        delayed_production *= flux.flatten()
        delayed_production.shape = (self.shape_nxyz, self.nd, self.ng)
        if integrated:
            delayed_production = delayed_production.sum(axis=(0,2))
        else:
            delayed_production = delayed_production.sum(axis=(2,))

        prompt_production = self.prompt_production / self.k_crit
        shape = np.tile(self.shape, self.ng).reshape((self.shape_nxyz, self.ng, self.ng))
        prompt_production = prompt_production * shape
        pp_matrix = openmc.kinetics.block_diag(prompt_production)

        del_prod = self.delayed_production.sum(axis=1)
        shape = np.tile(self.shape, self.ng).reshape((self.shape_nxyz, self.ng, self.ng))
        del_prod *= shape / self.k_crit
        del_prod.shape = (self.shape_nxyz, self.ng, self.ng)
        dp_matrix = openmc.kinetics.block_diag(del_prod)

        prod_matrix = pp_matrix + dp_matrix
        amp = self.amplitude
        amp = openmc.kinetics.map_array(amp, self.amplitude_zyxg, self.shape_zyxg, True)
        amp.shape = (self.shape_nxyz, self.ng)

        production = prod_matrix * amp.flatten()
        production.shape = (self.shape_nxyz, self.ng)
        if integrated:
            production = production.sum(axis=(0,1))
            production = np.repeat(production, self.nd)
        else:
            production = production.sum(axis=(1,))
            production = np.repeat(production, self.nd)
            production.shape = (self.shape_nxyz, self.nd)

        return delayed_production / production

class OuterState(State):

    def __init__(self, states):
        super().__init__(states)

        # Initialize Solver class attributes
        self._mgxs_lib = None
        self._chi_delayed_by_delayed_group = False
        self._chi_delayed_by_mesh = False
        self._method = 'ADIABATIC'
        self._mgxs_loaded = False
        self._shape = None

    @property
    def shape(self):
        self._shape.shape = (self.shape_nxyz, self.ng)
        return self._shape

    @property
    def mgxs_lib(self):
        return self._mgxs_lib

    @property
    def chi_delayed_by_delayed_group(self):
        return self._chi_delayed_by_delayed_group

    @property
    def chi_delayed_by_mesh(self):
        return self._chi_delayed_by_mesh

    @property
    def use_agd(self):
        return self._use_agd

    @property
    def use_pcmfd(self):
        return self._use_pcmfd

    @property
    def mgxs_loaded(self):
        return self._mgxs_loaded

    @shape.setter
    def shape(self, shape):
        self._shape = copy.deepcopy(shape)
        self._shape.shape = (self.shape_nxyz, self.ng)

    @mgxs_lib.setter
    def mgxs_lib(self, mgxs_lib):
        self._mgxs_lib = mgxs_lib

    @chi_delayed_by_delayed_group.setter
    def chi_delayed_by_delayed_group(self, chi_delayed_by_delayed_group):
        self._chi_delayed_by_delayed_group = chi_delayed_by_delayed_group

    @chi_delayed_by_mesh.setter
    def chi_delayed_by_mesh(self, chi_delayed_by_mesh):
        self._chi_delayed_by_mesh = chi_delayed_by_mesh

    @use_agd.setter
    def use_agd(self, use_agd):
        self._use_agd = use_agd

    @use_pcmfd.setter
    def use_pcmfd(self, use_pcmfd):
        self._use_pcmfd = use_pcmfd

    @mgxs_loaded.setter
    def mgxs_loaded(self, mgxs_loaded):
        self._mgxs_loaded = mgxs_loaded

    @property
    def inscatter(self):
        if not self.mgxs_loaded:
            self._inscatter = self.mgxs_lib['nu-scatter matrix'].get_condensed_xs(self.energy_groups).get_xs(row_column='outin')

        self._inscatter.shape = (self.shape_nxyz, self.ng, self.ng)
        return self._inscatter

    @property
    def outscatter(self):
        return self.inscatter.sum(axis=1)

    @property
    def absorption(self):
        if not self.mgxs_loaded:
            self._absorption = self.mgxs_lib['absorption'].get_condensed_xs(self.energy_groups).get_xs()

        self._absorption.shape = (self.shape_nxyz, self.ng)
        return self._absorption

    @property
    def kappa_fission(self):
        if not self.mgxs_loaded:
            self._kappa_fission = self.mgxs_lib['kappa-fission'].get_condensed_xs(self.energy_groups).get_xs()

        self._kappa_fission.shape = (self.shape_nxyz, self.ng)
        return self._kappa_fission

    @property
    def chi_prompt(self):
        if not self.mgxs_loaded:
            self._chi_prompt = self.mgxs_lib['chi-prompt'].get_condensed_xs(self.energy_groups).get_xs()

        self._chi_prompt.shape = (self.shape_nxyz, self.ng)
        return self._chi_prompt

    @property
    def prompt_nu_fission(self):
        if not self.mgxs_loaded:
            self._prompt_nu_fission = self.mgxs_lib['prompt-nu-fission'].get_condensed_xs(self.energy_groups).get_xs()

        self._prompt_nu_fission.shape = (self.shape_nxyz, self.ng)
        return self._prompt_nu_fission

    @property
    def nu_fission(self):
        return self.delayed_nu_fission.sum(axis=1) + self.prompt_nu_fission

    @property
    def chi_delayed(self):

        if not self.mgxs_loaded:
            self._chi_delayed = self.mgxs_lib['chi-delayed'].get_condensed_xs(self.energy_groups).get_xs()

            if self.chi_delayed_by_mesh:
                if not self.chi_delayed_by_delayed_group:
                    self._chi_delayed.shape = (self.shape_nxyz, self.ng)
                    self._chi_delayed = np.tile(self._chi_delayed, self.nd)
            else:
                if self.chi_delayed_by_delayed_group:
                    self._chi_delayed = np.tile(self._chi_delayed.flatten(), self.shape_nxyz)
                else:
                    self._chi_delayed = np.tile(self._chi_delayed.flatten(), self.shape_nxyz)
                    self._chi_delayed.shape = (self.shape_nxyz, self.ng)
                    self._chi_delayed = np.tile(self._chi_delayed, self.nd)

        self._chi_delayed.shape = (self.shape_nxyz, self.nd, self.ng)
        return self._chi_delayed

    @property
    def delayed_nu_fission(self):
        if not self.mgxs_loaded:
            self._delayed_nu_fission = self.mgxs_lib['delayed-nu-fission'].get_condensed_xs(self.energy_groups).get_xs()

        self._delayed_nu_fission.shape = (self.shape_nxyz, self.nd, self.ng)
        return self._delayed_nu_fission

    @property
    def inverse_velocity(self):
        if not self.mgxs_loaded:
            self._inverse_velocity = self.mgxs_lib['inverse-velocity'].get_condensed_xs(self.energy_groups).get_xs()

        self._inverse_velocity.shape = (self.shape_nxyz, self.ng)
        return self._inverse_velocity

    @property
    def decay_rate(self):
        if not self.mgxs_loaded:
            self._decay_rate = self.mgxs_lib['decay-rate'].get_xs()
            self._decay_rate[self._decay_rate < 1.e-5] = 0.

        self._decay_rate.shape = (self.shape_nxyz, self.nd)
        return self._decay_rate

    @property
    def flux_tallied(self):
        if not self.mgxs_loaded:
            self._flux_tallied = self.mgxs_lib['absorption'].get_condensed_xs(self.energy_groups).tallies['flux'].get_values()
            self._flux_tallied.shape = (self.shape_nxyz, self.ng)
            self._flux_tallied = self._flux_tallied[:, ::-1]

        self._flux_tallied.shape = (self.shape_nxyz, self.ng)
        return self._flux_tallied

    @property
    def amplitude_flux_tallied(self):
        if not self.mgxs_loaded:
            amp_flux = self.mgxs_lib['amplitude-kappa-fission'].get_condensed_xs(self.energy_groups).tallies['flux'].get_values()
            amp_flux.shape = self.tally_zyxg
            amp_flux = amp_flux[..., ::-1]
            self._amplitude_flux_tallied = openmc.kinetics.map_array(amp_flux, self.tally_zyxg, self.amplitude_zyxg, False)

        self._amplitude_flux_tallied.shape = (self.amplitude_nxyz, self.ng)
        return self._amplitude_flux_tallied

    @property
    def current_tallied(self):
        if not self.mgxs_loaded:
            current = self.mgxs_lib['current'].get_condensed_xs(self.energy_groups).get_xs() 

            # Condense down the the coarse mesh
            tally_shape = self.tally_zyxg + (12,)
            
            amp_shape = self.amplitude_zyxg + (12,)

            self._current_tallied = openmc.kinetics.surface_integral(current, tally_shape, amp_shape)

        self._current_tallied.shape = (self.amplitude_nxyz, self.ng, 12)

        return self._current_tallied

    @property
    def diffusion_coefficient(self):
        if not self.mgxs_loaded:
            dif_coef = self.mgxs_lib['diffusion-coefficient']\
                .get_condensed_xs(self.energy_groups)
            flux = dif_coef.tallies['flux (tracklength)'].get_values()
            dif_coef = dif_coef.get_xs()
            dif_coef.shape = self.tally_zyxg
            flux.shape = self.tally_zyxg
            flux = flux[..., ::-1]

            # Spatially condense the transport cross section
            sig_tr = 1. / (3. * dif_coef)
            sig_tr = openmc.kinetics.nan_inf_to_zero(sig_tr)
            sig_tr = openmc.kinetics.map_array(sig_tr * flux, self.tally_zyxg, self.amplitude_zyxg, False)
            flux   = openmc.kinetics.map_array(         flux, self.tally_zyxg, self.amplitude_zyxg, False)
            sig_tr /= flux
            self._diffusion_coefficient = openmc.kinetics.nan_inf_to_zero(1.0 / (3.0 * sig_tr))

        self._diffusion_coefficient.shape = (self.amplitude_nxyz, self.ng)
        return self._diffusion_coefficient

    def extract_shape(self):

        # Get the current power
        power = self.core_power_density

        # Get the tallied pin-wise flux
        flux = self.flux_tallied
        flux.shape = self.shape_zyxg

        # Expand the coarse mesh amplitude to the pin mesh
        amp = openmc.kinetics.map_array(self.amplitude, self.amplitude_zyxg, self.shape_zyxg, True)

        # Compute the pin shape
        self.shape = flux / amp
        self.shape *= power / self.core_power_density

    @property
    def dump_to_log_file(self):

        time_point = str(self.clock.times[self.time_point])
        f = h5py.File(self._log_file, 'a')
        if time_point not in f['OUTER_STEPS'].keys():
            f['OUTER_STEPS'].require_group(time_point)

        if 'shape' not in f['OUTER_STEPS'][time_point].keys():
            f['OUTER_STEPS'][time_point].create_dataset('shape', data=self.shape)
            f['OUTER_STEPS'][time_point].create_dataset('power', data=self.power)
            f['OUTER_STEPS'][time_point].create_dataset('kappa_fission', data=self.kappa_fission)
        else:
            shape = f['OUTER_STEPS'][time_point]['shape']
            shape[...] = self.shape
            power = f['OUTER_STEPS'][time_point]['power']
            power[...] = self.power
            kappa_fission = f['OUTER_STEPS'][time_point]['kappa_fission']
            kappa_fission[...] = self.kappa_fission

        f.close()

    def compute_initial_precursor_concentration(self):
        flux = np.tile(self.flux, self.nd).flatten()
        del_fis_rate = self.delayed_nu_fission.flatten() * flux
        del_fis_rate.shape = (self.shape_nxyz, self.nd, self.ng)
        precursors = del_fis_rate.sum(axis=2) / self.decay_rate / self.k_crit
        self.precursors = openmc.kinetics.nan_inf_to_zero(precursors)

    def load_mgxs(self):
        self.mgxs_loaded = False
        self.delayed_nu_fission
        self.inscatter
        self.absorption
        self.chi_prompt
        self.prompt_nu_fission
        self.chi_delayed
        self.kappa_fission
        self.inverse_velocity
        self.decay_rate
        self.flux_tallied
        self.amplitude_flux_tallied
        self.current_tallied
        self.diffusion_coefficient
        self.mgxs_loaded = True

    def initialize_mgxs(self):
        """Initialize all the tallies for the problem.

        """

        # Instantiate a list of the delayed groups
        delayed_groups = list(range(1,self.nd + 1))

        # Create elements and ordered dicts and initialize to None
        self._mgxs_lib = OrderedDict()

        mgxs_types = ['absorption', 'diffusion-coefficient', 'decay-rate',
                      'kappa-fission', 'chi-prompt', 'chi-delayed', 'inverse-velocity',
                      'prompt-nu-fission', 'current', 'delayed-nu-fission']


        self._mgxs_lib['amplitude-kappa-fission'] = openmc.mgxs.MGXS.get_mgxs(
            'kappa-fission', domain=self.tally_mesh, domain_type='mesh',
            energy_groups=self.tally_groups, by_nuclide=False,
            name= self.time_point + ' - amplitude-kappa-fission')

        mgxs_types.append('nu-scatter matrix')

        # Populate the MGXS in the MGXS lib
        for mgxs_type in mgxs_types:
            mesh = self.shape_mesh
            if mgxs_type == 'diffusion-coefficient':
                self._mgxs_lib[mgxs_type] = openmc.mgxs.MGXS.get_mgxs(
                    mgxs_type, domain=self.tally_mesh, domain_type='mesh',
                    energy_groups=self.fine_groups, by_nuclide=False,
                    name= self.time_point + ' - ' + mgxs_type)
            elif mgxs_type == 'current':
                self._mgxs_lib[mgxs_type] = openmc.mgxs.MGXS.get_mgxs(
                    mgxs_type, domain=self.tally_mesh, domain_type='mesh',
                    energy_groups=self.tally_groups, by_nuclide=False,
                    name= self.time_point + ' - ' + mgxs_type)
            elif 'nu-scatter matrix' in mgxs_type:
                self._mgxs_lib[mgxs_type] = openmc.mgxs.MGXS.get_mgxs(
                    mgxs_type, domain=mesh, domain_type='mesh',
                    energy_groups=self.tally_groups, by_nuclide=False,
                    name= self.time_point + ' - ' + mgxs_type)
                self._mgxs_lib[mgxs_type].correction = None
            elif mgxs_type == 'decay-rate':
                self._mgxs_lib[mgxs_type] = openmc.mgxs.MDGXS.get_mgxs(
                    mgxs_type, domain=mesh, domain_type='mesh',
                    energy_groups=self.one_group,
                    delayed_groups=delayed_groups, by_nuclide=False,
                    name= self.time_point + ' - ' + mgxs_type)
            elif mgxs_type == 'chi-prompt':
                self._mgxs_lib[mgxs_type] = openmc.mgxs.MGXS.get_mgxs(
                    mgxs_type, domain=mesh, domain_type='mesh',
                    energy_groups=self.tally_groups, by_nuclide=False,
                    name= self.time_point + ' - ' + mgxs_type)
            elif mgxs_type in openmc.mgxs.MGXS_TYPES:
                self._mgxs_lib[mgxs_type] = openmc.mgxs.MGXS.get_mgxs(
                    mgxs_type, domain=mesh, domain_type='mesh',
                    energy_groups=self.tally_groups, by_nuclide=False,
                    name= self.time_point + ' - ' + mgxs_type)
            elif mgxs_type == 'chi-delayed':
                if self.chi_delayed_by_delayed_group:
                    if self.chi_delayed_by_mesh:
                        self._mgxs_lib[mgxs_type] = openmc.mgxs.MDGXS.get_mgxs(
                            mgxs_type, domain=mesh, domain_type='mesh',
                            energy_groups=self.tally_groups,
                            delayed_groups=delayed_groups, by_nuclide=False,
                            name= self.time_point + ' - ' + mgxs_type)
                    else:
                        self._mgxs_lib[mgxs_type] = openmc.mgxs.MDGXS.get_mgxs(
                            mgxs_type, domain=self.unity_mesh, domain_type='mesh',
                            energy_groups=self.tally_groups,
                            delayed_groups=delayed_groups, by_nuclide=False,
                            name= self.time_point + ' - ' + mgxs_type)
                else:
                    if self.chi_delayed_by_mesh:
                        self._mgxs_lib[mgxs_type] = openmc.mgxs.MDGXS.get_mgxs(
                            mgxs_type, domain=mesh, domain_type='mesh',
                            energy_groups=self.tally_groups, by_nuclide=False,
                            name= self.time_point + ' - ' + mgxs_type)
                    else:
                        self._mgxs_lib[mgxs_type] = openmc.mgxs.MDGXS.get_mgxs(
                            mgxs_type, domain=self.unity_mesh, domain_type='mesh',
                            energy_groups=self.tally_groups, by_nuclide=False,
                            name= self.time_point + ' - ' + mgxs_type)
            elif mgxs_type in openmc.mgxs.MDGXS_TYPES:
                self._mgxs_lib[mgxs_type] = openmc.mgxs.MDGXS.get_mgxs(
                    mgxs_type, domain=mesh, domain_type='mesh',
                    energy_groups=self.tally_groups,
                    delayed_groups=delayed_groups, by_nuclide=False,
                    name= self.time_point + ' - ' + mgxs_type)

        self.mgxs_loaded = False

    @property
    def coupling_terms(self):

        # Get the dimensions of the mesh
        nz , ny , nx  = self.amplitude_dimension
        if self.amplitude_mesh.width:
            dx , dy , dz  = self.amplitude_mesh.width
        else:
            dx , dy , dz  = [(i-j)/k for i,j,k in zip(
                self.amplitude_mesh.upper_right, 
                self.amplitude_mesh.lower_left, 
                self.amplitude_mesh.dimension)]
        ng            = self.ng

        dx /= nx
        dy /= ny
        dz /= nz

        # Get the array of the surface-integrated surface net currents
        partial_current = copy.deepcopy(self.current_tallied)
        partial_current = partial_current.reshape(int(np.prod(partial_current.shape) / 12), 12)
        net_current = partial_current[..., 0:12:2] - partial_current[..., 1:12:2]
        net_current[:, 0:6:2] = -net_current[:, 0:6:2]
        net_current.shape = self.amplitude_zyxg + (6,)
        partial_current.shape = self.amplitude_zyxg + (12,)

        # Get the flux
        flux = copy.deepcopy(self.amplitude_flux_tallied)
        flux.shape = self.amplitude_zyxg

        # Create an array of the neighbor cell fluxes
        flux_nbr = np.zeros((nz, ny, nx, ng, 6))
        flux_nbr[:  , :  , 1: , :, 0] = flux[:  , :  , :-1, :]
        flux_nbr[:  , :  , :-1, :, 1] = flux[:  , :  , 1: , :]
        flux_nbr[:  , 1: , :  , :, 2] = flux[:  , :-1, :  , :]
        flux_nbr[:  , :-1, :  , :, 3] = flux[:  , 1: , :  , :]
        flux_nbr[1: , :  , :  , :, 4] = flux[:-1, :  , :  , :]
        flux_nbr[:-1, :  , :  , :, 5] = flux[1: , :  , :  , :]

        # Get the diffusion coefficients tally
        dc       = copy.deepcopy(self.diffusion_coefficient)
        dc.shape = self.amplitude_zyxg
        dc_nbr   = np.zeros(self.amplitude_zyxg + (6,))

        # Create array of neighbor cell diffusion coefficients
        dc_nbr[:  , :  , 1: , :, 0] = dc[:  , :  , :-1, :]
        dc_nbr[:  , :  , :-1, :, 1] = dc[:  , :  , 1: , :]
        dc_nbr[:  , 1: , :  , :, 2] = dc[:  , :-1, :  , :]
        dc_nbr[:  , :-1, :  , :, 3] = dc[:  , 1: , :  , :]
        dc_nbr[1: , :  , :  , :, 4] = dc[:-1, :  , :  , :]
        dc_nbr[:-1, :  , :  , :, 5] = dc[1: , :  , :  , :]

        # Compute the linear finite difference diffusion term for interior surfaces
        dc_linear = np.zeros((nz, ny, nx, ng, 6))
        dc_linear[:  , :  , 1: , :, 0] = 2 * dc_nbr[:  , :  , 1: , :, 0] * dc[:  , :  , 1: , :] / (dc_nbr[:  , :  , 1: , :, 0] * dx + dc[:  , :  , 1: , :] * dx)
        dc_linear[:  , :  , :-1, :, 1] = 2 * dc_nbr[:  , :  , :-1, :, 1] * dc[:  , :  , :-1, :] / (dc_nbr[:  , :  , :-1, :, 1] * dx + dc[:  , :  , :-1, :] * dx)
        dc_linear[:  , 1: , :  , :, 2] = 2 * dc_nbr[:  , 1: , :  , :, 2] * dc[:  , 1: , :  , :] / (dc_nbr[:  , 1: , :  , :, 2] * dy + dc[:  , 1: , :  , :] * dy)
        dc_linear[:  , :-1, :  , :, 3] = 2 * dc_nbr[:  , :-1, :  , :, 3] * dc[:  , :-1, :  , :] / (dc_nbr[:  , :-1, :  , :, 3] * dy + dc[:  , :-1, :  , :] * dy)
        dc_linear[1: , :  , :  , :, 4] = 2 * dc_nbr[1: , :  , :  , :, 4] * dc[1: , :  , :  , :] / (dc_nbr[1: , :  , :  , :, 4] * dz + dc[1: , :  , :  , :] * dz)
        dc_linear[:-1, :  , :  , :, 5] = 2 * dc_nbr[:-1, :  , :  , :, 5] * dc[:-1, :  , :  , :] / (dc_nbr[:-1, :  , :  , :, 5] * dz + dc[:-1, :  , :  , :] * dz)

        if self.use_agd:
            dc_linear[:  , :  , 1: , :, 0] += 0.25 * dx
            dc_linear[:  , :  , :-1, :, 1] += 0.25 * dx
            dc_linear[:  , 1: , :  , :, 2] += 0.25 * dy
            dc_linear[:  , :-1, :  , :, 3] += 0.25 * dy
            dc_linear[1: , :  , :  , :, 4] += 0.25 * dz
            dc_linear[:-1, :  , :  , :, 5] += 0.25 * dz

        # Make any cells that have no dif coef or flux tally highly diffusive
        dc_linear[np.isnan(dc_linear)] = 0.0
        dc_linear[dc_linear == 0.] = 0.0

        # Compute the non-linear finite difference diffusion term for interior surfaces
        if self.use_pcmfd:
            dc_nonlinear = np.zeros((nz, ny, nx, ng, 12))
            dc_nonlinear[..., 0 ] = (-dc_linear[..., 0] * (-flux_nbr[..., 0] + flux) + 2 * partial_current[..., 0 ]) / (2 * flux)
            dc_nonlinear[..., 1 ] = ( dc_linear[..., 0] * (-flux_nbr[..., 0] + flux) + 2 * partial_current[..., 1 ]) / (2 * flux_nbr[..., 0])
            dc_nonlinear[..., 2 ] = ( dc_linear[..., 1] * ( flux_nbr[..., 1] - flux) + 2 * partial_current[..., 2 ]) / (2 * flux)
            dc_nonlinear[..., 3 ] = (-dc_linear[..., 1] * ( flux_nbr[..., 1] - flux) + 2 * partial_current[..., 3 ]) / (2 * flux_nbr[..., 1])
            dc_nonlinear[..., 4 ] = (-dc_linear[..., 2] * (-flux_nbr[..., 2] + flux) + 2 * partial_current[..., 4 ]) / (2 * flux)
            dc_nonlinear[..., 5 ] = ( dc_linear[..., 2] * (-flux_nbr[..., 2] + flux) + 2 * partial_current[..., 5 ]) / (2 * flux_nbr[..., 2])
            dc_nonlinear[..., 6 ] = ( dc_linear[..., 3] * ( flux_nbr[..., 3] - flux) + 2 * partial_current[..., 6 ]) / (2 * flux)
            dc_nonlinear[..., 7 ] = (-dc_linear[..., 3] * ( flux_nbr[..., 3] - flux) + 2 * partial_current[..., 7 ]) / (2 * flux_nbr[..., 3])
            dc_nonlinear[..., 8 ] = (-dc_linear[..., 4] * (-flux_nbr[..., 4] + flux) + 2 * partial_current[..., 8 ]) / (2 * flux)
            dc_nonlinear[..., 9 ] = ( dc_linear[..., 4] * (-flux_nbr[..., 4] + flux) + 2 * partial_current[..., 9 ]) / (2 * flux_nbr[..., 4])
            dc_nonlinear[..., 10] = ( dc_linear[..., 5] * ( flux_nbr[..., 5] - flux) + 2 * partial_current[..., 10]) / (2 * flux)
            dc_nonlinear[..., 11] = (-dc_linear[..., 5] * ( flux_nbr[..., 5] - flux) + 2 * partial_current[..., 11]) / (2 * flux_nbr[..., 5])
        else:
            dc_nonlinear = np.zeros((nz, ny, nx, ng, 6))
            dc_nonlinear[..., 0] = (-dc_linear[..., 0] * (-flux_nbr[..., 0] + flux) - net_current[..., 0]) / (flux_nbr[..., 0] + flux)
            dc_nonlinear[..., 1] = (-dc_linear[..., 1] * ( flux_nbr[..., 1] - flux) - net_current[..., 1]) / (flux_nbr[..., 1] + flux)
            dc_nonlinear[..., 2] = (-dc_linear[..., 2] * (-flux_nbr[..., 2] + flux) - net_current[..., 2]) / (flux_nbr[..., 2] + flux)
            dc_nonlinear[..., 3] = (-dc_linear[..., 3] * ( flux_nbr[..., 3] - flux) - net_current[..., 3]) / (flux_nbr[..., 3] + flux)
            dc_nonlinear[..., 4] = (-dc_linear[..., 4] * (-flux_nbr[..., 4] + flux) - net_current[..., 4]) / (flux_nbr[..., 4] + flux)
            dc_nonlinear[..., 5] = (-dc_linear[..., 5] * ( flux_nbr[..., 5] - flux) - net_current[..., 5]) / (flux_nbr[..., 5] + flux)

        # Ensure there are no nans
        dc_nonlinear[np.isnan(dc_nonlinear)] = 0.

        # Check for diagonal dominance
        if False:

            # Make a mask of the location of all terms that need to be corrected (dd_mask) and
            # terms that don't need to be corrected (nd_mask)
            dd_mask = (np.abs(dc_nonlinear) > dc_linear)
            nd_mask = (dd_mask == False)

            # Save arrays as to whether the correction term is positive or negative.
            sign = np.abs(dc_nonlinear) / dc_nonlinear
            sign_pos  = (dc_nonlinear > 0.)
            if self.use_pcmfd:
                sense_pos = np.zeros((nz, ny, nx, ng, 12))
                sense_pos[..., 0:12:4] = False
                sense_pos[..., 1:12:4] = False
                sense_pos[..., 2:12:4] = True
                sense_pos[..., 3:12:4] = True
            else:
                sense_pos = np.zeros((nz, ny, nx, ng, 6))
                sense_pos[..., 0:6:2] = False
                sense_pos[..., 1:6:2] = True

            sign_sense = (sign_pos == sense_pos)
            not_sign_sense = (sign_sense == False)

            # Correct dc_linear
            if self.use_pcmfd:
                dc_linear[:  , :  , 1: , :, 0] = nd_mask[:  , :  , 1: , :, 0] * dc_linear[:  , :  , 1: , :, 0] + dd_mask[:  , :  , 1: , :, 0] * (sign_sense[:  , :  , 1: , :, 0] * np.abs(net_current[:  , :  , 1: , :, 0] / (2 * flux_nbr[:  , :  , 1: , :, 0])) + not_sign_sense[:  , :  , 1: , :, 0] * np.abs(net_current[:  , :  , 1: , :, 0] / (2 * flux[:  , :  , 1: , :])))
                dc_linear[:  , :  , :-1, :, 1] = nd_mask[:  , :  , :-1, :, 1] * dc_linear[:  , :  , :-1, :, 1] + dd_mask[:  , :  , :-1, :, 1] * (sign_sense[:  , :  , :-1, :, 1] * np.abs(net_current[:  , :  , :-1, :, 1] / (2 * flux_nbr[:  , :  , :-1, :, 1])) + not_sign_sense[:  , :  , :-1, :, 1] * np.abs(net_current[:  , :  , :-1, :, 1] / (2 * flux[:  , :  , :-1, :])))
                dc_linear[:  , 1: , :  , :, 2] = nd_mask[:  , 1: , :  , :, 2] * dc_linear[:  , 1: , :  , :, 2] + dd_mask[:  , 1: , :  , :, 2] * (sign_sense[:  , 1: , :  , :, 2] * np.abs(net_current[:  , 1: , :  , :, 2] / (2 * flux_nbr[:  , 1: , :  , :, 2])) + not_sign_sense[:  , 1: , :  , :, 2] * np.abs(net_current[:  , 1: , :  , :, 2] / (2 * flux[:  , 1: , :  , :])))
                dc_linear[:  , :-1, :  , :, 3] = nd_mask[:  , :-1, :  , :, 3] * dc_linear[:  , :-1, :  , :, 3] + dd_mask[:  , :-1, :  , :, 3] * (sign_sense[:  , :-1, :  , :, 3] * np.abs(net_current[:  , :-1, :  , :, 3] / (2 * flux_nbr[:  , :-1, :  , :, 3])) + not_sign_sense[:  , :-1, :  , :, 3] * np.abs(net_current[:  , :-1, :  , :, 3] / (2 * flux[:  , :-1, :  , :])))
                dc_linear[1: , :  , :  , :, 4] = nd_mask[1: , :  , :  , :, 4] * dc_linear[1: , :  , :  , :, 4] + dd_mask[1: , :  , :  , :, 4] * (sign_sense[1: , :  , :  , :, 4] * np.abs(net_current[1: , :  , :  , :, 4] / (2 * flux_nbr[1: , :  , :  , :, 4])) + not_sign_sense[1: , :  , :  , :, 4] * np.abs(net_current[1: , :  , :  , :, 4] / (2 * flux[1: , :  , :  , :])))
                dc_linear[:-1, :  , :  , :, 5] = nd_mask[:-1, :  , :  , :, 5] * dc_linear[:-1, :  , :  , :, 5] + dd_mask[:-1, :  , :  , :, 5] * (sign_sense[:-1, :  , :  , :, 5] * np.abs(net_current[:-1, :  , :  , :, 5] / (2 * flux_nbr[:-1, :  , :  , :, 5])) + not_sign_sense[:-1, :  , :  , :, 5] * np.abs(net_current[:-1, :  , :  , :, 5] / (2 * flux[:-1, :  , :  , :])))

                dc_nonlinear[:  , :  , 1: , :, 0] = nd_mask[:  , :  , 1: , :, 0] * dc_nonlinear[:  , :  , 1: , :, 0] + dd_mask[:  , :  , 1: , :, 0] * sign[:  , :  , 1: , :, 0] * dc_linear[:  , :  , 1: , :, 0]
                dc_nonlinear[:  , :  , :-1, :, 1] = nd_mask[:  , :  , :-1, :, 1] * dc_nonlinear[:  , :  , :-1, :, 1] + dd_mask[:  , :  , :-1, :, 1] * sign[:  , :  , :-1, :, 1] * dc_linear[:  , :  , :-1, :, 1]
                dc_nonlinear[:  , 1: , :  , :, 2] = nd_mask[:  , 1: , :  , :, 2] * dc_nonlinear[:  , 1: , :  , :, 2] + dd_mask[:  , 1: , :  , :, 2] * sign[:  , 1: , :  , :, 2] * dc_linear[:  , 1: , :  , :, 2]
                dc_nonlinear[:  , :-1, :  , :, 3] = nd_mask[:  , :-1, :  , :, 3] * dc_nonlinear[:  , :-1, :  , :, 3] + dd_mask[:  , :-1, :  , :, 3] * sign[:  , :-1, :  , :, 3] * dc_linear[:  , :-1, :  , :, 3]
                dc_nonlinear[1: , :  , :  , :, 4] = nd_mask[1: , :  , :  , :, 4] * dc_nonlinear[1: , :  , :  , :, 4] + dd_mask[1: , :  , :  , :, 4] * sign[1: , :  , :  , :, 4] * dc_linear[1: , :  , :  , :, 4]
                dc_nonlinear[:-1, :  , :  , :, 5] = nd_mask[:-1, :  , :  , :, 5] * dc_nonlinear[:-1, :  , :  , :, 5] + dd_mask[:-1, :  , :  , :, 5] * sign[:-1, :  , :  , :, 5] * dc_linear[:-1, :  , :  , :, 5]

            else:
                dc_linear[:  , :  , 1: , :, 0] = nd_mask[:  , :  , 1: , :, 0] * dc_linear[:  , :  , 1: , :, 0] + dd_mask[:  , :  , 1: , :, 0] * (sign_sense[:  , :  , 1: , :, 0] * np.abs(net_current[:  , :  , 1: , :, 0] / (2 * flux_nbr[:  , :  , 1: , :, 0])) + not_sign_sense[:  , :  , 1: , :, 0] * np.abs(net_current[:  , :  , 1: , :, 0] / (2 * flux[:  , :  , 1: , :])))
                dc_linear[:  , :  , :-1, :, 1] = nd_mask[:  , :  , :-1, :, 1] * dc_linear[:  , :  , :-1, :, 1] + dd_mask[:  , :  , :-1, :, 1] * (sign_sense[:  , :  , :-1, :, 1] * np.abs(net_current[:  , :  , :-1, :, 1] / (2 * flux_nbr[:  , :  , :-1, :, 1])) + not_sign_sense[:  , :  , :-1, :, 1] * np.abs(net_current[:  , :  , :-1, :, 1] / (2 * flux[:  , :  , :-1, :])))
                dc_linear[:  , 1: , :  , :, 2] = nd_mask[:  , 1: , :  , :, 2] * dc_linear[:  , 1: , :  , :, 2] + dd_mask[:  , 1: , :  , :, 2] * (sign_sense[:  , 1: , :  , :, 2] * np.abs(net_current[:  , 1: , :  , :, 2] / (2 * flux_nbr[:  , 1: , :  , :, 2])) + not_sign_sense[:  , 1: , :  , :, 2] * np.abs(net_current[:  , 1: , :  , :, 2] / (2 * flux[:  , 1: , :  , :])))
                dc_linear[:  , :-1, :  , :, 3] = nd_mask[:  , :-1, :  , :, 3] * dc_linear[:  , :-1, :  , :, 3] + dd_mask[:  , :-1, :  , :, 3] * (sign_sense[:  , :-1, :  , :, 3] * np.abs(net_current[:  , :-1, :  , :, 3] / (2 * flux_nbr[:  , :-1, :  , :, 3])) + not_sign_sense[:  , :-1, :  , :, 3] * np.abs(net_current[:  , :-1, :  , :, 3] / (2 * flux[:  , :-1, :  , :])))
                dc_linear[1: , :  , :  , :, 4] = nd_mask[1: , :  , :  , :, 4] * dc_linear[1: , :  , :  , :, 4] + dd_mask[1: , :  , :  , :, 4] * (sign_sense[1: , :  , :  , :, 4] * np.abs(net_current[1: , :  , :  , :, 4] / (2 * flux_nbr[1: , :  , :  , :, 4])) + not_sign_sense[1: , :  , :  , :, 4] * np.abs(net_current[1: , :  , :  , :, 4] / (2 * flux[1: , :  , :  , :])))
                dc_linear[:-1, :  , :  , :, 5] = nd_mask[:-1, :  , :  , :, 5] * dc_linear[:-1, :  , :  , :, 5] + dd_mask[:-1, :  , :  , :, 5] * (sign_sense[:-1, :  , :  , :, 5] * np.abs(net_current[:-1, :  , :  , :, 5] / (2 * flux_nbr[:-1, :  , :  , :, 5])) + not_sign_sense[:-1, :  , :  , :, 5] * np.abs(net_current[:-1, :  , :  , :, 5] / (2 * flux[:-1, :  , :  , :])))

                dc_nonlinear[:  , :  , 1: , :, 0] = nd_mask[:  , :  , 1: , :, 0] * dc_nonlinear[:  , :  , 1: , :, 0] + dd_mask[:  , :  , 1: , :, 0] * sign[:  , :  , 1: , :, 0] * dc_linear[:  , :  , 1: , :, 0]
                dc_nonlinear[:  , :  , :-1, :, 1] = nd_mask[:  , :  , :-1, :, 1] * dc_nonlinear[:  , :  , :-1, :, 1] + dd_mask[:  , :  , :-1, :, 1] * sign[:  , :  , :-1, :, 1] * dc_linear[:  , :  , :-1, :, 1]
                dc_nonlinear[:  , 1: , :  , :, 2] = nd_mask[:  , 1: , :  , :, 2] * dc_nonlinear[:  , 1: , :  , :, 2] + dd_mask[:  , 1: , :  , :, 2] * sign[:  , 1: , :  , :, 2] * dc_linear[:  , 1: , :  , :, 2]
                dc_nonlinear[:  , :-1, :  , :, 3] = nd_mask[:  , :-1, :  , :, 3] * dc_nonlinear[:  , :-1, :  , :, 3] + dd_mask[:  , :-1, :  , :, 3] * sign[:  , :-1, :  , :, 3] * dc_linear[:  , :-1, :  , :, 3]
                dc_nonlinear[1: , :  , :  , :, 4] = nd_mask[1: , :  , :  , :, 4] * dc_nonlinear[1: , :  , :  , :, 4] + dd_mask[1: , :  , :  , :, 4] * sign[1: , :  , :  , :, 4] * dc_linear[1: , :  , :  , :, 4]
                dc_nonlinear[:-1, :  , :  , :, 5] = nd_mask[:-1, :  , :  , :, 5] * dc_nonlinear[:-1, :  , :  , :, 5] + dd_mask[:-1, :  , :  , :, 5] * sign[:-1, :  , :  , :, 5] * dc_linear[:-1, :  , :  , :, 5]


        # Reshape the diffusion coefficient array
        dc_linear.shape    = (nx*ny*nz*ng, 6)
        if self.use_pcmfd:
            dc_nonlinear.shape = (nx*ny*nz*ng, 12)
        else:
            dc_nonlinear.shape = (nx*ny*nz*ng, 6)

        # reshape the flux shape
        shape       = openmc.kinetics.map_array(self.shape, self.shape_zyxg, self.amplitude_zyxg, False)
        shape.shape = (nz, ny, nx, ng)
        shape_nbr = np.zeros((nz, ny, nx, ng, 6))
        shape_nbr[:  , :  , 1: , :, 0] = shape[:  , :  , :-1, :]
        shape_nbr[:  , :  , :-1, :, 1] = shape[:  , :  , 1: , :]
        shape_nbr[:  , 1: , :  , :, 2] = shape[:  , :-1, :  , :]
        shape_nbr[:  , :-1, :  , :, 3] = shape[:  , 1: , :  , :]
        shape_nbr[1: , :  , :  , :, 4] = shape[:-1, :  , :  , :]
        shape_nbr[:-1, :  , :  , :, 5] = shape[1: , :  , :  , :]
        shape.shape     = (nx*ny*nz*ng,)
        shape_nbr.shape = (nx*ny*nz*ng, 6)

        # Set the diagonal
        dc_linear_diag    =  dc_linear   [:, 1:6:2].sum(axis=1) + dc_linear   [:, 0:6:2].sum(axis=1)

        if self.use_pcmfd:
            dc_nonlinear_diag =  dc_nonlinear[:, 2:12:4].sum(axis=1) + dc_nonlinear[:, 0:12:4].sum(axis=1)
        else:
            dc_nonlinear_diag = -dc_nonlinear[:, 1:6:2].sum(axis=1) + dc_nonlinear[:, 0:6:2].sum(axis=1)

        dc_linear_data    = [dc_linear_diag * shape]
        dc_nonlinear_data = [dc_nonlinear_diag * shape]
        diags             = [0]

        # Zero boundary dc_nonlinear
        if self.use_pcmfd:
            dc_nonlinear.shape = (nz, ny, nx, ng, 12)
            dc_nonlinear[:  ,  :,  0, :, 0 ] = 0.
            dc_nonlinear[:  ,  :,  0, :, 1 ] = 0.
            dc_nonlinear[:  ,  :, -1, :, 2 ] = 0.
            dc_nonlinear[:  ,  :, -1, :, 3 ] = 0.
            dc_nonlinear[:  ,  0,  :, :, 4 ] = 0.
            dc_nonlinear[:  ,  0,  :, :, 5 ] = 0.
            dc_nonlinear[:  , -1,  :, :, 6 ] = 0.
            dc_nonlinear[:  , -1,  :, :, 7 ] = 0.
            dc_nonlinear[0  ,  :,  :, :, 8 ] = 0.
            dc_nonlinear[0  ,  :,  :, :, 9 ] = 0.
            dc_nonlinear[-1 ,  :,  :, :, 10] = 0.
            dc_nonlinear[-1 ,  :,  :, :, 11] = 0.
            dc_nonlinear.shape = (nz*ny*nx*ng, 12)
        else:
            dc_nonlinear.shape = (nz, ny, nx, ng, 6)
            dc_nonlinear[:  ,  :,  0, :, 0] = 0.
            dc_nonlinear[:  ,  :, -1, :, 1] = 0.
            dc_nonlinear[:  ,  0,  :, :, 2] = 0.
            dc_nonlinear[:  , -1,  :, :, 3] = 0.
            dc_nonlinear[0  ,  :,  :, :, 4] = 0.
            dc_nonlinear[-1 ,  :,  :, :, 5] = 0.
            dc_nonlinear.shape = (nz*ny*nx*ng, 6)

        # Set the off-diagonals
        if nx > 1:
            dc_linear_data.append(-dc_linear[ng: , 0] * shape_nbr[ng:, 0])
            dc_linear_data.append(-dc_linear[:-ng, 1] * shape_nbr[:-ng, 1])

            if self.use_pcmfd:
                dc_nonlinear_data.append(-dc_nonlinear[ng: , 1] * shape_nbr[ng:, 0])
                dc_nonlinear_data.append(-dc_nonlinear[:-ng, 3] * shape_nbr[:-ng, 1])
            else:
                dc_nonlinear_data.append( dc_nonlinear[ng: , 0] * shape_nbr[ng:, 0])
                dc_nonlinear_data.append(-dc_nonlinear[:-ng, 1] * shape_nbr[:-ng, 1])

            diags.append(-ng)
            diags.append(ng)
        if ny > 1:
            dc_linear_data.append(-dc_linear[nx*ng: , 2] * shape_nbr[nx*ng:, 2])
            dc_linear_data.append(-dc_linear[:-nx*ng, 3] * shape_nbr[:-nx*ng, 3])

            if self.use_pcmfd:
                dc_nonlinear_data.append(-dc_nonlinear[nx*ng: , 5] * shape_nbr[nx*ng:, 2])
                dc_nonlinear_data.append(-dc_nonlinear[:-nx*ng, 7] * shape_nbr[:-nx*ng, 3])
            else:
                dc_nonlinear_data.append( dc_nonlinear[nx*ng: , 2] * shape_nbr[nx*ng:, 2])
                dc_nonlinear_data.append(-dc_nonlinear[:-nx*ng, 3] * shape_nbr[:-nx*ng, 3])

            diags.append(-nx*ng)
            diags.append(nx*ng)
        if nz > 1:
            dc_linear_data.append(-dc_linear[nx*ny*ng: , 4] * shape_nbr[ny*nx*ng:, 4])
            dc_linear_data.append(-dc_linear[:-nx*ny*ng, 5] * shape_nbr[:-ny*nx*ng, 5])

            if self.use_pcmfd:
                dc_nonlinear_data.append(-dc_nonlinear[nx*ny*ng: , 9] * shape_nbr[ny*nx*ng:, 4])
                dc_nonlinear_data.append(-dc_nonlinear[:-nx*ny*ng,11] * shape_nbr[:-ny*nx*ng, 5])
            else:
                dc_nonlinear_data.append( dc_nonlinear[nx*ny*ng: , 4] * shape_nbr[ny*nx*ng:, 4])
                dc_nonlinear_data.append(-dc_nonlinear[:-nx*ny*ng, 5] * shape_nbr[:-ny*nx*ng, 5])

            diags.append(-nx*ny*ng)
            diags.append(nx*ny*ng)

        return diags, dc_linear_data, dc_nonlinear_data


class InnerState(State):

    def __init__(self, states):
        super().__init__(states)

        # Initialize Solver class attributes
        self.fwd_state = states['FORWARD_OUTER']
        self.pre_state = states['PREVIOUS_OUTER']

    @property
    def fwd_state(self):
        return self._fwd_state

    @property
    def pre_state(self):
        return self._pre_state

    @fwd_state.setter
    def fwd_state(self, fwd_state):
        self._fwd_state = fwd_state

    @pre_state.setter
    def pre_state(self, pre_state):
        self._pre_state = pre_state

    @property
    def weight(self):
        time_point = self.clock.times[self.time_point]
        fwd_time = self.clock.times['FORWARD_OUTER']
        weight = 1 - (fwd_time - time_point) / self.clock.dt_outer
        return weight

    @property
    def inscatter(self):
        wgt = self.weight
        inscatter_fwd  = self.fwd_state.inscatter
        inscatter_prev = self.pre_state.inscatter
        inscatter = inscatter_fwd * wgt + inscatter_prev * (1 - wgt)
        inscatter[inscatter < 0.] = 0.
        return inscatter

    @property
    def outscatter(self):
        return self.inscatter.sum(axis=1)

    @property
    def absorption(self):
        wgt = self.weight
        absorption_fwd  = self.fwd_state.absorption
        absorption_prev = self.pre_state.absorption
        absorption = absorption_fwd * wgt + absorption_prev * (1 - wgt)
        absorption[absorption < 0.] = 0.
        return absorption

    @property
    def kappa_fission(self):
        wgt = self.weight
        kappa_fission_fwd  = self.fwd_state.kappa_fission
        kappa_fission_prev = self.pre_state.kappa_fission
        kappa_fission = kappa_fission_fwd * wgt + kappa_fission_prev * (1 - wgt)
        kappa_fission[kappa_fission < 0.] = 0.
        return kappa_fission

    @property
    def chi_prompt(self):
        wgt = self.weight
        chi_prompt_fwd  = self.fwd_state.chi_prompt
        chi_prompt_prev = self.pre_state.chi_prompt
        chi_prompt = chi_prompt_fwd * wgt + chi_prompt_prev * (1 - wgt)
        chi_prompt[chi_prompt < 0.] = 0.
        return chi_prompt

    @property
    def prompt_nu_fission(self):
        wgt = self.weight
        prompt_nu_fission_fwd  = self.fwd_state.prompt_nu_fission
        prompt_nu_fission_prev = self.pre_state.prompt_nu_fission
        prompt_nu_fission = prompt_nu_fission_fwd * wgt + prompt_nu_fission_prev * (1 - wgt)
        prompt_nu_fission[prompt_nu_fission < 0.] = 0.
        return prompt_nu_fission

    @property
    def chi_delayed(self):
        wgt = self.weight
        chi_delayed_fwd  = self.fwd_state.chi_delayed
        chi_delayed_prev = self.pre_state.chi_delayed
        chi_delayed = chi_delayed_fwd * wgt + chi_delayed_prev * (1 - wgt)
        chi_delayed[chi_delayed < 0.] = 0.
        return chi_delayed

    @property
    def delayed_nu_fission(self):
        wgt = self.weight
        delayed_nu_fission_fwd  = self.fwd_state.delayed_nu_fission
        delayed_nu_fission_prev = self.pre_state.delayed_nu_fission
        delayed_nu_fission = delayed_nu_fission_fwd * wgt + delayed_nu_fission_prev * (1 - wgt)
        delayed_nu_fission[delayed_nu_fission < 0.] = 0.
        return delayed_nu_fission

    @property
    def inverse_velocity(self):
        wgt = self.weight
        inverse_velocity_fwd  = self.fwd_state.inverse_velocity
        inverse_velocity_prev = self.pre_state.inverse_velocity
        inverse_velocity = inverse_velocity_fwd * wgt + inverse_velocity_prev * (1 - wgt)
        inverse_velocity[inverse_velocity < 0.] = 0.
        return inverse_velocity

    @property
    def decay_rate(self):
        wgt = self.weight
        decay_rate_fwd  = self.fwd_state.decay_rate
        decay_rate_prev = self.pre_state.decay_rate
        decay_rate = decay_rate_fwd * wgt + decay_rate_prev * (1 - wgt)
        decay_rate[decay_rate < 0.] = 0.
        decay_rate[decay_rate < 1.e-5] = 0.
        return decay_rate

    @property
    def shape(self):
        wgt = self.weight
        shape_fwd  = self.fwd_state.shape
        shape_prev = self.pre_state.shape
        shape = shape_fwd * wgt + shape_prev * (1 - wgt)
        shape[shape < 0.] = 0.
        return shape

    @property
    def shape_deriv(self):
        shape_fwd  = self.fwd_state.shape
        shape_prev = self.pre_state.shape
        deriv = (shape_fwd - shape_prev) / self.dt_outer
        return deriv

    @property
    def coupling_terms(self):
        wgt = self.weight
        diag, dc_lin_fwd, dc_nonlin_fwd = self.fwd_state.coupling_terms
        diag, dc_lin_pre, dc_nonlin_pre = self.pre_state.coupling_terms
        dc_lin = []
        dc_nonlin = []
        for i in range(len(dc_lin_fwd)):
            dc_lin   .append(dc_lin_fwd[i]    * wgt + dc_lin_pre[i]    * (1 - wgt))
            dc_nonlin.append(dc_nonlin_fwd[i] * wgt + dc_nonlin_pre[i] * (1 - wgt))
        return diag, dc_lin, dc_nonlin

    @property
    def time_removal_source(self):
        state = self.states['PREVIOUS_INNER']
        time_removal = self.inverse_velocity / self.dt_inner * self.shape
        time_removal = openmc.kinetics.map_array(time_removal, self.shape_zyxg, self.amplitude_zyxg, False)
        time_removal.shape = (self.amplitude_nxyz, self.ng)
        time_removal *= state.amplitude

        return time_removal

    @property
    def k1_source(self):
        state = self.states['PREVIOUS_INNER']
        decay_source = self.decay_rate * state.k1 * state.precursors
        decay_source = np.repeat(decay_source, self.ng)
        decay_source.shape = (self.shape_nxyz, self.nd, self.ng)
        decay_source *= self.chi_delayed
        decay_source = openmc.kinetics.map_array(decay_source.sum(axis=1), self.shape_zyxg, self.amplitude_zyxg, False)
        decay_source.shape = (self.amplitude_nxyz, self.ng)

        return decay_source

    @property
    def k3_source(self):
        state = self.states['PREVIOUS_INNER']
        k3_source = self.k3_source_matrix * state.amplitude.flatten()
        k3_source.shape = (self.amplitude_nxyz, self.ng)

        return k3_source

    @property
    def transient_source(self):
        k3_src = self.k3_source.flatten()
        k1_src = self.k1_source.flatten()
        time_src = self.time_removal_source.flatten()

        return time_src + k1_src - k3_src

    @property
    def transient_matrix(self):
        time_removal = self.inverse_velocity / self.dt_inner * self.shape
        time_removal = openmc.kinetics.map_array(time_removal, self.shape_zyxg, self.amplitude_zyxg, False)

        shape_deriv = self.inverse_velocity * self.shape_deriv
        shape_deriv = openmc.kinetics.map_array(shape_deriv, self.shape_zyxg, self.amplitude_zyxg, False)

        return sps.diags(time_removal.flatten()) + sps.diags(shape_deriv.flatten()) \
            - self.prompt_production_matrix + self.destruction_matrix \
            - self.k2_source_matrix

    @property
    def k1(self):
        return np.exp(- self.dt_inner * self.decay_rate)

    @property
    def k2(self):
        # Compute k2 / (lambda * k_crit)
        k2 = 1. - (1. - self.k1) / (self.dt_inner * self.decay_rate)
        return openmc.kinetics.nan_inf_to_zero(k2 / (self.decay_rate * self.k_crit))

    @property
    def k3(self):
        # Compute k3 / (lambda * k_crit)
        k3 = self.k1 - (1. - self.k1) / (self.dt_inner * self.decay_rate)
        return openmc.kinetics.nan_inf_to_zero(k3 / (self.decay_rate * self.k_crit))

    @property
    def k2_source_matrix(self):

        k2 = np.repeat(self.decay_rate * self.k2, self.ng * self.ng)
        k2.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)

        chi = np.repeat(self.chi_delayed, self.ng)
        chi.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)

        shape = np.tile(self.shape, self.nd)
        shape.shape = (self.shape_nxyz, self.nd, self.ng)
        del_fis_rate = self.delayed_nu_fission * shape
        del_fis_rate = np.tile(del_fis_rate, self.ng)
        del_fis_rate.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)

        term_k2 = (chi * k2 * del_fis_rate).sum(axis=1)
        term_k2 = openmc.kinetics.map_array(term_k2, self.shape_zyxgg, self.amplitude_zyxgg, False)
        term_k2.shape = (self.amplitude_nxyz, self.ng, self.ng)

        return openmc.kinetics.block_diag(term_k2)

    @property
    def k3_source_matrix(self):

        state = self.states['PREVIOUS_INNER']

        k3 = np.repeat(self.decay_rate * state.k3, self.ng * self.ng)
        k3.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)

        chi = np.repeat(self.chi_delayed, self.ng)
        chi.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)

        shape = np.tile(state.shape, self.nd)
        shape.shape = (self.shape_nxyz, self.nd, self.ng)
        del_fis_rate = state.delayed_nu_fission * shape
        del_fis_rate = np.tile(del_fis_rate, self.ng)
        del_fis_rate.shape = (self.shape_nxyz, self.nd, self.ng, self.ng)

        term_k3 = (chi * k3 * del_fis_rate).sum(axis=1)
        term_k3 = openmc.kinetics.map_array(term_k3, self.shape_zyxgg, self.amplitude_zyxgg, False)
        term_k3.shape = (self.amplitude_nxyz, self.ng, self.ng)

        return openmc.kinetics.block_diag(term_k3)

    @property
    def propagate_precursors(self):

        state = self.states['PREVIOUS_INNER']

        # Contribution from current precursors
        term_k1 = state.k1 * state.precursors

        # Contribution from generation at current time point
        flux = np.tile(self.flux, self.nd)
        flux.shape = (self.shape_nxyz, self.nd, self.ng)
        term_k2 = self.k2 * (self.delayed_nu_fission * flux).sum(axis=2)

        # Contribution from generation at previous time step
        flux = np.tile(state.flux, state.nd)
        flux.shape = (state.shape_nxyz, state.nd, state.ng)
        term_k3 = state.k3 * (state.delayed_nu_fission * flux).sum(axis=2)

        self._precursors = term_k1 + term_k2 - term_k3

    @property
    def dump_to_log_file(self):

        time_point = str(self.clock.times[self.time_point])
        f = h5py.File(self._log_file, 'a')
        if time_point not in f['INNER_STEPS'].keys():
            f['INNER_STEPS'].require_group(time_point)

        if 'amplitude' not in f['INNER_STEPS'][time_point].keys():
            f['INNER_STEPS'][time_point].create_dataset('amplitude', data=self.amplitude)
            f['INNER_STEPS'][time_point].create_dataset('flux', data=self.flux)
            f['INNER_STEPS'][time_point].create_dataset('power', data=self.power)
        else:
            amp = f['INNER_STEPS'][time_point]['amplitude']
            amp[...] = self.amplitude
            flux = f['INNER_STEPS'][time_point]['flux']
            flux[...] = self.flux
            power = f['INNER_STEPS'][time_point]['power']
            power[...] = self.power

        f['INNER_STEPS'][time_point].attrs['reactivity'] = self.reactivity
        f['INNER_STEPS'][time_point].attrs['beta_eff'] = self.beta_eff.sum()
        f['INNER_STEPS'][time_point].attrs['core_power_density'] = self.core_power_density

        f.close()
