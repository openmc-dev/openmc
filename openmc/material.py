from __future__ import annotations
from collections import defaultdict, namedtuple, Counter
from collections.abc import Iterable
from copy import deepcopy
from numbers import Real
from pathlib import Path
import re
import sys
import tempfile
from typing import Sequence, Dict
import warnings

import lxml.etree as ET
import numpy as np
import h5py

import openmc
import openmc.data
import openmc.checkvalue as cv
from ._xml import clean_indentation, get_elem_list, get_text
from .mixin import IDManagerMixin
from .utility_funcs import input_path
from . import waste
from openmc.checkvalue import PathLike
from openmc.stats import Univariate, Discrete, Mixture
from openmc.data.data import _get_element_symbol


# Units for density supported by OpenMC
DENSITY_UNITS = ('g/cm3', 'g/cc', 'kg/m3', 'atom/b-cm', 'atom/cm3', 'sum',
                 'macro')

# Smallest normalized floating point number
_SMALLEST_NORMAL = sys.float_info.min

_BECQUEREL_PER_CURIE = 3.7e10

NuclideTuple = namedtuple('NuclideTuple', ['name', 'percent', 'percent_type'])


class Material(IDManagerMixin):
    """A material composed of a collection of nuclides/elements.

    To create a material, one should create an instance of this class, add
    nuclides or elements with :meth:`Material.add_nuclide` or
    :meth:`Material.add_element`, respectively, and set the total material
    density with :meth:`Material.set_density()`. Alternatively, you can use
    :meth:`Material.add_components()` to pass a dictionary containing all the
    component information. The material can then be assigned to a cell using the
    :attr:`Cell.fill` attribute.

    Parameters
    ----------
    material_id : int, optional
        Unique identifier for the material. If not specified, an identifier will
        automatically be assigned.
    name : str, optional
        Name of the material. If not specified, the name will be the empty
        string.
    temperature : float, optional
        Temperature of the material in Kelvin. If not specified, the material
        inherits the default temperature applied to the model.

    Attributes
    ----------
    id : int
        Unique identifier for the material
    temperature : float
        Temperature of the material in Kelvin.
    density : float
        Density of the material (units defined separately)
    density_units : str
        Units used for `density`. Can be one of 'g/cm3', 'g/cc', 'kg/m3',
        'atom/b-cm', 'atom/cm3', 'sum', or 'macro'.  The 'macro' unit only
        applies in the case of a multi-group calculation.
    depletable : bool
        Indicate whether the material is depletable.
    nuclides : list of namedtuple
        List in which each item is a namedtuple consisting of a nuclide string,
        the percent density, and the percent type ('ao' or 'wo'). The namedtuple
        has field names ``name``, ``percent``, and ``percent_type``.
    isotropic : list of str
        Nuclides for which elastic scattering should be treated as though it
        were isotropic in the laboratory system.
    average_molar_mass : float
        The average molar mass of nuclides in the material in units of grams per
        mol.  For example, UO2 with 3 nuclides will have an average molar mass
        of 270 / 3 = 90 g / mol.
    volume : float
        Volume of the material in cm^3. This can either be set manually or
        calculated in a stochastic volume calculation and added via the
        :meth:`Material.add_volume_information` method.
    paths : list of str
        The paths traversed through the CSG tree to reach each material
        instance. This property is initialized by calling the
        :meth:`Geometry.determine_paths` method.
    num_instances : int
        The number of instances of this material throughout the geometry. This
        property is initialized by calling the :meth:`Geometry.determine_paths`
        method.
    fissionable_mass : float
        Mass of fissionable nuclides in the material in [g]. Requires that the
        :attr:`volume` attribute is set.
    ncrystal_cfg : str
        NCrystal configuration string

        .. versionadded:: 0.13.3

    """

    next_id = 1
    used_ids = set()

    def __init__(self, material_id=None, name='', temperature=None):
        # Initialize class attributes
        self.id = material_id
        self.name = name
        self.temperature = temperature
        self._density = None
        self._density_units = 'sum'
        self._depletable = False
        self._paths = None
        self._num_instances = None
        self._volume = None
        self._atoms = {}
        self._isotropic = []
        self._ncrystal_cfg = None

        # A list of tuples (nuclide, percent, percent type)
        self._nuclides = []

        # The single instance of Macroscopic data present in this material
        # (only one is allowed, hence this is different than _nuclides, etc)
        self._macroscopic = None

        # If specified, a list of table names
        self._sab = []

    def __repr__(self) -> str:
        string = 'Material\n'
        string += '{: <16}=\t{}\n'.format('\tID', self._id)
        string += '{: <16}=\t{}\n'.format('\tName', self._name)
        string += '{: <16}=\t{}\n'.format('\tTemperature', self._temperature)

        string += '{: <16}=\t{}'.format('\tDensity', self._density)
        string += f' [{self._density_units}]\n'

        string += '{: <16}=\t{} [cm^3]\n'.format('\tVolume', self._volume)
        string += '{: <16}=\t{}\n'.format('\tDepletable', self._depletable)

        string += '{: <16}\n'.format('\tS(a,b) Tables')

        if self._ncrystal_cfg:
            string += '{: <16}=\t{}\n'.format('\tNCrystal conf', self._ncrystal_cfg)

        for sab in self._sab:
            string += '{: <16}=\t{}\n'.format('\tS(a,b)', sab)

        string += '{: <16}\n'.format('\tNuclides')

        for nuclide, percent, percent_type in self._nuclides:
            string += '{: <16}'.format('\t{}'.format(nuclide))
            string += f'=\t{percent: <12} [{percent_type}]\n'

        if self._macroscopic is not None:
            string += '{: <16}\n'.format('\tMacroscopic Data')
            string += '{: <16}'.format('\t{}'.format(self._macroscopic))

        return string

    @property
    def name(self) -> str | None:
        return self._name

    @name.setter
    def name(self, name: str | None):
        if name is not None:
            cv.check_type(f'name for Material ID="{self._id}"',
                          name, str)
            self._name = name
        else:
            self._name = ''

    @property
    def temperature(self) -> float | None:
        return self._temperature

    @temperature.setter
    def temperature(self, temperature: Real | None):
        cv.check_type(f'Temperature for Material ID="{self._id}"',
                      temperature, (Real, type(None)))
        self._temperature = temperature

    @property
    def density(self) -> float | None:
        return self._density

    @property
    def density_units(self) -> str:
        return self._density_units

    @property
    def depletable(self) -> bool:
        return self._depletable

    @depletable.setter
    def depletable(self, depletable: bool):
        cv.check_type(f'Depletable flag for Material ID="{self._id}"',
                      depletable, bool)
        self._depletable = depletable

    @property
    def paths(self) -> list[str]:
        if self._paths is None:
            raise ValueError('Material instance paths have not been determined. '
                             'Call the Geometry.determine_paths() method.')
        return self._paths

    @property
    def num_instances(self) -> int:
        if self._num_instances is None:
            raise ValueError(
                'Number of material instances have not been determined. Call '
                'the Geometry.determine_paths() method.')
        return self._num_instances

    @property
    def nuclides(self) -> list[namedtuple]:
        return self._nuclides

    @property
    def isotropic(self) -> list[str]:
        return self._isotropic

    @isotropic.setter
    def isotropic(self, isotropic: Iterable[str]):
        cv.check_iterable_type('Isotropic scattering nuclides', isotropic,
                               str)
        self._isotropic = list(isotropic)

    @property
    def average_molar_mass(self) -> float:
        # Using the sum of specified atomic or weight amounts as a basis, sum
        # the mass and moles of the material
        mass = 0.
        moles = 0.
        for nuc in self.nuclides:
            if nuc.percent_type == 'ao':
                mass += nuc.percent * openmc.data.atomic_mass(nuc.name)
                moles += nuc.percent
            else:
                moles += nuc.percent / openmc.data.atomic_mass(nuc.name)
                mass += nuc.percent

        # Compute and return the molar mass
        return mass / moles

    @property
    def volume(self) -> float | None:
        return self._volume

    @volume.setter
    def volume(self, volume: Real):
        if volume is not None:
            cv.check_type('material volume', volume, Real)
        self._volume = volume

    @property
    def ncrystal_cfg(self) -> str | None:
        return self._ncrystal_cfg

    @property
    def fissionable_mass(self) -> float:
        if self.volume is None:
            raise ValueError("Volume must be set in order to determine mass.")
        density = 0.0
        for nuc, atoms_per_bcm in self.get_nuclide_atom_densities().items():
            Z = openmc.data.zam(nuc)[0]
            if Z >= 90:
                density += 1e24 * atoms_per_bcm * openmc.data.atomic_mass(nuc) \
                           / openmc.data.AVOGADRO
        return density*self.volume

    @property
    def decay_photon_energy(self) -> Univariate | None:
        warnings.warn(
            "The 'decay_photon_energy' property has been replaced by the "
            "get_decay_photon_energy() method and will be removed in a future "
            "version.", FutureWarning)
        return self.get_decay_photon_energy(0.0)

    def get_decay_photon_energy(
        self,
        clip_tolerance: float = 1e-6,
        units: str = 'Bq',
        volume: float | None = None,
        exclude_nuclides: list[str] | None = None,
        include_nuclides: list[str] | None = None
    ) -> Univariate | None:
        r"""Return energy distribution of decay photons from unstable nuclides.

        .. versionadded:: 0.14.0

        Parameters
        ----------
        clip_tolerance : float
            Maximum fraction of :math:`\sum_i x_i p_i` for discrete distributions
            that will be discarded.
        units : {'Bq', 'Bq/g', 'Bq/kg', 'Bq/cm3'}
            Specifies the units on the integral of the distribution.
        volume : float, optional
            Volume of the material. If not passed, defaults to using the
            :attr:`Material.volume` attribute.
        exclude_nuclides : list of str, optional
            Nuclides to exclude from the photon source calculation.
        include_nuclides : list of str, optional
            Nuclides to include in the photon source calculation. If specified,
            only these nuclides are used.

        Returns
        -------
        Univariate or None
            Decay photon energy distribution. The integral of this distribution is
            the total intensity of the photon source in the requested units.

        """
        cv.check_value('units', units, {'Bq', 'Bq/g', 'Bq/kg', 'Bq/cm3'})

        if exclude_nuclides is not None and include_nuclides is not None:
            raise ValueError("Cannot specify both exclude_nuclides and include_nuclides")

        if units == 'Bq':
            multiplier = volume if volume is not None else self.volume
            if multiplier is None:
                raise ValueError("volume must be specified if units='Bq'")
        elif units == 'Bq/cm3':
            multiplier = 1
        elif units == 'Bq/g':
            multiplier = 1.0 / self.get_mass_density()
        elif units == 'Bq/kg':
            multiplier = 1000.0 / self.get_mass_density()

        dists = []
        probs = []
        for nuc, atoms_per_bcm in self.get_nuclide_atom_densities().items():
            if exclude_nuclides is not None and nuc in exclude_nuclides:
                continue
            if include_nuclides is not None and nuc not in include_nuclides:
                continue

            source_per_atom = openmc.data.decay_photon_energy(nuc)
            if source_per_atom is not None and atoms_per_bcm > 0.0:
                dists.append(source_per_atom)
                probs.append(1e24 * atoms_per_bcm * multiplier)

        # If no photon sources, exit early
        if not dists:
            return None

        # Get combined distribution, clip low-intensity values in discrete spectra
        combined = openmc.data.combine_distributions(dists, probs)
        if isinstance(combined, (Discrete, Mixture)):
            combined.clip(clip_tolerance, inplace=True)

        # If clipping resulted in a single distribution within a mixture, pick
        # out that single distribution
        if isinstance(combined, Mixture) and len(combined.distribution) == 1:
            combined = combined.distribution[0]

        return combined


    def get_decay_photon_dose(
        self,
        n_samples: int = 10000,
        distance: float | None = None
    ) -> float:
        """Return an approximation of decay photon dose rate from unstable nuclides.

        Assumes a point source in a vacuum with no attenuation.

        .. versionadded:: 0.15.3

        Parameters
        ----------
        n_samples : int, optional
            Number of energy samples to use for Monte Carlo integration.
            For Discrete distributions, exact calculation is used instead of
            sampling. Defaults to 10000.
        distance : float, optional
            Distance from the source in cm at which to calculate dose rate.
            If not specified, the dose is calculated at the surface of an
            equivalent sphere based on the material's volume.

        Returns
        -------
        float
            Dose rate in pSv/s at the specified distance or at the surface
            of an equivalent sphere

        """
        cv.check_type('n_samples', n_samples, int)
        cv.check_greater_than('n_samples', n_samples, 0)
        if distance is not None:
            cv.check_type('distance', distance, Real)
            cv.check_greater_than('distance', distance, 0.0)

        # Get the decay photon energy distribution
        dist = self.get_decay_photon_energy(units='Bq')

        # If no photon sources, return zero dose
        if dist is None:
            return 0.0

        # Get dose coefficients for isotropic geometry
        energy_bins, dose_coeffs = openmc.data.dose_coefficients(
            particle="photon",
            geometry="AP"
        )

        # For pure Discrete distributions, compute mean dose coefficient exactly
        # without sampling to improve performance
        if isinstance(dist, Discrete):
            # Check if maximum energy exceeds the dose coefficient range
            max_energy = np.max(dist.x)
            if max_energy > energy_bins[-1]:
                raise ValueError(
                    f"Maximum photon energy {max_energy:.2e} eV exceeds "
                    f"dose coefficient upper limit of {energy_bins[-1]:.2e} eV"
                )

            # Interpolate dose coefficients at discrete energy points
            interpolated_coeffs = np.interp(dist.x, energy_bins, dose_coeffs)

            # Calculate weighted mean dose coefficient
            mean_dose_coeff = np.sum(interpolated_coeffs * dist.p) / dist.integral()

        else:
            # For other distributions (Mixture, etc.), use Monte Carlo sampling
            # Sample energies from the distribution with fixed seed for reproducibility
            energies = dist.sample(n_samples=n_samples, seed=1)

            # Check if maximum sampled energy exceeds the dose coefficient range
            max_energy = np.max(energies)
            if max_energy > energy_bins[-1]:
                raise ValueError(
                    f"Maximum sampled photon energy {max_energy:.2e} eV exceeds "
                    f"dose coefficient upper limit of {energy_bins[-1]:.2e} eV"
                )

            # Interpolate dose coefficients at sampled energies
            # dose_coeffs are in pSv·cm²
            interpolated_coeffs = np.interp(energies, energy_bins, dose_coeffs)

            # Calculate mean dose coefficient
            mean_dose_coeff = np.mean(interpolated_coeffs)

        # Total activity in Bq (photons/s)
        total_activity = dist.integral()

        # If distance not specified, use surface of equivalent sphere
        if distance is None:
            if self.volume is None:
                raise ValueError(
                    "Either distance must be specified or material volume "
                    "must be set to calculate dose rate."
                )
            # Calculate radius of equivalent sphere: V = (4/3)πr³
            # r = (3V/4π)^(1/3)
            distance = (3.0 * self.volume / (4.0 * np.pi)) ** (1.0 / 3.0)
        # Dose rate at specified distance in pSv/s
        # dose_rate = (coeff * activity) / (4π·r²)
        dose_rate_pSv_per_s = mean_dose_coeff * total_activity / (4 * np.pi * distance**2)

        return dose_rate_pSv_per_s


    @classmethod
    def from_hdf5(cls, group: h5py.Group) -> Material:
        """Create material from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file

        Returns
        -------
        openmc.Material
            Material instance

        """
        mat_id = int(group.name.split('/')[-1].lstrip('material '))

        name = group['name'][()].decode() if 'name' in group else ''
        density = group['atom_density'][()]
        if 'nuclide_densities' in group:
            nuc_densities = group['nuclide_densities'][()]

        # Create the Material
        material = cls(mat_id, name)
        material.depletable = bool(group.attrs['depletable'])
        if 'volume' in group.attrs:
            material.volume = group.attrs['volume']
        if "temperature" in group.attrs:
            material.temperature = group.attrs["temperature"]

        # Read the names of the S(a,b) tables for this Material and add them
        if 'sab_names' in group:
            sab_tables = group['sab_names'][()]
            for sab_table in sab_tables:
                name = sab_table.decode()
                material.add_s_alpha_beta(name)

        # Set the Material's density to atom/b-cm as used by OpenMC
        material.set_density(density=density, units='atom/b-cm')

        if 'nuclides' in group:
            nuclides = group['nuclides'][()]
            # Add all nuclides to the Material
            for fullname, density in zip(nuclides, nuc_densities):
                name = fullname.decode().strip()
                material.add_nuclide(name, percent=density, percent_type='ao')
        if 'macroscopics' in group:
            macroscopics = group['macroscopics'][()]
            # Add all macroscopics to the Material
            for fullname in macroscopics:
                name = fullname.decode().strip()
                material.add_macroscopic(name)

        return material

    @classmethod
    def from_ncrystal(cls, cfg, **kwargs) -> Material:
        """Create material from NCrystal configuration string.

        Density, temperature, and material composition, and (ultimately) thermal
        neutron scattering will be automatically be provided by NCrystal based
        on this string. The name and material_id parameters are simply passed on
        to the Material constructor.

        .. versionadded:: 0.13.3

        Parameters
        ----------
        cfg : str
            NCrystal configuration string
        **kwargs
            Keyword arguments passed to :class:`openmc.Material`

        Returns
        -------
        openmc.Material
            Material instance

        """

        try:
            import NCrystal
        except ModuleNotFoundError as e:
            raise RuntimeError('The .from_ncrystal method requires'
                               ' NCrystal to be installed.') from e
        nc_mat = NCrystal.createInfo(cfg)

        def openmc_natabund(Z):
            #nc_mat.getFlattenedComposition might need natural abundancies.
            #This call-back function is used so NCrystal can flatten composition
            #using OpenMC's natural abundancies. In practice this function will
            #only get invoked in the unlikely case where a material is specified
            #by referring both to natural elements and specific isotopes of the
            #same element.
            elem_name = openmc.data.ATOMIC_SYMBOL[Z]
            return [
                (int(iso_name[len(elem_name):]), abund)
                for iso_name, abund in openmc.data.isotopes(elem_name)
            ]

        flat_compos = nc_mat.getFlattenedComposition(
            preferNaturalElements=True, naturalAbundProvider=openmc_natabund)

        # Create the Material
        material = cls(temperature=nc_mat.getTemperature(), **kwargs)

        for Z, A_vals in flat_compos:
            elemname = openmc.data.ATOMIC_SYMBOL[Z]
            for A, frac in A_vals:
                if A:
                    material.add_nuclide(f'{elemname}{A}', frac)
                else:
                    material.add_element(elemname, frac)

        material.set_density('g/cm3', nc_mat.getDensity())
        material._ncrystal_cfg = NCrystal.normaliseCfg(cfg)

        return material

    def add_volume_information(self, volume_calc):
        """Add volume information to a material.

        Parameters
        ----------
        volume_calc : openmc.VolumeCalculation
            Results from a stochastic volume calculation

        """
        if volume_calc.domain_type == 'material':
            if self.id in volume_calc.volumes:
                self._volume = volume_calc.volumes[self.id].n
                self._atoms = volume_calc.atoms[self.id]
            else:
                raise ValueError('No volume information found for material ID={}.'
                                 .format(self.id))
        else:
            raise ValueError(f'No volume information found for material ID={self.id}.')

    def set_density(self, units: str, density: float | None = None):
        """Set the density of the material

        Parameters
        ----------
        units : {'g/cm3', 'g/cc', 'kg/m3', 'atom/b-cm', 'atom/cm3', 'sum', 'macro'}
            Physical units of density.
        density : float, optional
            Value of the density. Must be specified unless units is given as
            'sum'.

        """

        cv.check_value('density units', units, DENSITY_UNITS)
        self._density_units = units

        if units == 'sum':
            if density is not None:
                msg = 'Density "{}" for Material ID="{}" is ignored ' \
                      'because the unit is "sum"'.format(density, self.id)
                warnings.warn(msg)
        else:
            if density is None:
                msg = 'Unable to set the density for Material ID="{}" ' \
                      'because a density value must be given when not using ' \
                      '"sum" unit'.format(self.id)
                raise ValueError(msg)

            cv.check_type(f'the density for Material ID="{self.id}"',
                          density, Real)
            self._density = density

    def add_nuclide(self, nuclide: str, percent: float, percent_type: str = 'ao'):
        """Add a nuclide to the material

        Parameters
        ----------
        nuclide : str
            Nuclide to add, e.g., 'Mo95'
        percent : float
            Atom or weight percent
        percent_type : {'ao', 'wo'}
            'ao' for atom percent and 'wo' for weight percent

        """
        cv.check_type('nuclide', nuclide, str)
        cv.check_type('percent', percent, Real)
        cv.check_value('percent type', percent_type, {'ao', 'wo'})
        cv.check_greater_than('percent', percent, 0, equality=True)

        if self._macroscopic is not None:
            msg = 'Unable to add a Nuclide to Material ID="{}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if self._ncrystal_cfg is not None:
            raise ValueError("Cannot add nuclides to NCrystal material")

        # If nuclide name doesn't look valid, give a warning
        try:
            Z, _, _ = openmc.data.zam(nuclide)
        except ValueError as e:
            warnings.warn(str(e))
        else:
            # For actinides, have the material be depletable by default
            if Z >= 89:
                self.depletable = True

        self._nuclides.append(NuclideTuple(nuclide, percent, percent_type))

    def add_components(self, components: dict, percent_type: str = 'ao'):
        """ Add multiple elements or nuclides to a material

        .. versionadded:: 0.13.1

        Parameters
        ----------
        components : dict of str to float or dict
            Dictionary mapping element or nuclide names to their atom or weight
            percent. To specify enrichment of an element, the entry of
            ``components`` for that element must instead be a dictionary
            containing the keyword arguments as well as a value for
            ``'percent'``
        percent_type : {'ao', 'wo'}
            'ao' for atom percent and 'wo' for weight percent

        Examples
        --------
        >>> mat = openmc.Material()
        >>> components  = {'Li': {'percent': 1.0,
        >>>                       'enrichment': 60.0,
        >>>                       'enrichment_target': 'Li7'},
        >>>                'Fl': 1.0,
        >>>                'Be6': 0.5}
        >>> mat.add_components(components)

        """

        for component, params in components.items():
            cv.check_type('component', component, str)
            if isinstance(params, Real):
                params = {'percent': params}

            else:
                cv.check_type('params', params, dict)
                if 'percent' not in params:
                    raise ValueError("An entry in the dictionary does not have "
                                     "a required key: 'percent'")

            params['percent_type'] = percent_type

            # check if nuclide
            if not component.isalpha():
                self.add_nuclide(component, **params)
            else:
                self.add_element(component, **params)

    def remove_nuclide(self, nuclide: str):
        """Remove a nuclide from the material

        Parameters
        ----------
        nuclide : str
            Nuclide to remove

        """
        cv.check_type('nuclide', nuclide, str)

        # If the Material contains the Nuclide, delete it
        for nuc in reversed(self.nuclides):
            if nuclide == nuc.name:
                self.nuclides.remove(nuc)

    def remove_element(self, element):
        """Remove an element from the material

        .. versionadded:: 0.13.1

        Parameters
        ----------
        element : str
            Element to remove

        """
        cv.check_type('element', element, str)

        # If the Material contains the element, delete it
        for nuc in reversed(self.nuclides):
            element_name = re.split(r'\d+', nuc.name)[0]
            if element_name == element:
                self.nuclides.remove(nuc)

    def add_macroscopic(self, macroscopic: str):
        """Add a macroscopic to the material.  This will also set the
        density of the material to 1.0, unless it has been otherwise set,
        as a default for Macroscopic cross sections.

        Parameters
        ----------
        macroscopic : str
            Macroscopic to add

        """

        # Ensure no nuclides, elements, or sab are added since these would be
        # incompatible with macroscopics
        if self._nuclides or self._sab:
            msg = 'Unable to add a Macroscopic data set to Material ID="{}" ' \
                  'with a macroscopic value "{}" as an incompatible data ' \
                  'member (i.e., nuclide or S(a,b) table) ' \
                  'has already been added'.format(self._id, macroscopic)
            raise ValueError(msg)

        if not isinstance(macroscopic, str):
            msg = 'Unable to add a Macroscopic to Material ID="{}" with a ' \
                  'non-string value "{}"'.format(self._id, macroscopic)
            raise ValueError(msg)

        if self._macroscopic is None:
            self._macroscopic = macroscopic
        else:
            msg = 'Unable to add a Macroscopic to Material ID="{}". ' \
                  'Only one Macroscopic allowed per ' \
                  'Material.'.format(self._id)
            raise ValueError(msg)

        # Generally speaking, the density for a macroscopic object will
        # be 1.0. Therefore, lets set density to 1.0 so that the user
        # doesn't need to set it unless its needed.
        # Of course, if the user has already set a value of density,
        # then we will not override it.
        if self._density is None:
            self.set_density('macro', 1.0)

    def remove_macroscopic(self, macroscopic: str):
        """Remove a macroscopic from the material

        Parameters
        ----------
        macroscopic : str
            Macroscopic to remove

        """

        if not isinstance(macroscopic, str):
            msg = 'Unable to remove a Macroscopic "{}" in Material ID="{}" ' \
                  'since it is not a string'.format(self._id, macroscopic)
            raise ValueError(msg)

        # If the Material contains the Macroscopic, delete it
        if macroscopic == self._macroscopic:
            self._macroscopic = None

    def add_element(self, element: str, percent: float, percent_type: str = 'ao',
                    enrichment: float | None = None,
                    enrichment_target: str | None = None,
                    enrichment_type: str | None = None,
                    cross_sections: str | None = None):
        """Add a natural element to the material

        Parameters
        ----------
        element : str
            Element to add, e.g., 'Zr' or 'Zirconium'
        percent : float
            Atom or weight percent
        percent_type : {'ao', 'wo'}, optional
            'ao' for atom percent and 'wo' for weight percent. Defaults to atom
            percent.
        enrichment : float, optional
            Enrichment of an enrichment_target nuclide in percent (ao or wo).
            If enrichment_target is not supplied then it is enrichment for U235
            in weight percent. For example, input 4.95 for 4.95 weight percent
            enriched U.
            Default is None (natural composition).
        enrichment_target: str, optional
            Single nuclide name to enrich from a natural composition (e.g., 'O16')

            .. versionadded:: 0.12
        enrichment_type: {'ao', 'wo'}, optional
            'ao' for enrichment as atom percent and 'wo' for weight percent.
            Default is: 'ao' for two-isotope enrichment; 'wo' for U enrichment

            .. versionadded:: 0.12
        cross_sections : str, optional
            Location of cross_sections.xml file.

        Notes
        -----
        General enrichment procedure is allowed only for elements composed of
        two isotopes. If `enrichment_target` is given without `enrichment`
        natural composition is added to the material.

        """

        cv.check_type('nuclide', element, str)
        cv.check_type('percent', percent, Real)
        cv.check_greater_than('percent', percent, 0, equality=True)
        cv.check_value('percent type', percent_type, {'ao', 'wo'})

        # Make sure element name is just that
        if not element.isalpha():
            raise ValueError("Element name should be given by the "
                             "element's symbol or name, e.g., 'Zr', 'zirconium'")

        if self._ncrystal_cfg is not None:
            raise ValueError("Cannot add elements to NCrystal material")

        # Allow for element identifier to be given as a symbol or name
        if len(element) > 2:
            el = element.lower()
            element = openmc.data.ELEMENT_SYMBOL.get(el)
            if element is None:
                msg = f'Element name "{el}" not recognised'
                raise ValueError(msg)
        else:
            if element[0].islower():
                msg = f'Element name "{element}" should start with an uppercase letter'
                raise ValueError(msg)
            if len(element) == 2 and element[1].isupper():
                msg = f'Element name "{element}" should end with a lowercase letter'
                raise ValueError(msg)
            # skips the first entry of ATOMIC_SYMBOL which is n for neutron
            if element not in list(openmc.data.ATOMIC_SYMBOL.values())[1:]:
                msg = f'Element name "{element}" not recognised'
                raise ValueError(msg)

        if self._macroscopic is not None:
            msg = 'Unable to add an Element to Material ID="{}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if enrichment is not None and enrichment_target is None:
            if not isinstance(enrichment, Real):
                msg = 'Unable to add an Element to Material ID="{}" with a ' \
                      'non-floating point enrichment value "{}"'\
                      .format(self._id, enrichment)
                raise ValueError(msg)

            elif element != 'U':
                msg = 'Unable to use enrichment for element {} which is not ' \
                      'uranium for Material ID="{}"'.format(element, self._id)
                raise ValueError(msg)

            # Check that the enrichment is in the valid range
            cv.check_less_than('enrichment', enrichment, 100./1.008)
            cv.check_greater_than('enrichment', enrichment, 0., equality=True)

            if enrichment > 5.0:
                msg = 'A uranium enrichment of {} was given for Material ID='\
                      '"{}". OpenMC assumes the U234/U235 mass ratio is '\
                      'constant at 0.008, which is only valid at low ' \
                      'enrichments. Consider setting the isotopic ' \
                      'composition manually for enrichments over 5%.'.\
                      format(enrichment, self._id)
                warnings.warn(msg)

        # Add naturally-occuring isotopes
        element = openmc.Element(element)
        for nuclide in element.expand(percent,
                                      percent_type,
                                      enrichment,
                                      enrichment_target,
                                      enrichment_type,
                                      cross_sections):
            self.add_nuclide(*nuclide)

    def add_elements_from_formula(self, formula: str, percent_type: str = 'ao',
                                  enrichment: float | None = None,
                                  enrichment_target: str | None = None,
                                  enrichment_type: str | None = None):
        """Add a elements from a chemical formula to the material.

        .. versionadded:: 0.12

        Parameters
        ----------
        formula : str
            Formula to add, e.g., 'C2O', 'C6H12O6', or (NH4)2SO4.
            Note this is case sensitive, elements must start with an uppercase
            character. Multiplier numbers must be integers.
        percent_type : {'ao', 'wo'}, optional
            'ao' for atom percent and 'wo' for weight percent. Defaults to atom
            percent.
        enrichment : float, optional
            Enrichment of an enrichment_target nuclide in percent (ao or wo).
            If enrichment_target is not supplied then it is enrichment for U235
            in weight percent. For example, input 4.95 for 4.95 weight percent
            enriched U. Default is None (natural composition).
        enrichment_target : str, optional
            Single nuclide name to enrich from a natural composition (e.g., 'O16')
        enrichment_type : {'ao', 'wo'}, optional
            'ao' for enrichment as atom percent and 'wo' for weight percent.
            Default is: 'ao' for two-isotope enrichment; 'wo' for U enrichment

        Notes
        -----
        General enrichment procedure is allowed only for elements composed of
        two isotopes. If `enrichment_target` is given without `enrichment`
        natural composition is added to the material.

        """
        cv.check_type('formula', formula, str)

        if '.' in formula:
            msg = 'Non-integer multiplier values are not accepted. The ' \
                  'input formula {} contains a "." character.'.format(formula)
            raise ValueError(msg)

        # Tokenizes the formula and check validity of tokens
        tokens = re.findall(r"([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)", formula)
        for row in tokens:
            for token in row:
                if token.isalpha():
                    if token == "n" or token not in openmc.data.ATOMIC_NUMBER:
                        msg = f'Formula entry {token} not an element symbol.'
                        raise ValueError(msg)
                elif token not in ['(', ')', ''] and not token.isdigit():
                        msg = 'Formula must be made from a sequence of ' \
                              'element symbols, integers, and brackets. ' \
                              '{} is not an allowable entry.'.format(token)
                        raise ValueError(msg)

        # Checks that the number of opening and closing brackets are equal
        if formula.count('(') != formula.count(')'):
            msg = 'Number of opening and closing brackets is not equal ' \
                  'in the input formula {}.'.format(formula)
            raise ValueError(msg)

        # Checks that every part of the original formula has been tokenized
        for row in tokens:
            for token in row:
                formula = formula.replace(token, '', 1)
        if len(formula) != 0:
            msg = 'Part of formula was not successfully parsed as an ' \
                  'element symbol, bracket or integer. {} was not parsed.' \
                  .format(formula)
            raise ValueError(msg)

        # Works through the tokens building a stack
        mat_stack = [Counter()]
        for symbol, multi1, opening_bracket, closing_bracket, multi2 in tokens:
            if symbol:
                mat_stack[-1][symbol] += int(multi1 or 1)
            if opening_bracket:
                mat_stack.append(Counter())
            if closing_bracket:
                stack_top = mat_stack.pop()
                for symbol, value in stack_top.items():
                    mat_stack[-1][symbol] += int(multi2 or 1) * value

        # Normalizing percentages
        percents = mat_stack[0].values()
        norm_percents = [float(i) / sum(percents) for i in percents]
        elements = mat_stack[0].keys()

        # Adds each element and percent to the material
        for element, percent in zip(elements, norm_percents):
            if enrichment_target is not None and element == re.sub(r'\d+$', '', enrichment_target):
                self.add_element(element, percent, percent_type, enrichment,
                                 enrichment_target, enrichment_type)
            elif enrichment is not None and enrichment_target is None and element == 'U':
                self.add_element(element, percent, percent_type, enrichment)
            else:
                self.add_element(element, percent, percent_type)

    def add_s_alpha_beta(self, name: str, fraction: float = 1.0):
        r"""Add an :math:`S(\alpha,\beta)` table to the material

        Parameters
        ----------
        name : str
            Name of the :math:`S(\alpha,\beta)` table
        fraction : float
            The fraction of relevant nuclei that are affected by the
            :math:`S(\alpha,\beta)` table.  For example, if the material is a
            block of carbon that is 60% graphite and 40% amorphous then add a
            graphite :math:`S(\alpha,\beta)` table with fraction=0.6.

        """

        if self._macroscopic is not None:
            msg = 'Unable to add an S(a,b) table to Material ID="{}" as a ' \
                  'macroscopic data-set has already been added'.format(self._id)
            raise ValueError(msg)

        if not isinstance(name, str):
            msg = 'Unable to add an S(a,b) table to Material ID="{}" with a ' \
                        'non-string table name "{}"'.format(self._id, name)
            raise ValueError(msg)

        cv.check_type('S(a,b) fraction', fraction, Real)
        cv.check_greater_than('S(a,b) fraction', fraction, 0.0, True)
        cv.check_less_than('S(a,b) fraction', fraction, 1.0, True)
        self._sab.append((name, fraction))

    def make_isotropic_in_lab(self):
        self.isotropic = [x.name for x in self._nuclides]

    def get_elements(self) -> list[str]:
        """Returns all elements in the material

        .. versionadded:: 0.12

        Returns
        -------
        elements : list of str
            List of element names

        """

        return sorted({re.split(r'(\d+)', i)[0] for i in self.get_nuclides()})

    def get_nuclides(self, element: str | None = None) -> list[str]:
        """Returns a list of all nuclides in the material, if the element
        argument is specified then just nuclides of that element are returned.

        Parameters
        ----------
        element : str
            Specifies the element to match when searching through the nuclides

            .. versionadded:: 0.13.2

        Returns
        -------
        nuclides : list of str
            List of nuclide names
        """

        matching_nuclides = []
        if element:
            for nuclide in self._nuclides:
                if re.split(r'(\d+)', nuclide.name)[0] == element:
                    if nuclide.name not in matching_nuclides:
                        matching_nuclides.append(nuclide.name)
        else:
            for nuclide in self._nuclides:
                if nuclide.name not in matching_nuclides:
                    matching_nuclides.append(nuclide.name)

        return matching_nuclides

    def get_nuclide_densities(self) -> dict[str, tuple]:
        """Returns all nuclides in the material and their densities

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are 3-tuples of
            (nuclide, density percent, density percent type)

        """

        nuclides = {}

        for nuclide in self._nuclides:
            nuclides[nuclide.name] = nuclide

        return nuclides

    def get_nuclide_atom_densities(self, nuclide: str | None = None) -> dict[str, float]:
        """Returns one or all nuclides in the material and their atomic
        densities in units of atom/b-cm

        .. versionchanged:: 0.13.1
            The values in the dictionary were changed from a tuple containing
            the nuclide name and the density to just the density.

        Parameters
        ----------
        nuclides : str, optional
            Nuclide for which atom density is desired. If not specified, the
            atom density for each nuclide in the material is given.

            .. versionadded:: 0.13.2

        Returns
        -------
        nuclides : dict
            Dictionary whose keys are nuclide names and values are densities in
            [atom/b-cm]

        """

        sum_density = False
        if self.density_units == 'sum':
            sum_density = True
            density = 0.
        elif self.density_units == 'macro':
            density = self.density
        elif self.density_units == 'g/cc' or self.density_units == 'g/cm3':
            density = -self.density
        elif self.density_units == 'kg/m3':
            density = -0.001 * self.density
        elif self.density_units == 'atom/b-cm':
            density = self.density
        elif self.density_units == 'atom/cm3' or self.density_units == 'atom/cc':
            density = 1.e-24 * self.density

        # For ease of processing split out nuc, nuc_density,
        # and nuc_density_type into separate arrays
        nucs = []
        nuc_densities = []
        nuc_density_types = []

        for nuc in self.nuclides:
            nucs.append(nuc.name)
            nuc_densities.append(nuc.percent)
            nuc_density_types.append(nuc.percent_type)

        nuc_densities = np.array(nuc_densities)
        nuc_density_types = np.array(nuc_density_types)

        if sum_density:
            density = np.sum(nuc_densities)

        percent_in_atom = np.all(nuc_density_types == 'ao')
        density_in_atom = density > 0.
        sum_percent = 0.

        # Convert the weight amounts to atomic amounts
        if not percent_in_atom:
            for n, nuc in enumerate(nucs):
                nuc_densities[n] *= self.average_molar_mass / \
                                    openmc.data.atomic_mass(nuc)

        # Now that we have the atomic amounts, lets finish calculating densities
        sum_percent = np.sum(nuc_densities)
        nuc_densities = nuc_densities / sum_percent

        # Convert the mass density to an atom density
        if not density_in_atom:
            density = -density / self.average_molar_mass * 1.e-24 \
                      * openmc.data.AVOGADRO

        nuc_densities = density * nuc_densities

        nuclides = {}
        for n, nuc in enumerate(nucs):
            if nuclide is None or nuclide == nuc:
                nuclides[nuc] = nuc_densities[n]

        return nuclides

    def get_element_atom_densities(self, element: str | None = None) -> dict[str, float]:
        """Returns one or all elements in the material and their atomic
        densities in units of atom/b-cm

        .. versionadded:: 0.15.1

        Parameters
        ----------
        element : str, optional
            Element for which atom density is desired. If not specified, the
            atom density for each element in the material is given.

        Returns
        -------
        elements : dict
            Dictionary whose keys are element names and values are densities in
            [atom/b-cm]

        """
        if element is not None:
            element = _get_element_symbol(element)

        nuc_densities = self.get_nuclide_atom_densities()

        # Initialize an empty dictionary for summed values
        densities = {}

        # Accumulate densities for each nuclide
        for nuclide, density in nuc_densities.items():
            nuc_element = openmc.data.ATOMIC_SYMBOL[openmc.data.zam(nuclide)[0]]
            if element is None or element == nuc_element:
                if nuc_element not in densities:
                    densities[nuc_element] = 0.0
                densities[nuc_element] += float(density)

        # If specific element was requested, make sure it is present
        if element is not None and element not in densities:
                raise ValueError(f'Element {element} not found in material.')

        return densities


    def get_activity(self, units: str = 'Bq/cm3', by_nuclide: bool = False,
                     volume: float | None = None) -> dict[str, float] | float:
        """Returns the activity of the material or of each nuclide within.

        .. versionadded:: 0.13.1

        Parameters
        ----------
        units : {'Bq', 'Bq/g', 'Bq/kg', 'Bq/cm3', 'Ci', 'Ci/m3'}
            Specifies the type of activity to return, options include total
            activity [Bq,Ci], specific [Bq/g, Bq/kg] or volumetric activity
            [Bq/cm3,Ci/m3]. Default is volumetric activity [Bq/cm3].
        by_nuclide : bool
            Specifies if the activity should be returned for the material as a
            whole or per nuclide. Default is False.
        volume : float, optional
            Volume of the material. If not passed, defaults to using the
            :attr:`Material.volume` attribute.

            .. versionadded:: 0.13.3

        Returns
        -------
        Union[dict, float]
            If by_nuclide is True then a dictionary whose keys are nuclide
            names and values are activity is returned. Otherwise the activity
            of the material is returned as a float.
        """

        cv.check_value('units', units, {'Bq', 'Bq/g', 'Bq/kg', 'Bq/cm3', 'Ci', 'Ci/m3'})
        cv.check_type('by_nuclide', by_nuclide, bool)

        if volume is None:
            volume = self.volume

        if units == 'Bq':
            multiplier = volume
        elif units == 'Bq/cm3':
            multiplier = 1
        elif units == 'Bq/g':
            multiplier = 1.0 / self.get_mass_density()
        elif units == 'Bq/kg':
            multiplier = 1000.0 / self.get_mass_density()
        elif units == 'Ci':
            multiplier = volume / _BECQUEREL_PER_CURIE
        elif units == 'Ci/m3':
            multiplier = 1e6 / _BECQUEREL_PER_CURIE

        activity = {}
        for nuclide, atoms_per_bcm in self.get_nuclide_atom_densities().items():
            inv_seconds = openmc.data.decay_constant(nuclide)
            activity[nuclide] = inv_seconds * 1e24 * atoms_per_bcm * multiplier

        return activity if by_nuclide else sum(activity.values())

    def get_decay_heat(self, units: str = 'W', by_nuclide: bool = False,
                       volume: float | None = None) -> dict[str, float] | float:
        """Returns the decay heat of the material or for each nuclide in the
        material in units of [W], [W/g], [W/kg] or [W/cm3].

        .. versionadded:: 0.13.3

        Parameters
        ----------
        units : {'W', 'W/g', 'W/kg', 'W/cm3'}
            Specifies the units of decay heat to return. Options include total
            heat [W], specific [W/g, W/kg] or volumetric heat [W/cm3].
            Default is total heat [W].
        by_nuclide : bool
            Specifies if the decay heat should be returned for the material as a
            whole or per nuclide. Default is False.
        volume : float, optional
            Volume of the material. If not passed, defaults to using the
            :attr:`Material.volume` attribute.

            .. versionadded:: 0.13.3

        Returns
        -------
        Union[dict, float]
            If `by_nuclide` is True then a dictionary whose keys are nuclide
            names and values are decay heat is returned. Otherwise the decay heat
            of the material is returned as a float.
        """

        cv.check_value('units', units, {'W', 'W/g', 'W/kg', 'W/cm3'})
        cv.check_type('by_nuclide', by_nuclide, bool)

        if units == 'W':
            multiplier = volume if volume is not None else self.volume
        elif units == 'W/cm3':
            multiplier = 1
        elif units == 'W/g':
            multiplier = 1.0 / self.get_mass_density()
        elif units == 'W/kg':
            multiplier = 1000.0 / self.get_mass_density()

        decayheat = {}
        for nuclide, atoms_per_bcm in self.get_nuclide_atom_densities().items():
            decay_erg = openmc.data.decay_energy(nuclide)
            inv_seconds = openmc.data.decay_constant(nuclide)
            decay_erg *= openmc.data.JOULE_PER_EV
            decayheat[nuclide] = inv_seconds * decay_erg * 1e24 * atoms_per_bcm * multiplier

        return decayheat if by_nuclide else sum(decayheat.values())

    def get_nuclide_atoms(self, volume: float | None = None) -> dict[str, float]:
        """Return number of atoms of each nuclide in the material

        .. versionadded:: 0.13.1

        Parameters
        ----------
        volume : float, optional
            Volume of the material. If not passed, defaults to using the
            :attr:`Material.volume` attribute.

            .. versionadded:: 0.13.3

        Returns
        -------
        dict
            Dictionary whose keys are nuclide names and values are number of
            atoms present in the material.

        """
        if volume is None:
            volume = self.volume
        if volume is None:
            raise ValueError("Volume must be set in order to determine atoms.")
        atoms = {}
        for nuclide, atom_per_bcm in self.get_nuclide_atom_densities().items():
            atoms[nuclide] = 1.0e24 * atom_per_bcm * volume
        return atoms

    def get_mass_density(self, nuclide: str | None = None) -> float:
        """Return mass density of one or all nuclides

        Parameters
        ----------
        nuclides : str, optional
            Nuclide for which density is desired. If not specified, the density
            for the entire material is given.

        Returns
        -------
        float
            Density of the nuclide/material in [g/cm^3]

        """
        mass_density = 0.0
        for nuc, atoms_per_bcm in self.get_nuclide_atom_densities(nuclide=nuclide).items():
            density_i = 1e24 * atoms_per_bcm * openmc.data.atomic_mass(nuc) \
                        / openmc.data.AVOGADRO
            mass_density += density_i
        return mass_density

    def get_mass(self, nuclide: str | None = None, volume: float | None = None) -> float:
        """Return mass of one or all nuclides.

        Note that this method requires that the :attr:`Material.volume` has
        already been set.

        Parameters
        ----------
        nuclides : str, optional
            Nuclide for which mass is desired. If not specified, the density
            for the entire material is given.
        volume : float, optional
            Volume of the material. If not passed, defaults to using the
            :attr:`Material.volume` attribute.

            .. versionadded:: 0.13.3


        Returns
        -------
        float
            Mass of the nuclide/material in [g]

        """
        if volume is None:
            volume = self.volume
        if volume is None:
            raise ValueError("Volume must be set in order to determine mass.")
        return volume*self.get_mass_density(nuclide)

    def waste_classification(self, metal: bool = False) -> str:
        """Classify the material for near-surface waste disposal.

        This method determines a waste classification for the material based on
        the NRC regulations (10 CFR 61.55). Note that the NRC regulations do not
        consider many long-lived radionuclides relevant to fusion systems; for
        fusion applications, it is recommended to calculate a waste disposal
        rating based on limits by Fetter et al. using the
        :meth:`~openmc.Material.waste_disposal_rating` method.

        Parameters
        ----------
        metal : bool, optional
            Whether or not the material is in metal form.

        Returns
        -------
        str
            The waste disposal classification, which can be "Class A", "Class
            B", "Class C", or "GTCC" (greater than class C).

        """
        return waste._waste_classification(self, metal=metal)

    def waste_disposal_rating(
        self,
        limits: str | dict[str, float] = 'Fetter',
        metal: bool = False,
        by_nuclide: bool = False,
    ) -> float | dict[str, float]:
        """Return the waste disposal rating for the material.

        This method returns a waste disposal rating for the material based on a
        set of specific activity limits. The waste disposal rating is a single
        number that represents the sum of the ratios of the specific activity
        for each radionuclide in the material against a nuclide-specific limit.
        A value less than 1.0 indicates that the material "meets" the limits
        whereas a value greater than 1.0 exceeds the limits.

        Note that the limits for NRC do not consider many long-lived
        radionuclides relevant to fusion systems. A paper by `Fetter et al.
        <https://doi.org/10.1016/0920-3796(90)90104-E>`_ applies the NRC
        methodology to calculate specific activity limits for an expanded set of
        radionuclides.

        Parameters
        ----------
        limits : str or dict, optional
            The name of a predefined set of specific activity limits or a
            dictionary that contains specific activity limits for radionuclides,
            where keys are nuclide names and values are activities in units of
            [Ci/m3]. The predefined options are:

            - 'Fetter': Uses limits from Fetter et al. (1990)
            - 'NRC_long': Uses the 10 CFR 61.55 limits for long-lived
              radionuclides
            - 'NRC_short_A': Uses the 10 CFR 61.55 class A limits for
              short-lived radionuclides
            - 'NRC_short_B': Uses the 10 CFR 61.55 class B limits for
              short-lived radionuclides
            - 'NRC_short_C': Uses the 10 CFR 61.55 class C limits for
              short-lived radionuclides
        metal : bool, optional
            Whether or not the material is in metal form (only applicable for
            NRC based limits)
        by_nuclide : bool, optional
            Whether to return the waste disposal rating for each nuclide in the
            material. If True, a dictionary is returned where the keys are the
            nuclide names and the values are the waste disposal ratings for each
            nuclide. If False, a single float value is returned that represents
            the overall waste disposal rating for the material.

        Returns
        -------
        float or dict
            The waste disposal rating for the material or its constituent
            nuclides.

        See also
        --------
        Material.waste_classification()

        """
        return waste._waste_disposal_rating(self, limits, metal, by_nuclide)

    def clone(self, memo: dict | None = None) -> Material:
        """Create a copy of this material with a new unique ID.

        Parameters
        ----------
        memo : dict or None
            A nested dictionary of previously cloned objects. This parameter
            is used internally and should not be specified by the user.

        Returns
        -------
        clone : openmc.Material
            The clone of this material

        """

        if memo is None:
            memo = {}

        # If no nemoize'd clone exists, instantiate one
        if self not in memo:
            # Temporarily remove paths -- this is done so that when the clone is
            # made, it doesn't create a copy of the paths (which are specific to
            # an instance)
            paths = self._paths
            self._paths = None

            clone = deepcopy(self)
            clone.id = None
            clone._num_instances = None

            # Restore paths on original instance
            self._paths = paths

            # Memoize the clone
            memo[self] = clone

        return memo[self]

    def _get_nuclide_xml(self, nuclide: NuclideTuple) -> ET.Element:
        xml_element = ET.Element("nuclide")
        xml_element.set("name", nuclide.name)

        # Prevent subnormal numbers from being written to XML, which causes an
        # exception on the C++ side when calling std::stod
        val = nuclide.percent
        if abs(val) < _SMALLEST_NORMAL:
            val = 0.0

        if nuclide.percent_type == 'ao':
            xml_element.set("ao", str(val))
        else:
            xml_element.set("wo", str(val))

        return xml_element

    def _get_macroscopic_xml(self, macroscopic: str) -> ET.Element:
        xml_element = ET.Element("macroscopic")
        xml_element.set("name", macroscopic)

        return xml_element

    def _get_nuclides_xml(
            self, nuclides: Iterable[NuclideTuple],
            nuclides_to_ignore: Iterable[str] | None = None)-> list[ET.Element]:
        xml_elements = []

        # Remove any nuclides to ignore from the XML export
        if nuclides_to_ignore:
            nuclides = [nuclide for nuclide in nuclides if nuclide.name not in nuclides_to_ignore]

        xml_elements = [self._get_nuclide_xml(nuclide) for nuclide in nuclides]

        return xml_elements

    def to_xml_element(
            self, nuclides_to_ignore: Iterable[str] | None = None) -> ET.Element:
        """Return XML representation of the material

        Parameters
        ----------
        nuclides_to_ignore : list of str
            Nuclides to ignore when exporting to XML.

        Returns
        -------
        element : lxml.etree._Element
            XML element containing material data

        """

        # Create Material XML element
        element = ET.Element("material")
        element.set("id", str(self._id))

        if len(self._name) > 0:
            element.set("name", str(self._name))

        if self._depletable:
            element.set("depletable", "true")

        if self._volume:
            element.set("volume", str(self._volume))

        if self._ncrystal_cfg:
            if self._sab:
                raise ValueError("NCrystal materials are not compatible with S(a,b).")
            if self._macroscopic is not None:
                raise ValueError("NCrystal materials are not compatible with macroscopic cross sections.")

            element.set("cfg", str(self._ncrystal_cfg))

        # Create temperature XML subelement
        if self.temperature is not None:
            element.set("temperature", str(self.temperature))

        # Create density XML subelement
        if self._density is not None or self._density_units == 'sum':
            subelement = ET.SubElement(element, "density")
            if self._density_units != 'sum':
                subelement.set("value", str(self._density))
            subelement.set("units", self._density_units)
        else:
            raise ValueError(f'Density has not been set for material {self.id}!')

        if self._macroscopic is None:
            # Create nuclide XML subelements
            subelements = self._get_nuclides_xml(self._nuclides,
                                                 nuclides_to_ignore=nuclides_to_ignore)
            for subelement in subelements:
                element.append(subelement)
        else:
            # Create macroscopic XML subelements
            subelement = self._get_macroscopic_xml(self._macroscopic)
            element.append(subelement)

        if self._sab:
            for sab in self._sab:
                subelement = ET.SubElement(element, "sab")
                subelement.set("name", sab[0])
                if sab[1] != 1.0:
                    subelement.set("fraction", str(sab[1]))

        if self._isotropic:
            subelement = ET.SubElement(element, "isotropic")
            subelement.text = ' '.join(self._isotropic)

        return element

    @classmethod
    def mix_materials(cls, materials, fracs: Iterable[float],
                      percent_type: str = 'ao', **kwargs) -> Material:
        """Mix materials together based on atom, weight, or volume fractions

        .. versionadded:: 0.12

        Parameters
        ----------
        materials : Iterable of openmc.Material
            Materials to combine
        fracs : Iterable of float
            Fractions of each material to be combined
        percent_type : {'ao', 'wo', 'vo'}
            Type of percentage, must be one of 'ao', 'wo', or 'vo', to signify atom
            percent (molar percent), weight percent, or volume percent,
            optional. Defaults to 'ao'
        **kwargs
            Keyword arguments passed to :class:`openmc.Material`

        Returns
        -------
        openmc.Material
            Mixture of the materials

        """

        cv.check_type('materials', materials, Iterable, Material)
        cv.check_type('fracs', fracs, Iterable, Real)
        cv.check_value('percent type', percent_type, {'ao', 'wo', 'vo'})

        fracs = np.asarray(fracs)
        void_frac = 1. - np.sum(fracs)

        # Warn that fractions don't add to 1, set remainder to void, or raise
        # an error if percent_type isn't 'vo'
        if not np.isclose(void_frac, 0.):
            if percent_type in ('ao', 'wo'):
                msg = ('A non-zero void fraction is not acceptable for '
                       'percent_type: {}'.format(percent_type))
                raise ValueError(msg)
            else:
                msg = ('Warning: sum of fractions do not add to 1, void '
                       'fraction set to {}'.format(void_frac))
                warnings.warn(msg)

        # Calculate appropriate weights which are how many cc's of each
        # material are found in 1cc of the composite material
        amms = np.asarray([mat.average_molar_mass for mat in materials])
        mass_dens = np.asarray([mat.get_mass_density() for mat in materials])
        if percent_type == 'ao':
            wgts = fracs * amms / mass_dens
            wgts /= np.sum(wgts)
        elif percent_type == 'wo':
            wgts = fracs / mass_dens
            wgts /= np.sum(wgts)
        elif percent_type == 'vo':
            wgts = fracs

        # If any of the involved materials contain S(a,b) tables raise an error
        sab_names = set(sab[0] for mat in materials for sab in mat._sab)
        if sab_names:
            msg = ('Currently we do not support mixing materials containing '
                   'S(a,b) tables')
            raise NotImplementedError(msg)

        # Add nuclide densities weighted by appropriate fractions
        nuclides_per_cc = defaultdict(float)
        mass_per_cc = defaultdict(float)
        for mat, wgt in zip(materials, wgts):
            for nuc, atoms_per_bcm in mat.get_nuclide_atom_densities().items():
                nuc_per_cc = wgt*1.e24*atoms_per_bcm
                nuclides_per_cc[nuc] += nuc_per_cc
                mass_per_cc[nuc] += nuc_per_cc*openmc.data.atomic_mass(nuc) / \
                                    openmc.data.AVOGADRO

        # Create the new material with the desired name
        if "name" not in kwargs:
            kwargs["name"] = '-'.join([f'{m.name}({f})' for m, f in
                             zip(materials, fracs)])

        new_mat = cls(**kwargs)

        # Compute atom fractions of nuclides and add them to the new material
        tot_nuclides_per_cc = np.sum([dens for dens in nuclides_per_cc.values()])
        for nuc, atom_dens in nuclides_per_cc.items():
            new_mat.add_nuclide(nuc, atom_dens/tot_nuclides_per_cc, 'ao')

        # Compute mass density for the new material and set it
        new_density = np.sum([dens for dens in mass_per_cc.values()])
        new_mat.set_density('g/cm3', new_density)

        # If any of the involved materials is depletable, the new material is
        # depletable
        new_mat.depletable = any(mat.depletable for mat in materials)

        return new_mat

    @classmethod
    def from_xml_element(cls, elem: ET.Element) -> Material:
        """Generate material from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.Material
            Material generated from XML element

        """
        mat_id = int(get_text(elem, 'id'))

        # Add NCrystal material from cfg string
        cfg = get_text(elem, "cfg")
        if cfg is not None:
            return Material.from_ncrystal(cfg, material_id=mat_id)

        mat = cls(mat_id)
        mat.name = get_text(elem, 'name')

        temperature = get_text(elem, "temperature")
        if temperature is not None:
            mat.temperature = float(temperature)

        volume = get_text(elem, "volume")
        if volume is not None:
            mat.volume = float(volume)

        # Get each nuclide
        for nuclide in elem.findall('nuclide'):
            name = get_text(nuclide, "name")
            if 'ao' in nuclide.attrib:
                mat.add_nuclide(name, float(nuclide.attrib['ao']))
            elif 'wo' in nuclide.attrib:
                mat.add_nuclide(name, float(nuclide.attrib['wo']), 'wo')

        # Get depletable attribute
        depletable = get_text(elem, "depletable")
        mat.depletable = depletable in ('true', '1')

        # Get each S(a,b) table
        for sab in elem.findall('sab'):
            fraction = float(get_text(sab, "fraction", 1.0))
            name = get_text(sab, "name")
            mat.add_s_alpha_beta(name, fraction)

        # Get total material density
        density = elem.find('density')
        units = get_text(density, "units")
        if units == 'sum':
            mat.set_density(units)
        else:
            value = float(get_text(density, 'value'))
            mat.set_density(units, value)

        # Check for isotropic scattering nuclides
        isotropic = get_elem_list(elem, "isotropic", str)
        if isotropic is not None:
            mat.isotropic = isotropic

        return mat

    def deplete(
        self,
        multigroup_flux: Sequence[float],
        energy_group_structure: Sequence[float] | str,
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's',
        chain_file: cv.PathLike | "openmc.deplete.Chain" | None = None,
        reactions: Sequence[str] | None = None,
    ) -> list[openmc.Material]:
        """Depletes that material, evolving the nuclide densities

        .. versionadded:: 0.15.3

        Parameters
        ----------
        multigroup_flux: Sequence[float]
            Energy-dependent multigroup flux values, where each sublist corresponds
            to a specific material. Will be normalized so that it sums to 1.
        energy_group_structure : Sequence[float] | str
            Energy group boundaries in [eV] or the name of the group structure.
        timesteps : iterable of float or iterable of tuple
            Array of timesteps. Note that values are not cumulative. The units are
            specified by the `timestep_units` argument when `timesteps` is an
            iterable of float. Alternatively, units can be specified for each step
            by passing an iterable of (value, unit) tuples.
        source_rates : float or iterable of float, optional
            Source rate in [neutron/sec] or neutron flux in [neutron/s-cm^2] for
            each interval in :attr:`timesteps`
        timestep_units : {'s', 'min', 'h', 'd', 'a', 'MWd/kg'}
            Units for values specified in the `timesteps` argument. 's' means
            seconds, 'min' means minutes, 'h' means hours, 'a' means Julian years
            and 'MWd/kg' indicates that the values are given in burnup (MW-d of
            energy deposited per kilogram of initial heavy metal).
        chain_file : PathLike or Chain
            Path to the depletion chain XML file or instance of openmc.deplete.Chain.
            Defaults to ``openmc.config['chain_file']``.
        reactions : list of str, optional
            Reactions to get cross sections for. If not specified, all neutron
            reactions listed in the depletion chain file are used.

        Returns
        -------
        list of openmc.Material, one for each timestep

        """

        materials = openmc.Materials([self])

        depleted_materials_dict = materials.deplete(
            multigroup_fluxes=[multigroup_flux],
            energy_group_structures=[energy_group_structure],
            timesteps=timesteps,
            source_rates=source_rates,
            timestep_units=timestep_units,
            chain_file=chain_file,
            reactions=reactions,
        )

        return depleted_materials_dict[self.id]


    def mean_free_path(self, energy: float) -> float:
        """Calculate the mean free path of neutrons in the material at a given
        energy.

        .. versionadded:: 0.15.3

        Parameters
        ----------
        energy : float
            Neutron energy in eV

        Returns
        -------
        float
            Mean free path in cm

        """
        from openmc.plotter import _calculate_cexs_elem_mat

        energy_grid, cexs = _calculate_cexs_elem_mat(
            this=self,
            types=["total"],
        )
        total_cexs = cexs[0]

        interpolated_cexs = float(np.interp(energy, energy_grid, total_cexs))

        return 1.0 / interpolated_cexs


class Materials(cv.CheckedList):
    """Collection of Materials used for an OpenMC simulation.

    This class corresponds directly to the materials.xml input file. It can be
    thought of as a normal Python list where each member is a :class:`Material`.
    It behaves like a list as the following example demonstrates:

    >>> fuel = openmc.Material()
    >>> clad = openmc.Material()
    >>> water = openmc.Material()
    >>> m = openmc.Materials([fuel])
    >>> m.append(water)
    >>> m += [clad]

    Parameters
    ----------
    materials : Iterable of openmc.Material
        Materials to add to the collection

    Attributes
    ----------
    cross_sections : str or path-like
        Indicates the path to an XML cross section listing file (usually named
        cross_sections.xml). If it is not set, the
        :envvar:`OPENMC_CROSS_SECTIONS` environment variable will be used for
        continuous-energy calculations and :envvar:`OPENMC_MG_CROSS_SECTIONS`
        will be used for multi-group calculations to find the path to the HDF5
        cross section file.

    """

    def __init__(self, materials=None):
        super().__init__(Material, 'materials collection')
        self._cross_sections = None

        if materials is not None:
            self += materials

    @property
    def cross_sections(self) -> Path | None:
        return self._cross_sections

    @cross_sections.setter
    def cross_sections(self, cross_sections):
        if cross_sections is not None:
            self._cross_sections = input_path(cross_sections)

    def append(self, material):
        """Append material to collection

        Parameters
        ----------
        material : openmc.Material
            Material to append

        """
        super().append(material)

    def insert(self, index: int, material):
        """Insert material before index

        Parameters
        ----------
        index : int
            Index in list
        material : openmc.Material
            Material to insert

        """
        super().insert(index, material)

    def make_isotropic_in_lab(self):
        for material in self:
            material.make_isotropic_in_lab()

    def _write_xml(self, file, header=True, level=0, spaces_per_level=2,
                   trailing_indent=True, nuclides_to_ignore=None):
        """Writes XML content of the materials to an open file handle.

        Parameters
        ----------
        file : IOTextWrapper
            Open file handle to write content into.
        header : bool
            Whether or not to write the XML header
        level : int
            Indentation level of materials element
        spaces_per_level : int
            Number of spaces per indentation
        trailing_indentation : bool
            Whether or not to write a trailing indentation for the materials element
        nuclides_to_ignore : list of str
            Nuclides to ignore when exporting to XML.

        """
        indentation = level*spaces_per_level*' '
        # Write the header and the opening tag for the root element.
        if header:
            file.write("<?xml version='1.0' encoding='utf-8'?>\n")
        file.write(indentation+'<materials>\n')

        # Write the <cross_sections> element.
        if self.cross_sections is not None:
            element = ET.Element('cross_sections')
            element.text = str(self.cross_sections)
            clean_indentation(element, level=level+1)
            element.tail = element.tail.strip(' ')
            file.write((level+1)*spaces_per_level*' ')
            file.write(ET.tostring(element, encoding="unicode"))

        # Write the <material> elements.
        for material in sorted(set(self), key=lambda x: x.id):
            element = material.to_xml_element(nuclides_to_ignore=nuclides_to_ignore)
            clean_indentation(element, level=level+1)
            element.tail = element.tail.strip(' ')
            file.write((level+1)*spaces_per_level*' ')
            file.write(ET.tostring(element, encoding="unicode"))

        # Write the closing tag for the root element.
        file.write(indentation+'</materials>\n')

        # Write a trailing indentation for the next element
        # at this level if needed
        if trailing_indent:
            file.write(indentation)

    def export_to_xml(self, path: PathLike = 'materials.xml',
                      nuclides_to_ignore: Iterable[str] | None = None):
        """Export material collection to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'materials.xml'.
        nuclides_to_ignore : list of str
            Nuclides to ignore when exporting to XML.

        """
        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'materials.xml'

        # Write materials to the file one-at-a-time.  This significantly reduces
        # memory demand over allocating a complete ElementTree and writing it in
        # one go.
        with open(str(p), 'w', encoding='utf-8',
                  errors='xmlcharrefreplace') as fh:
            self._write_xml(fh, nuclides_to_ignore=nuclides_to_ignore)

    @classmethod
    def from_xml_element(cls, elem) -> Materials:
        """Generate materials collection from XML file

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.Materials
            Materials collection

        """
        # Generate each material
        materials = cls()
        for material in elem.findall('material'):
            materials.append(Material.from_xml_element(material))

        # Check for cross sections settings
        xs = get_text(elem, "cross_sections")
        if xs is not None:
            materials.cross_sections = xs

        return materials

    @classmethod
    def from_xml(cls, path: PathLike = 'materials.xml') -> Materials:
        """Generate materials collection from XML file

        Parameters
        ----------
        path : str
            Path to materials XML file

        Returns
        -------
        openmc.Materials
            Materials collection

        """
        parser = ET.XMLParser(huge_tree=True)
        tree = ET.parse(path, parser=parser)
        root = tree.getroot()

        return cls.from_xml_element(root)


    def deplete(
        self,
        multigroup_fluxes: Sequence[Sequence[float]],
        energy_group_structures: Sequence[Sequence[float] | str],
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's',
        chain_file: cv.PathLike | "openmc.deplete.Chain" | None = None,
        reactions: Sequence[str] | None = None,
    ) -> Dict[int, list[openmc.Material]]:
        """Depletes that material, evolving the nuclide densities

        .. versionadded:: 0.15.3

        Parameters
        ----------
        multigroup_fluxes: Sequence[Sequence[float]]
            Energy-dependent multigroup flux values, where each sublist corresponds
            to a specific material. Will be normalized so that it sums to 1.
        energy_group_structures': Sequence[Sequence[float] | str]
            Energy group boundaries in [eV] or the name of the group structure.
        timesteps : iterable of float or iterable of tuple
            Array of timesteps. Note that values are not cumulative. The units are
            specified by the `timestep_units` argument when `timesteps` is an
            iterable of float. Alternatively, units can be specified for each step
            by passing an iterable of (value, unit) tuples.
        source_rates : float or iterable of float, optional
            Source rate in [neutron/sec] or neutron flux in [neutron/s-cm^2] for
            each interval in :attr:`timesteps`
        timestep_units : {'s', 'min', 'h', 'd', 'a', 'MWd/kg'}
            Units for values specified in the `timesteps` argument. 's' means
            seconds, 'min' means minutes, 'h' means hours, 'a' means Julian years
            and 'MWd/kg' indicates that the values are given in burnup (MW-d of
            energy deposited per kilogram of initial heavy metal).
        chain_file : PathLike or Chain
            Path to the depletion chain XML file or instance of openmc.deplete.Chain.
            Defaults to ``openmc.config['chain_file']``.
        reactions : list of str, optional
            Reactions to get cross sections for. If not specified, all neutron
            reactions listed in the depletion chain file are used.

        Returns
        -------
        list of openmc.Material, one for each timestep

        """

        import openmc.deplete
        from .deplete.chain import _get_chain

        # setting all materials to be depletable
        for mat in self:
            mat.depletable = True

        chain = _get_chain(chain_file)

        # Create MicroXS objects for all materials
        micros = []
        fluxes = []

        with openmc.lib.TemporarySession():
            for material, flux, energy in zip(
                self, multigroup_fluxes, energy_group_structures
            ):
                temperature = material.temperature or 293.6
                micro_xs = openmc.deplete.MicroXS.from_multigroup_flux(
                    energies=energy,
                    multigroup_flux=flux,
                    chain_file=chain,
                    temperature=temperature,
                    reactions=reactions,
                )
                micros.append(micro_xs)
                fluxes.append(material.volume)

        # Create a single operator for all materials
        operator = openmc.deplete.IndependentOperator(
            materials=self,
            fluxes=fluxes,
            micros=micros,
            normalization_mode="source-rate",
            chain_file=chain,
        )

        integrator = openmc.deplete.PredictorIntegrator(
            operator=operator,
            timesteps=timesteps,
            source_rates=source_rates,
            timestep_units=timestep_units,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Run integrator
            results_path = Path(tmpdir) / "depletion_results.h5"
            integrator.integrate(path=results_path)

            # Load depletion results
            results = openmc.deplete.Results(results_path)

            # For each material, get activated composition at each timestep
            all_depleted_materials = {
                material.id: [
                    result.get_material(str(material.id))
                    for result in results
                ]
                for material in self
            }

        return all_depleted_materials
