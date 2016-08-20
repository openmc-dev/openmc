from __future__ import division, unicode_literals
from collections import Iterable, Callable
from copy import deepcopy
from numbers import Real
from warnings import warn

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Uniform
from .angle_distribution import AngleDistribution
from .angle_energy import AngleEnergy
from .function import Tabulated1D, Polynomial
from .data import REACTION_NAME
from .product import Product
from .uncorrelated import UncorrelatedAngleEnergy


def _get_fission_products(ace):
    """Generate fission products from an ACE table

    Parameters
    ----------
    ace : openmc.data.ace.Table
        ACE table to read from

    Returns
    -------
    products : list of openmc.data.Product
        Prompt and delayed fission neutrons
    derived_products : list of openmc.data.Product
        "Total" fission neutron

    """
    # No NU block
    if ace.jxs[2] == 0:
        return None, None

    products = []
    derived_products = []

    # Either prompt nu or total nu is given
    if ace.xss[ace.jxs[2]] > 0:
        whichnu = 'prompt' if ace.jxs[24] > 0 else 'total'

        neutron = Product('neutron')
        neutron.emission_mode = whichnu

        idx = ace.jxs[2]
        LNU = int(ace.xss[idx])
        if LNU == 1:
            # Polynomial function form of nu
            NC = int(ace.xss[idx+1])
            coefficients = ace.xss[idx+2 : idx+2+NC]
            neutron.yield_ = Polynomial(coefficients)
        elif LNU == 2:
            # Tabular data form of nu
            neutron.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        products.append(neutron)

    # Both prompt nu and total nu
    elif ace.xss[ace.jxs[2]] < 0:
        # Read prompt neutron yield
        prompt_neutron = Product('neutron')
        prompt_neutron.emission_mode = 'prompt'

        idx = ace.jxs[2] + 1
        LNU = int(ace.xss[idx])
        if LNU == 1:
            # Polynomial function form of nu
            NC = int(ace.xss[idx+1])
            coefficients = ace.xss[idx+2 : idx+2+NC]
            prompt_neutron.yield_ = Polynomial(coefficients)
        elif LNU == 2:
            # Tabular data form of nu
            prompt_neutron.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        # Read total neutron yield
        total_neutron = Product('neutron')
        total_neutron.emission_mode = 'total'

        idx = ace.jxs[2] + int(abs(ace.xss[ace.jxs[2]])) + 1
        LNU = int(ace.xss[idx])

        if LNU == 1:
            # Polynomial function form of nu
            NC = int(ace.xss[idx+1])
            coefficients = ace.xss[idx+2 : idx+2+NC]
            total_neutron.yield_ = Polynomial(coefficients)
        elif LNU == 2:
            # Tabular data form of nu
            total_neutron.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        products.append(prompt_neutron)
        derived_products.append(total_neutron)

    # Check for delayed nu data
    if ace.jxs[24] > 0:
        yield_delayed = Tabulated1D.from_ace(ace, ace.jxs[24] + 1)

        # Delayed neutron precursor distribution
        idx = ace.jxs[25]
        n_group = ace.nxs[8]
        total_group_probability = 0.
        for group in range(n_group):
            delayed_neutron = Product('neutron')
            delayed_neutron.emission_mode = 'delayed'
            delayed_neutron.decay_rate = ace.xss[idx]

            group_probability = Tabulated1D.from_ace(ace, idx + 1)
            if np.all(group_probability.y == group_probability.y[0]):
                delayed_neutron.yield_ = deepcopy(yield_delayed)
                delayed_neutron.yield_.y *= group_probability.y[0]
                total_group_probability += group_probability.y[0]
            else:
                # Get union energy grid and ensure energies are within
                # interpolable range of both functions
                max_energy = min(yield_delayed.x[-1], group_probability.x[-1])
                energy = np.union1d(yield_delayed.x, group_probability.x)
                energy = energy[energy <= max_energy]

                # Calculate group yield
                group_yield = yield_delayed(energy) * group_probability(energy)
                delayed_neutron.yield_ = Tabulated1D(energy, group_yield)

            # Advance position
            nr = int(ace.xss[idx + 1])
            ne = int(ace.xss[idx + 2 + 2*nr])
            idx += 3 + 2*nr + 2*ne

            # Energy distribution for delayed fission neutrons
            location_start = int(ace.xss[ace.jxs[26] + group])
            delayed_neutron.distribution.append(
                AngleEnergy.from_ace(ace, ace.jxs[27], location_start))

            products.append(delayed_neutron)

        # Renormalize delayed neutron yields to reflect fact that in ACE
        # file, the sum of the group probabilities is not exactly one
        for product in products[1:]:
            if total_group_probability > 0.:
                product.yield_.y /= total_group_probability

    return products, derived_products


def _get_photon_products(ace, rx):
    """Generate photon products from an ACE table

    Parameters
    ----------
    ace : openmc.data.ace.Table
        ACE table to read from
    rx : openmc.data.Reaction
        Reaction that generates photons

    Returns
    -------
    photons : list of openmc.Products
        Photons produced from reaction with given MT

    """
    n_photon_reactions = ace.nxs[6]
    photon_mts = ace.xss[ace.jxs[13]:ace.jxs[13] +
                         n_photon_reactions].astype(int)

    photons = []
    for i in range(n_photon_reactions):
        # Determine corresponding reaction
        neutron_mt = photon_mts[i] // 1000

        # Restrict to photons that match the requested MT. Note that if the
        # photon is assigned to MT=18 but the file splits fission into
        # MT=19,20,21,38, we assign the photon product to each of the individual
        # reactions
        if neutron_mt == 18:
            if rx.mt not in (18, 19, 20, 21, 38):
                continue
        elif neutron_mt != rx.mt:
            continue

        # Create photon product and assign to reactions
        photon = Product('photon')

        # ==================================================================
        # Photon yield / production cross section

        loca = int(ace.xss[ace.jxs[14] + i])
        idx = ace.jxs[15] + loca - 1
        mftype = int(ace.xss[idx])
        idx += 1

        if mftype in (12, 16):
            # Yield data taken from ENDF File 12 or 6
            mtmult = int(ace.xss[idx])
            assert mtmult == neutron_mt

            # Read photon yield as function of energy
            photon.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        elif mftype == 13:
            # Cross section data from ENDF File 13

            # Energy grid index at which data starts
            threshold_idx = int(ace.xss[idx]) - 1
            n_energy = int(ace.xss[idx + 1])
            energy = ace.xss[ace.jxs[1] + threshold_idx:
                             ace.jxs[1] + threshold_idx + n_energy]

            # Get photon production cross section
            photon_prod_xs = ace.xss[idx + 2:idx + 2 + n_energy]
            neutron_xs = rx.xs(energy)
            idx = np.where(neutron_xs > 0.)

            # Calculate photon yield
            yield_ = np.zeros_like(photon_prod_xs)
            yield_[idx] = photon_prod_xs[idx] / neutron_xs[idx]
            photon.yield_ = Tabulated1D(energy, yield_)

        else:
            raise ValueError("MFTYPE must be 12, 13, 16. Got {0}".format(
                mftype))

        # ==================================================================
        # Photon energy distribution

        location_start = int(ace.xss[ace.jxs[18] + i])
        distribution = AngleEnergy.from_ace(ace, ace.jxs[19], location_start)
        assert isinstance(distribution, UncorrelatedAngleEnergy)

        # ==================================================================
        # Photon angular distribution
        loc = int(ace.xss[ace.jxs[16] + i])

        if loc == 0:
            # No angular distribution data are given for this reaction,
            # isotropic scattering is asssumed in LAB
            energy = np.array([photon.yield_.x[0], photon.yield_.x[-1]])
            mu_isotropic = Uniform(-1., 1.)
            distribution.angle = AngleDistribution(
                energy, [mu_isotropic, mu_isotropic])
        else:
            distribution.angle = AngleDistribution.from_ace(ace, ace.jxs[17], loc)

        # Add to list of distributions
        photon.distribution.append(distribution)
        photons.append(photon)

    return photons


class Reaction(EqualityMixin):
    """A nuclear reaction

    A Reaction object represents a single reaction channel for a nuclide with
    an associated cross section and, if present, a secondary angle and energy
    distribution.

    Parameters
    ----------
    mt : int
        The ENDF MT number for this reaction. On occasion, MCNP uses MT numbers
        that don't correspond exactly to the ENDF specification.

    Attributes
    ----------
    center_of_mass : bool
        Indicates whether scattering kinematics should be performed in the
        center-of-mass or laboratory reference frame.
        grid above the threshold value in barns.
    mt : int
        The ENDF MT number for this reaction.
    q_value : float
        The Q-value of this reaction in MeV.
    table : openmc.data.ace.Table
        The ACE table which contains this reaction.
    threshold : float
        Threshold of the reaction in MeV
    threshold_idx : int
        The index on the energy grid corresponding to the threshold of this
        reaction.
    xs : callable
        Microscopic cross section for this reaction as a function of incident
        energy
    products : Iterable of openmc.data.Product
        Reaction products
    derived_products : Iterable of openmc.data.Product
        Derived reaction products. Used for 'total' fission neutron data when
        prompt/delayed data also exists.

    """

    def __init__(self, mt):
        self.center_of_mass = True
        self.mt = mt
        self.q_value = 0.
        self.threshold_idx = 0
        self._xs = None
        self.products = []
        self.derived_products = []

    def __repr__(self):
        if self.mt in REACTION_NAME:
            return "<Reaction: MT={} {}>".format(self.mt, REACTION_NAME[self.mt])
        else:
            return "<Reaction: MT={}>".format(self.mt)

    @property
    def center_of_mass(self):
        return self._center_of_mass

    @property
    def q_value(self):
        return self._q_value

    @property
    def products(self):
        return self._products

    @property
    def threshold(self):
        return self.xs.x[0]

    @property
    def xs(self):
        return self._xs

    @center_of_mass.setter
    def center_of_mass(self, center_of_mass):
        cv.check_type('center of mass', center_of_mass, (bool, np.bool_))
        self._center_of_mass = center_of_mass

    @q_value.setter
    def q_value(self, q_value):
        cv.check_type('Q value', q_value, Real)
        self._q_value = q_value

    @products.setter
    def products(self, products):
        cv.check_type('reaction products', products, Iterable, Product)
        self._products = products

    @xs.setter
    def xs(self, xs):
        cv.check_type('reaction cross section', xs, Callable)
        if isinstance(xs, Tabulated1D):
            for y in xs.y:
                cv.check_greater_than('reaction cross section', y, 0.0, True)
        self._xs = xs

    def to_hdf5(self, group):
        """Write reaction to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['mt'] = self.mt
        if self.mt in REACTION_NAME:
            group.attrs['label'] = np.string_(REACTION_NAME[self.mt])
        else:
            group.attrs['label'] = np.string_(self.mt)
        group.attrs['Q_value'] = self.q_value
        group.attrs['threshold_idx'] = self.threshold_idx + 1
        group.attrs['center_of_mass'] = 1 if self.center_of_mass else 0
        if self.xs is not None:
            group.create_dataset('xs', data=self.xs.y)
        for i, p in enumerate(self.products):
            pgroup = group.create_group('product_{}'.format(i))
            p.to_hdf5(pgroup)

    @classmethod
    def from_hdf5(cls, group, energy):
        """Generate reaction from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        energy : Iterable of float
            Array of energies at which cross sections are tabulated at

        Returns
        -------
        openmc.data.ace.Reaction
            Reaction data

        """
        mt = group.attrs['mt']
        rx = cls(mt)
        rx.q_value = group.attrs['Q_value']
        rx.threshold_idx = group.attrs['threshold_idx'] - 1
        rx.center_of_mass = bool(group.attrs['center_of_mass'])

        # Read cross section
        if 'xs' in group:
            xs = group['xs'].value
            rx.xs = Tabulated1D(energy[rx.threshold_idx:], xs)

        # Determine number of products
        n_product = 0
        for name in group:
            if name.startswith('product_'):
                n_product += 1

        # Read reaction products
        for i in range(n_product):
            pgroup = group['product_{}'.format(i)]
            rx.products.append(Product.from_hdf5(pgroup))

        return rx

    @classmethod
    def from_ace(cls, ace, i_reaction):
        # Get nuclide energy grid
        n_grid = ace.nxs[3]
        grid = ace.xss[ace.jxs[1]:ace.jxs[1] + n_grid]

        if i_reaction > 0:
            mt = int(ace.xss[ace.jxs[3] + i_reaction - 1])
            rx = cls(mt)

            # Get Q-value of reaction
            rx.q_value = ace.xss[ace.jxs[4] + i_reaction - 1]

            # ==================================================================
            # CROSS SECTION

            # Get locator for cross-section data
            loc = int(ace.xss[ace.jxs[6] + i_reaction - 1])

            # Determine starting index on energy grid
            rx.threshold_idx = int(ace.xss[ace.jxs[7] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energy = int(ace.xss[ace.jxs[7] + loc])
            energy = grid[rx.threshold_idx:rx.threshold_idx + n_energy]

            # Read reaction cross section
            xs = ace.xss[ace.jxs[7] + loc + 1:ace.jxs[7] + loc + 1 + n_energy]

            # Fix negatives -- known issue for Y89 in JEFF 3.2
            if np.any(xs < 0.0):
                warn("Negative cross sections found for MT={} in {}. Setting "
                     "to zero.".format(rx.mt, ace.name))
                xs[xs < 0.0] = 0.0

            rx.xs = Tabulated1D(energy, xs)

            # ==================================================================
            # YIELD AND ANGLE-ENERGY DISTRIBUTION

            # Determine multiplicity
            ty = int(ace.xss[ace.jxs[5] + i_reaction - 1])
            rx.center_of_mass = (ty < 0)
            if i_reaction < ace.nxs[5] + 1:
                if ty != 19:
                    if abs(ty) > 100:
                        # Energy-dependent neutron yield
                        idx = ace.jxs[11] + abs(ty) - 101
                        yield_ = Tabulated1D.from_ace(ace, idx)
                    else:
                        # 0-order polynomial i.e. a constant
                        yield_ = Polynomial((abs(ty),))

                    neutron = Product('neutron')
                    neutron.yield_ = yield_
                    rx.products.append(neutron)
                else:
                    assert mt in (18, 19, 20, 21, 38)
                    rx.products, rx.derived_products = _get_fission_products(ace)

                    for p in rx.products:
                        if p.emission_mode in ('prompt', 'total'):
                            neutron = p
                            break
                    else:
                        raise Exception("Couldn't find prompt/total fission neutron")

                # Determine locator for ith energy distribution
                lnw = int(ace.xss[ace.jxs[10] + i_reaction - 1])
                while lnw > 0:
                    # Applicability of this distribution
                    neutron.applicability.append(Tabulated1D.from_ace(
                        ace, ace.jxs[11] + lnw + 2))

                    # Read energy distribution data
                    neutron.distribution.append(AngleEnergy.from_ace(
                        ace, ace.jxs[11], lnw, rx))

                    lnw = int(ace.xss[ace.jxs[11] + lnw - 1])

        else:
            # Elastic scattering
            mt = 2
            rx = cls(mt)

            # Get elastic cross section values
            elastic_xs = ace.xss[ace.jxs[1] + 3*n_grid:ace.jxs[1] + 4*n_grid]

            # Fix negatives -- known issue for Ti46,49,50 in JEFF 3.2
            if np.any(elastic_xs < 0.0):
                warn("Negative elastic scattering cross section found for {}. "
                     "Setting to zero.".format(ace.name))
                elastic_xs[elastic_xs < 0.0] = 0.0

            rx.xs = Tabulated1D(grid, elastic_xs)

            # No energy distribution for elastic scattering
            neutron = Product('neutron')
            neutron.distribution.append(UncorrelatedAngleEnergy())
            rx.products.append(neutron)

        # ======================================================================
        # ANGLE DISTRIBUTION (FOR UNCORRELATED)

        if i_reaction < ace.nxs[5] + 1:
            # Check if angular distribution data exist
            loc = int(ace.xss[ace.jxs[8] + i_reaction])
            if loc <= 0:
                # Angular distribution is either given as part of a product
                # angle-energy distribution or is not given at all (in which
                # case isotropic scattering is assumed)
                angle_dist = None
            else:
                angle_dist = AngleDistribution.from_ace(ace, ace.jxs[9], loc)

            # Apply angular distribution to each uncorrelated angle-energy
            # distribution
            if angle_dist is not None:
                for d in neutron.distribution:
                    d.angle = angle_dist

        # ======================================================================
        # PHOTON PRODUCTION

        rx.products += _get_photon_products(ace, rx)

        return rx
