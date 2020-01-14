from collections.abc import Iterable, Callable, MutableMapping
from copy import deepcopy
from numbers import Real, Integral
from warnings import warn
from io import StringIO

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Uniform, Tabular, Legendre
from .angle_distribution import AngleDistribution
from .angle_energy import AngleEnergy
from .correlated import CorrelatedAngleEnergy
from .data import ATOMIC_SYMBOL, K_BOLTZMANN, EV_PER_MEV
from .endf import get_head_record, get_tab1_record, get_list_record, \
    get_tab2_record, get_cont_record
from .energy_distribution import EnergyDistribution, LevelInelastic, \
    DiscretePhoton
from .function import Tabulated1D, Polynomial
from .kalbach_mann import KalbachMann
from .laboratory import LaboratoryAngleEnergy
from .nbody import NBodyPhaseSpace
from .product import Product
from .uncorrelated import UncorrelatedAngleEnergy


REACTION_NAME = {1: '(n,total)', 2: '(n,elastic)', 4: '(n,level)',
                 5: '(n,misc)', 11: '(n,2nd)', 16: '(n,2n)', 17: '(n,3n)',
                 18: '(n,fission)', 19: '(n,f)', 20: '(n,nf)', 21: '(n,2nf)',
                 22: '(n,na)', 23: '(n,n3a)', 24: '(n,2na)', 25: '(n,3na)',
                 27: '(n,absorption)', 28: '(n,np)', 29: '(n,n2a)',
                 30: '(n,2n2a)', 32: '(n,nd)', 33: '(n,nt)', 34: '(n,nHe-3)',
                 35: '(n,nd2a)', 36: '(n,nt2a)', 37: '(n,4n)', 38: '(n,3nf)',
                 41: '(n,2np)', 42: '(n,3np)', 44: '(n,n2p)', 45: '(n,npa)',
                 91: '(n,nc)', 101: '(n,disappear)', 102: '(n,gamma)',
                 103: '(n,p)', 104: '(n,d)', 105: '(n,t)', 106: '(n,3He)',
                 107: '(n,a)', 108: '(n,2a)', 109: '(n,3a)', 111: '(n,2p)',
                 112: '(n,pa)', 113: '(n,t2a)', 114: '(n,d2a)', 115: '(n,pd)',
                 116: '(n,pt)', 117: '(n,da)', 152: '(n,5n)', 153: '(n,6n)',
                 154: '(n,2nt)', 155: '(n,ta)', 156: '(n,4np)', 157: '(n,3nd)',
                 158: '(n,nda)', 159: '(n,2npa)', 160: '(n,7n)', 161: '(n,8n)',
                 162: '(n,5np)', 163: '(n,6np)', 164: '(n,7np)', 165: '(n,4na)',
                 166: '(n,5na)', 167: '(n,6na)', 168: '(n,7na)', 169: '(n,4nd)',
                 170: '(n,5nd)', 171: '(n,6nd)', 172: '(n,3nt)', 173: '(n,4nt)',
                 174: '(n,5nt)', 175: '(n,6nt)', 176: '(n,2n3He)',
                 177: '(n,3n3He)', 178: '(n,4n3He)', 179: '(n,3n2p)',
                 180: '(n,3n3a)', 181: '(n,3npa)', 182: '(n,dt)',
                 183: '(n,npd)', 184: '(n,npt)', 185: '(n,ndt)',
                 186: '(n,np3He)', 187: '(n,nd3He)', 188: '(n,nt3He)',
                 189: '(n,nta)', 190: '(n,2n2p)', 191: '(n,p3He)',
                 192: '(n,d3He)', 193: '(n,3Hea)', 194: '(n,4n2p)',
                 195: '(n,4n2a)', 196: '(n,4npa)', 197: '(n,3p)',
                 198: '(n,n3p)', 199: '(n,3n2pa)', 200: '(n,5n2p)', 203: '(n,Xp)',
                 204: '(n,Xd)', 205: '(n,Xt)', 206: '(n,X3He)', 207: '(n,Xa)',
                 301: 'heating', 444: 'damage-energy',
                 649: '(n,pc)', 699: '(n,dc)', 749: '(n,tc)', 799: '(n,3Hec)',
                 849: '(n,ac)', 891: '(n,2nc)', 901: 'heating-local'}
REACTION_NAME.update({i: '(n,n{})'.format(i - 50) for i in range(50, 91)})
REACTION_NAME.update({i: '(n,p{})'.format(i - 600) for i in range(600, 649)})
REACTION_NAME.update({i: '(n,d{})'.format(i - 650) for i in range(650, 699)})
REACTION_NAME.update({i: '(n,t{})'.format(i - 700) for i in range(700, 749)})
REACTION_NAME.update({i: '(n,3He{})'.format(i - 750) for i in range(750, 799)})
REACTION_NAME.update({i: '(n,a{})'.format(i - 800) for i in range(800, 849)})
REACTION_NAME.update({i: '(n,2n{})'.format(i - 875) for i in range(875, 891)})


def _get_products(ev, mt):
    """Generate products from MF=6 in an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation
        ENDF evaluation to read from
    mt : int
        The MT value of the reaction to get products for

    Returns
    -------
    products : list of openmc.data.Product
        Products of the reaction

    """
    file_obj = StringIO(ev.section[6, mt])

    # Read HEAD record
    items = get_head_record(file_obj)
    reference_frame = {1: 'laboratory', 2: 'center-of-mass',
                       3: 'light-heavy', 4: 'breakup'}[items[3]]
    n_products = items[4]

    products = []
    for i in range(n_products):
        # Get yield for this product
        params, yield_ = get_tab1_record(file_obj)

        za = int(params[0])
        awr = params[1]
        lip = params[2]
        law = params[3]

        if za == 0:
            p = Product('photon')
        elif za == 1:
            p = Product('neutron')
        elif za == 1000:
            p = Product('electron')
        else:
            Z, A = divmod(za, 1000)
            p = Product('{}{}'.format(ATOMIC_SYMBOL[Z], A))

        p.yield_ = yield_

        """
        # Set reference frame
        if reference_frame == 'laboratory':
            p.center_of_mass = False
        elif reference_frame == 'center-of-mass':
            p.center_of_mass = True
        elif reference_frame == 'light-heavy':
            p.center_of_mass = (awr <= 4.0)
        """

        if law == 0:
            # No distribution given
            pass
        if law == 1:
            # Continuum energy-angle distribution

            # Peak ahead to determine type of distribution
            position = file_obj.tell()
            params = get_cont_record(file_obj)
            file_obj.seek(position)

            lang = params[2]
            if lang == 1:
                p.distribution = [CorrelatedAngleEnergy.from_endf(file_obj)]
            elif lang == 2:
                p.distribution = [KalbachMann.from_endf(file_obj)]

        elif law == 2:
            # Discrete two-body scattering
            params, tab2 = get_tab2_record(file_obj)
            ne = params[5]
            energy = np.zeros(ne)
            mu = []
            for i in range(ne):
                items, values = get_list_record(file_obj)
                energy[i] = items[1]
                lang = items[2]
                if lang == 0:
                    mu.append(Legendre(values))
                elif lang == 12:
                    mu.append(Tabular(values[::2], values[1::2]))
                elif lang == 14:
                    mu.append(Tabular(values[::2], values[1::2],
                                      'log-linear'))

            angle_dist = AngleDistribution(energy, mu)
            dist = UncorrelatedAngleEnergy(angle_dist)
            p.distribution = [dist]
            # TODO: Add level-inelastic info?

        elif law == 3:
            # Isotropic discrete emission
            p.distribution = [UncorrelatedAngleEnergy()]
            # TODO: Add level-inelastic info?

        elif law == 4:
            # Discrete two-body recoil
            pass

        elif law == 5:
            # Charged particle elastic scattering
            pass

        elif law == 6:
            # N-body phase-space distribution
            p.distribution = [NBodyPhaseSpace.from_endf(file_obj)]

        elif law == 7:
            # Laboratory energy-angle distribution
            p.distribution = [LaboratoryAngleEnergy.from_endf(file_obj)]

        products.append(p)

    return products


def _get_fission_products_ace(ace):
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
            coefficients = ace.xss[idx+2 : idx+2+NC].copy()
            for i in range(coefficients.size):
                coefficients[i] *= EV_PER_MEV**(-i)
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
            coefficients = ace.xss[idx+2 : idx+2+NC].copy()
            for i in range(coefficients.size):
                coefficients[i] *= EV_PER_MEV**(-i)
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
            coefficients = ace.xss[idx+2 : idx+2+NC].copy()
            for i in range(coefficients.size):
                coefficients[i] *= EV_PER_MEV**(-i)
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

            # Convert units of inverse shakes to inverse seconds
            delayed_neutron.decay_rate = ace.xss[idx] * 1.e8

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


def _get_fission_products_endf(ev):
    """Generate fission products from an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation

    Returns
    -------
    products : list of openmc.data.Product
        Prompt and delayed fission neutrons
    derived_products : list of openmc.data.Product
        "Total" fission neutron

    """
    products = []
    derived_products = []

    if (1, 456) in ev.section:
        prompt_neutron = Product('neutron')
        prompt_neutron.emission_mode = 'prompt'

        # Prompt nu values
        file_obj = StringIO(ev.section[1, 456])
        lnu = get_head_record(file_obj)[3]
        if lnu == 1:
            # Polynomial representation
            items, coefficients = get_list_record(file_obj)
            prompt_neutron.yield_ = Polynomial(coefficients)
        elif lnu == 2:
            # Tabulated representation
            params, prompt_neutron.yield_ = get_tab1_record(file_obj)

        products.append(prompt_neutron)

    if (1, 452) in ev.section:
        total_neutron = Product('neutron')
        total_neutron.emission_mode = 'total'

        # Total nu values
        file_obj = StringIO(ev.section[1, 452])
        lnu = get_head_record(file_obj)[3]
        if lnu == 1:
            # Polynomial representation
            items, coefficients = get_list_record(file_obj)
            total_neutron.yield_ = Polynomial(coefficients)
        elif lnu == 2:
            # Tabulated representation
            params, total_neutron.yield_ = get_tab1_record(file_obj)

        if (1, 456) in ev.section:
            derived_products.append(total_neutron)
        else:
            products.append(total_neutron)

    if (1, 455) in ev.section:
        file_obj = StringIO(ev.section[1, 455])

        # Determine representation of delayed nu data
        items = get_head_record(file_obj)
        ldg = items[2]
        lnu = items[3]

        if ldg == 0:
            # Delayed-group constants energy independent
            items, decay_constants = get_list_record(file_obj)
            for constant in decay_constants:
                delayed_neutron = Product('neutron')
                delayed_neutron.emission_mode = 'delayed'
                delayed_neutron.decay_rate = constant
                products.append(delayed_neutron)
        elif ldg == 1:
            # Delayed-group constants energy dependent
            raise NotImplementedError('Delayed neutron with energy-dependent '
                                      'group constants.')

        # In MF=1, MT=455, the delayed-group abundances are actually not
        # specified if the group constants are energy-independent. In this case,
        # the abundances must be inferred from MF=5, MT=455 where multiple
        # energy distributions are given.
        if lnu == 1:
            # Nu represented as polynomial
            items, coefficients = get_list_record(file_obj)
            yield_ = Polynomial(coefficients)
            for neutron in products[-6:]:
                neutron.yield_ = deepcopy(yield_)
        elif lnu == 2:
            # Nu represented by tabulation
            params, yield_ = get_tab1_record(file_obj)
            for neutron in products[-6:]:
                neutron.yield_ = deepcopy(yield_)

        if (5, 455) in ev.section:
            file_obj = StringIO(ev.section[5, 455])
            items = get_head_record(file_obj)
            nk = items[4]
            if nk != len(decay_constants):
                raise ValueError(
                    'Number of delayed neutron fission spectra ({}) does not '
                    'match number of delayed neutron precursors ({}).'.format(
                        nk, len(decay_constants)))
            for i in range(nk):
                params, applicability = get_tab1_record(file_obj)
                dist = UncorrelatedAngleEnergy()
                dist.energy = EnergyDistribution.from_endf(file_obj, params)

                delayed_neutron = products[1 + i]
                yield_ = delayed_neutron.yield_

                # Here we handle the fact that the delayed neutron yield is the
                # product of the total delayed neutron yield and the
                # "applicability" of the energy distribution law in file 5.
                if isinstance(yield_, Tabulated1D):
                    if np.all(applicability.y == applicability.y[0]):
                        yield_.y *= applicability.y[0]
                    else:
                        # Get union energy grid and ensure energies are within
                        # interpolable range of both functions
                        max_energy = min(yield_.x[-1], applicability.x[-1])
                        energy = np.union1d(yield_.x, applicability.x)
                        energy = energy[energy <= max_energy]

                        # Calculate group yield
                        group_yield = yield_(energy) * applicability(energy)
                        delayed_neutron.yield_ = Tabulated1D(energy, group_yield)
                elif isinstance(yield_, Polynomial):
                    if len(yield_) == 1:
                        delayed_neutron.yield_ = deepcopy(applicability)
                        delayed_neutron.yield_.y *= yield_.coef[0]
                    else:
                        if np.all(applicability.y == applicability.y[0]):
                            yield_.coef[0] *= applicability.y[0]
                        else:
                            raise NotImplementedError(
                                'Total delayed neutron yield and delayed group '
                                'probability are both energy-dependent.')

                delayed_neutron.distribution.append(dist)

    return products, derived_products


def _get_activation_products(ev, rx):
    """Generate activation products from an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation
        The ENDF evaluation
    rx : openmc.data.Reaction
        Reaction which generates activation products

    Returns
    -------
    products : list of openmc.data.Product
        Activation products

    """
    file_obj = StringIO(ev.section[8, rx.mt])

    # Determine total number of states and whether decay chain is given in a
    # decay sublibrary
    items = get_head_record(file_obj)
    n_states = items[4]
    decay_sublib = (items[5] == 1)

    # Determine if file 9/10 are present
    present = {9: False, 10: False}
    for i in range(n_states):
        if decay_sublib:
            items = get_cont_record(file_obj)
        else:
            items, values = get_list_record(file_obj)
        lmf = items[2]
        if lmf == 9:
            present[9] = True
        elif lmf == 10:
            present[10] = True

    products = []

    for mf in (9, 10):
        if not present[mf]:
            continue

        file_obj = StringIO(ev.section[mf, rx.mt])
        items = get_head_record(file_obj)
        n_states = items[4]
        for i in range(n_states):
            # Determine what the product is
            items, xs = get_tab1_record(file_obj)
            Z, A = divmod(items[2], 1000)
            excited_state = items[3]

            # Get GND name for product
            symbol = ATOMIC_SYMBOL[Z]
            if excited_state > 0:
                name = '{}{}_e{}'.format(symbol, A, excited_state)
            else:
                name = '{}{}'.format(symbol, A)

            p = Product(name)
            if mf == 9:
                p.yield_ = xs
            else:
                # Re-interpolate production cross section and neutron cross
                # section to union energy grid
                energy = np.union1d(xs.x, rx.xs['0K'].x)
                prod_xs = xs(energy)
                neutron_xs = rx.xs['0K'](energy)
                idx = np.where(neutron_xs > 0)

                # Calculate yield as ratio
                yield_ = np.zeros_like(energy)
                yield_[idx] = prod_xs[idx] / neutron_xs[idx]
                p.yield_ = Tabulated1D(energy, yield_)

            # Check if product already exists from MF=6 and if it does, just
            # overwrite the existing yield.
            for product in rx.products:
                if name == product.particle:
                    product.yield_ = p.yield_
                    break
            else:
                products.append(p)

    return products


def _get_photon_products_ace(ace, rx):
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

        if neutron_mt != rx.mt:
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
                             ace.jxs[1] + threshold_idx + n_energy]*EV_PER_MEV

            # Get photon production cross section
            photon_prod_xs = ace.xss[idx + 2:idx + 2 + n_energy]
            neutron_xs = list(rx.xs.values())[0](energy)
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


def _get_photon_products_endf(ev, rx):
    """Generate photon products from an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation
        ENDF evaluation to read from
    rx : openmc.data.Reaction
        Reaction that generates photons

    Returns
    -------
    products : list of openmc.Products
        Photons produced from reaction with given MT

    """
    products = []

    if (12, rx.mt) in ev.section:
        file_obj = StringIO(ev.section[12, rx.mt])

        items = get_head_record(file_obj)
        option = items[2]

        if option == 1:
            # Multiplicities given
            n_discrete_photon = items[4]
            if n_discrete_photon > 1:
                items, total_yield = get_tab1_record(file_obj)
            for k in range(n_discrete_photon):
                photon = Product('photon')

                # Get photon yield
                items, photon.yield_ = get_tab1_record(file_obj)

                # Get photon energy distribution
                law = items[3]
                dist = UncorrelatedAngleEnergy()
                if law == 1:
                    # TODO: Get file 15 distribution
                    pass
                elif law == 2:
                    energy = items[1]
                    primary_flag = items[2]
                    dist.energy = DiscretePhoton(primary_flag, energy,
                                                 ev.target['mass'])

                photon.distribution.append(dist)
                products.append(photon)

        elif option == 2:
            # Transition probability arrays given
            ppyield = {}
            ppyield['type'] = 'transition'
            ppyield['transition'] = transition = {}

            # Determine whether simple (LG=1) or complex (LG=2) transitions
            lg = items[3]

            # Get transition data
            items, values = get_list_record(file_obj)
            transition['energy_start'] = items[0]
            transition['energies'] = np.array(values[::lg + 1])
            transition['direct_probability'] = np.array(values[1::lg + 1])
            if lg == 2:
                # Complex case
                transition['conditional_probability'] = np.array(
                    values[2::lg + 1])

    elif (13, rx.mt) in ev.section:
        file_obj = StringIO(ev.section[13, rx.mt])

        # Determine option
        items = get_head_record(file_obj)
        n_discrete_photon = items[4]
        if n_discrete_photon > 1:
            items, total_xs = get_tab1_record(file_obj)
        for k in range(n_discrete_photon):
            photon = Product('photon')
            items, xs = get_tab1_record(file_obj)

            # Re-interpolate photon production cross section and neutron cross
            # section to union energy grid
            energy = np.union1d(xs.x, rx.xs['0K'].x)
            photon_prod_xs = xs(energy)
            neutron_xs = rx.xs['0K'](energy)
            idx = np.where(neutron_xs > 0)

            # Calculate yield as ratio
            yield_ = np.zeros_like(energy)
            yield_[idx] = photon_prod_xs[idx] / neutron_xs[idx]
            photon.yield_ = Tabulated1D(energy, yield_)

            # Get photon energy distribution
            law = items[3]
            dist = UncorrelatedAngleEnergy()
            if law == 1:
                # TODO: Get file 15 distribution
                pass
            elif law == 2:
                energy = items[1]
                primary_flag = items[2]
                dist.energy = DiscretePhoton(primary_flag, energy,
                                             ev.target['mass'])

            photon.distribution.append(dist)
            products.append(photon)

    return products


class Reaction(EqualityMixin):
    """A nuclear reaction

    A Reaction object represents a single reaction channel for a nuclide with
    an associated cross section and, if present, a secondary angle and energy
    distribution.

    Parameters
    ----------
    mt : int
        The ENDF MT number for this reaction.

    Attributes
    ----------
    center_of_mass : bool
        Indicates whether scattering kinematics should be performed in the
        center-of-mass or laboratory reference frame.
        grid above the threshold value in barns.
    redundant : bool
        Indicates whether or not this is a redundant reaction
    mt : int
        The ENDF MT number for this reaction.
    q_value : float
        The Q-value of this reaction in eV.
    xs : dict of str to openmc.data.Function1D
        Microscopic cross section for this reaction as a function of incident
        energy; these cross sections are provided in a dictionary where the key
        is the temperature of the cross section set.
    products : Iterable of openmc.data.Product
        Reaction products
    derived_products : Iterable of openmc.data.Product
        Derived reaction products. Used for 'total' fission neutron data when
        prompt/delayed data also exists.

    """

    def __init__(self, mt):
        self._center_of_mass = True
        self._redundant = False
        self._q_value = 0.
        self._xs = {}
        self._products = []
        self._derived_products = []

        self.mt = mt

    def __repr__(self):
        if self.mt in REACTION_NAME:
            return "<Reaction: MT={} {}>".format(self.mt, REACTION_NAME[self.mt])
        else:
            return "<Reaction: MT={}>".format(self.mt)

    @property
    def center_of_mass(self):
        return self._center_of_mass

    @property
    def redundant(self):
        return self._redundant

    @property
    def q_value(self):
        return self._q_value

    @property
    def products(self):
        return self._products

    @property
    def derived_products(self):
        return self._derived_products

    @property
    def xs(self):
        return self._xs

    @center_of_mass.setter
    def center_of_mass(self, center_of_mass):
        cv.check_type('center of mass', center_of_mass, (bool, np.bool_))
        self._center_of_mass = center_of_mass

    @redundant.setter
    def redundant(self, redundant):
        cv.check_type('redundant', redundant, (bool, np.bool_))
        self._redundant = redundant

    @q_value.setter
    def q_value(self, q_value):
        cv.check_type('Q value', q_value, Real)
        self._q_value = q_value

    @products.setter
    def products(self, products):
        cv.check_type('reaction products', products, Iterable, Product)
        self._products = products

    @derived_products.setter
    def derived_products(self, derived_products):
        cv.check_type('reaction derived products', derived_products,
                      Iterable, Product)
        self._derived_products = derived_products

    @xs.setter
    def xs(self, xs):
        cv.check_type('reaction cross section dictionary', xs, MutableMapping)
        for key, value in xs.items():
            cv.check_type('reaction cross section temperature', key, str)
            cv.check_type('reaction cross section', value, Callable)
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
        group.attrs['center_of_mass'] = 1 if self.center_of_mass else 0
        group.attrs['redundant'] = 1 if self.redundant else 0
        for T in self.xs:
            Tgroup = group.create_group(T)
            if self.xs[T] is not None:
                dset = Tgroup.create_dataset('xs', data=self.xs[T].y)
                threshold_idx = getattr(self.xs[T], '_threshold_idx', 0)
                dset.attrs['threshold_idx'] = threshold_idx
        for i, p in enumerate(self.products):
            pgroup = group.create_group('product_{}'.format(i))
            p.to_hdf5(pgroup)

    @classmethod
    def from_hdf5(cls, group, energy):
        """Generate reaction from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from
        energy : dict
            Dictionary whose keys are temperatures (e.g., '300K') and values are
            arrays of energies at which cross sections are tabulated at.

        Returns
        -------
        openmc.data.Reaction
            Reaction data

        """

        mt = group.attrs['mt']
        rx = cls(mt)
        rx.q_value = group.attrs['Q_value']
        rx.center_of_mass = bool(group.attrs['center_of_mass'])
        rx.redundant = bool(group.attrs.get('redundant', False))

        # Read cross section at each temperature
        for T, Tgroup in group.items():
            if T.endswith('K'):
                if 'xs' in Tgroup:
                    # Make sure temperature has associated energy grid
                    if T not in energy:
                        raise ValueError(
                            'Could not create reaction cross section for MT={} '
                            'at T={} because no corresponding energy grid '
                            'exists.'.format(mt, T))
                    xs = Tgroup['xs'][()]
                    threshold_idx = Tgroup['xs'].attrs['threshold_idx']
                    tabulated_xs = Tabulated1D(energy[T][threshold_idx:], xs)
                    tabulated_xs._threshold_idx = threshold_idx
                    rx.xs[T] = tabulated_xs

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
        grid = ace.xss[ace.jxs[1]:ace.jxs[1] + n_grid]*EV_PER_MEV

        # Convert data temperature to a "300.0K" number for indexing
        # temperature data
        strT = str(int(round(ace.temperature*EV_PER_MEV / K_BOLTZMANN))) + "K"

        if i_reaction > 0:
            mt = int(ace.xss[ace.jxs[3] + i_reaction - 1])
            rx = cls(mt)

            # Get Q-value of reaction
            rx.q_value = ace.xss[ace.jxs[4] + i_reaction - 1]*EV_PER_MEV

            # ==================================================================
            # CROSS SECTION

            # Get locator for cross-section data
            loc = int(ace.xss[ace.jxs[6] + i_reaction - 1])

            # Determine starting index on energy grid
            threshold_idx = int(ace.xss[ace.jxs[7] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energy = int(ace.xss[ace.jxs[7] + loc])
            energy = grid[threshold_idx:threshold_idx + n_energy]

            # Read reaction cross section
            xs = ace.xss[ace.jxs[7] + loc + 1:ace.jxs[7] + loc + 1 + n_energy]

            # For damage energy production, convert to eV
            if mt == 444:
                xs *= EV_PER_MEV

            # Fix negatives -- known issue for Y89 in JEFF 3.2
            if np.any(xs < 0.0):
                warn("Negative cross sections found for MT={} in {}. Setting "
                     "to zero.".format(rx.mt, ace.name))
                xs[xs < 0.0] = 0.0

            tabulated_xs = Tabulated1D(energy, xs)
            tabulated_xs._threshold_idx = threshold_idx
            rx.xs[strT] = tabulated_xs

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
                    rx.products, rx.derived_products = _get_fission_products_ace(ace)

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

            tabulated_xs = Tabulated1D(grid, elastic_xs)
            tabulated_xs._threshold_idx = 0
            rx.xs[strT] = tabulated_xs

            # No energy distribution for elastic scattering
            neutron = Product('neutron')
            neutron.distribution.append(UncorrelatedAngleEnergy())
            rx.products.append(neutron)

        # ======================================================================
        # ANGLE DISTRIBUTION (FOR UNCORRELATED)

        if i_reaction < ace.nxs[5] + 1:
            # Check if angular distribution data exist
            loc = int(ace.xss[ace.jxs[8] + i_reaction])
            if loc < 0:
                # Angular distribution is given as part of a product
                # angle-energy distribution
                angle_dist = None
            elif loc == 0:
                # Angular distribution is isotropic
                energy = [0.0, grid[-1]]
                mu = Uniform(-1., 1.)
                angle_dist = AngleDistribution(energy, [mu, mu])
            else:
                angle_dist = AngleDistribution.from_ace(ace, ace.jxs[9], loc)

            # Apply angular distribution to each uncorrelated angle-energy
            # distribution
            if angle_dist is not None:
                for d in neutron.distribution:
                    d.angle = angle_dist

        # ======================================================================
        # PHOTON PRODUCTION

        rx.products += _get_photon_products_ace(ace, rx)

        return rx

    @classmethod
    def from_endf(cls, ev, mt):
        """Generate a reaction from an ENDF evaluation

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF evaluation
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        rx : openmc.data.Reaction
            Reaction data

        """
        rx = Reaction(mt)

        # Integrated cross section
        if (3, mt) in ev.section:
            file_obj = StringIO(ev.section[3, mt])
            get_head_record(file_obj)
            params, rx.xs['0K'] = get_tab1_record(file_obj)
            rx.q_value = params[1]

        # Get fission product yields (nu) as well as delayed neutron energy
        # distributions
        if mt in (18, 19, 20, 21, 38):
            rx.products, rx.derived_products = _get_fission_products_endf(ev)

        if (6, mt) in ev.section:
            # Product angle-energy distribution
            for product in _get_products(ev, mt):
                if mt in (18, 19, 20, 21, 38) and product.particle == 'neutron':
                    rx.products[0].applicability = product.applicability
                    rx.products[0].distribution = product.distribution
                else:
                    rx.products.append(product)

        elif (4, mt) in ev.section or (5, mt) in ev.section:
            # Uncorrelated angle-energy distribution
            neutron = Product('neutron')

            # Note that the energy distribution for MT=455 is read in
            # _get_fission_products_endf rather than here
            if (5, mt) in ev.section:
                file_obj = StringIO(ev.section[5, mt])
                items = get_head_record(file_obj)
                nk = items[4]
                for i in range(nk):
                    params, applicability = get_tab1_record(file_obj)
                    dist = UncorrelatedAngleEnergy()
                    dist.energy = EnergyDistribution.from_endf(file_obj, params)

                    neutron.applicability.append(applicability)
                    neutron.distribution.append(dist)
            elif mt == 2:
                # Elastic scattering -- no energy distribution is given since it
                # can be calulcated analytically
                dist = UncorrelatedAngleEnergy()
                neutron.distribution.append(dist)
            elif mt >= 51 and mt < 91:
                # Level inelastic scattering -- no energy distribution is given
                # since it can be calculated analytically. Here we determine the
                # necessary parameters to create a LevelInelastic object
                dist = UncorrelatedAngleEnergy()

                A = ev.target['mass']
                threshold = (A + 1.)/A*abs(rx.q_value)
                mass_ratio = (A/(A + 1.))**2
                dist.energy = LevelInelastic(threshold, mass_ratio)

                neutron.distribution.append(dist)

            if (4, mt) in ev.section:
                for dist in neutron.distribution:
                    dist.angle = AngleDistribution.from_endf(ev, mt)

            if mt in (18, 19, 20, 21, 38) and (5, mt) in ev.section:
                # For fission reactions,
                rx.products[0].applicability = neutron.applicability
                rx.products[0].distribution = neutron.distribution
            else:
                rx.products.append(neutron)

        if (8, mt) in ev.section:
            rx.products += _get_activation_products(ev, rx)

        if (12, mt) in ev.section or (13, mt) in ev.section:
            rx.products += _get_photon_products_endf(ev, rx)

        return rx
