"""The SI-CE/LI CFQ4 integrator."""

import copy
from collections.abc import Iterable

from .cram import timed_deplete
from ..results import Results
from ..abc import OperatorResult
from .celi import _celi_f1, _celi_f2


def si_celi(operator, timesteps, power=None, power_density=None,
            print_out=True, m=10):
    r"""Deplete using the SI-CE/LI CFQ4 algorithm.

    Implements the Stochastic Implicit CE/LI Predictor-Corrector algorithm using
    the `fourth order commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    Detailed algorithm can be found in Section 3.2 in `Colin Josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        The operator object to simulate on.
    timesteps : iterable of float
        Array of timesteps in units of [s]. Note that values are not cumulative.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that the power is
        constant over all timesteps. An iterable indicates potentially different
        power levels for each timestep. For a 2D problem, the power can be given
        in [W/cm] as long as the "volume" assigned to a depletion material is
        actually an area in [cm^2]. Either `power` or `power_density` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by initial
        heavy metal inventory to get total power if `power` is not speficied.
    print_out : bool, optional
        Whether or not to print out time.
    m : int, optional
        Number of stages.
    """
    if power is None:
        if power_density is None:
            raise ValueError(
                "Neither power nor power density was specified.")
        if not isinstance(power_density, Iterable):
            power = power_density*operator.heavy_metal
        else:
            power = [i*operator.heavy_metal for i in power_density]

    if not isinstance(power, Iterable):
        power = [power]*len(timesteps)

    # Generate initial conditions
    with operator as vec:
        # Initialize time and starting index
        if operator.prev_res is None:
            t = 0.0
            i_res = 0
        else:
            t = operator.prev_res[-1].time[-1]
            i_res = len(operator.prev_res)

        # Get the concentrations and reaction rates for the first
        # beginning-of-timestep (BOS). Compute with m (stage number) times as
        # many neutrons as later simulations for statistics reasons if no
        # previous calculation results present
        if operator.prev_res is None:
            x = [copy.deepcopy(vec)]
            if hasattr(operator, "settings"):
                operator.settings.particles *= m
            op_results = [operator(x[0], power[0])]
            if hasattr(operator, "settings"):
                operator.settings.particles //= m
        else:
            # Get initial concentration
            x = [operator.prev_res[-1].data[0]]

            # Get rates
            op_results = [operator.prev_res[-1]]
            op_results[0].rates = op_results[0].rates[0]

            # Set first stage value of keff
            op_results[0].k = op_results[0].k[0]

            # Scale reaction rates by ratio of powers
            power_res = operator.prev_res[-1].power
            ratio_power = power[0] / power_res
            op_results[0].rates *= ratio_power[0]

        for i, (dt, p) in enumerate(zip(timesteps, power)):
            x, t, op_results = si_celi_inner(operator, x, op_results, p,
                                             i, i_res, t, dt, print_out, m)

        # Create results for last point, write to disk
        Results.save(
            operator, x, op_results, [t, t], p, i_res + len(timesteps))


def si_celi_inner(operator, x, op_results, p, i, i_res, t, dt, print_out, m=10):
    """ The inner loop of SI-CE/LI CFQ4.

    Parameters
    ----------
    operator : Operator
        The operator object to simulate on.
    x : list of nuclide vector
        Nuclide vector, beginning of time.
    op_results : list of OperatorResult
        Operator result at BOS.
    p : float
        Power of the reactor in [W]
    i : int
        Current iteration number.
    i_res : int
        Starting index, for restart calculation.
    t : float
        Time at start of step.
    dt : float
        Time step.
    print_out : bool
        Whether or not to print out time.
    m : int, optional
        Number of stages.

    Returns
    -------
    list of nuclide vector (numpy.array)
        Nuclide vector, end of time.
    float
        Next time
    list of OperatorResult
        Operator result at end of time.
    """

    chain = operator.chain

    # Deplete to end
    proc_time, x_new = timed_deplete(
        chain, x[0], op_results[0].rates, dt, print_out)
    x.append(x_new)

    for j in range(m + 1):
        op_res = operator(x_new, p)

        if j <= 1:
            op_res_bar = copy.deepcopy(op_res)
        else:
            rates = 1/j * op_res.rates + (1 - 1/j) * op_res_bar.rates
            k = 1/j * op_res.k + (1 - 1/j) * op_res_bar.k
            op_res_bar = OperatorResult(k, rates)

        rates = list(zip(op_results[0].rates, op_res_bar.rates))
        time_1, x_new = timed_deplete(
            chain, x[0], rates, dt, print_out, matrix_func=_celi_f1)
        time_2, x_new = timed_deplete(
            chain, x_new, rates, dt, print_out, matrix_func=_celi_f2)
        proc_time += time_1 + time_2

    # Create results, write to disk
    op_results.append(op_res_bar)
    Results.save(operator, x, op_results, [t, t+dt], p, i_res+i, proc_time)

    # return updated time and vectors
    return [x_new], t + dt, [op_res_bar]
