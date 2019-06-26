"""The SI-LE/QI CFQ4 integrator."""

import copy
from collections.abc import Iterable
from itertools import repeat

from .si_celi import si_celi_inner
from .leqi import _leqi_f1, _leqi_f2, _leqi_f3, _leqi_f4
from .cram import timed_deplete
from ..results import Results
from ..abc import OperatorResult


def si_leqi(operator, timesteps, power=None, power_density=None,
            print_out=True, m=10):
    r"""Deplete using the SI-LE/QI CFQ4 algorithm.

    Implements the Stochastic Implicit LE/QI Predictor-Corrector algorithm using
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

        chain = operator.chain

        for i, (dt, p) in enumerate(zip(timesteps, power)):
            # LE/QI needs the last step results to start
            # Perform SI-CE/LI CFQ4 or restore results for the first step
            if i == 0:
                dt_l = dt
                if i_res <= 1:
                    op_res_last = copy.deepcopy(op_results[0])
                    x, t, op_results = si_celi_inner(operator, x, op_results, p,
                                                     i, i_res, t, dt, print_out)
                    continue
                else:
                    dt_l = t - operator.prev_res[-2].time[0]
                    op_res_last = operator.prev_res[-2]
                    op_res_last.rates = op_res_last.rates[0]
                    x = [operator.prev_res[-1].data[0]]

            # Perform remaining LE/QI
            inputs = list(zip(op_res_last.rates, op_results[0].rates,
                              repeat(dt_l), repeat(dt)))
            proc_time, x_new = timed_deplete(
                chain, x[0], inputs, dt, print_out, matrix_func=_leqi_f1)
            time_1, x_new = timed_deplete(
                chain, x_new, inputs, dt, print_out, matrix_func=_leqi_f2)
            x.append(x_new)

            proc_time += time_1

            # Loop on inner
            for j in range(m + 1):
                op_res = operator(x_new, p)

                if j <= 1:
                    op_res_bar = copy.deepcopy(op_res)
                else:
                    rates = 1/j * op_res.rates + (1 - 1/j) * op_res_bar.rates
                    k = 1/j * op_res.k + (1 - 1/j) * op_res_bar.k
                    op_res_bar = OperatorResult(k, rates)

                inputs = list(zip(op_res_last.rates, op_results[0].rates,
                                  op_res_bar.rates, repeat(dt_l), repeat(dt)))
                time_1, x_new = timed_deplete(
                    chain, x[0], inputs, dt, print_out, matrix_func=_leqi_f3)
                time_2, x_new = timed_deplete(
                    chain, x_new, inputs, dt, print_out, matrix_func=_leqi_f4)

                proc_time += time_1 + time_2

            # Create results, write to disk
            op_results.append(op_res_bar)
            Results.save(
                operator, x, op_results, [t, t+dt], p, i_res+i, proc_time)

            # update results
            x = [x_new]
            op_res_last = copy.deepcopy(op_results[0])
            op_results = [op_res_bar]
            t += dt
            dt_l = dt

        # Create results for last point, write to disk
        Results.save(
                operator, x, op_results, [t, t], p, i_res+len(timesteps))
