"""The CF4 integrator."""

import copy
from collections.abc import Iterable

from .cram import timed_deplete
from ..results import Results


# Functions to form the special matrix for depletion
def _cf4_f1(chain, rates):
    return 1/2 * chain.form_matrix(rates)


def _cf4_f2(chain, rates):
    return -1/2 * chain.form_matrix(rates[0]) + \
                  chain.form_matrix(rates[1])


def _cf4_f3(chain, rates):
    return  1/4  * chain.form_matrix(rates[0]) + \
            1/6  * chain.form_matrix(rates[1]) + \
            1/6  * chain.form_matrix(rates[2]) + \
           -1/12 * chain.form_matrix(rates[3])


def _cf4_f4(chain, rates):
    return -1/12 * chain.form_matrix(rates[0]) + \
            1/6  * chain.form_matrix(rates[1]) + \
            1/6  * chain.form_matrix(rates[2]) + \
            1/4  * chain.form_matrix(rates[3])


def cf4(operator, timesteps, power=None, power_density=None, print_out=True):
    r"""Deplete using the CF4 algorithm.

    Implements the fourth order `commutator-free Lie algorithm
    <https://doi.org/10.1016/S0167-739X(02)00161-9>`_.
    This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        F_1 &= h A(y_0) \\
        y_1 &= \text{expm}(1/2 F_1) y_0 \\
        F_2 &= h A(y_1) \\
        y_2 &= \text{expm}(1/2 F_2) y_0 \\
        F_3 &= h A(y_2) \\
        y_3 &= \text{expm}(-1/2 F_1 + F_3) y_1 \\
        F_4 &= h A(y_3) \\
        y_4 &= \text{expm}( 1/4  F_1 + 1/6 F_2 + 1/6 F_3 - 1/12 F_4)
               \text{expm}(-1/12 F_1 + 1/6 F_2 + 1/6 F_3 + 1/4  F_4) y_0
        \end{aligned}

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

        chain = operator.chain

        for i, (dt, p) in enumerate(zip(timesteps, power)):
            # Get beginning-of-timestep concentrations and reaction rates
            # Avoid doing first transport run if already done in previous
            # calculation
            if i > 0 or operator.prev_res is None:
                x = [copy.deepcopy(vec)]
                op_results = [operator(x[0], p)]

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
                ratio_power = p / power_res
                op_results[0].rates *= ratio_power[0]

            # Step 1: deplete with matrix 1/2*A(y0)
            time_1, x_new = timed_deplete(
                chain, x[0], op_results[0].rates, dt, print_out,
                matrix_func=_cf4_f1)
            x.append(x_new)
            op_results.append(operator(x_new, p))

            # Step 2: deplete with matrix 1/2*A(y1)
            time_2, x_new = timed_deplete(
                chain, x[0], op_results[1].rates, dt, print_out,
                matrix_func=_cf4_f1)
            x.append(x_new)
            op_results.append(operator(x_new, p))

            # Step 3: deplete with matrix -1/2*A(y0)+A(y2)
            rates = list(zip(op_results[0].rates, op_results[2].rates))
            time_3, x_new = timed_deplete(
                chain, x[1], rates, dt, print_out, matrix_func=_cf4_f2)
            x.append(x_new)
            op_results.append(operator(x_new, p))

            # Step 4: deplete with two matrix exponentials
            rates = list(zip(op_results[0].rates, op_results[1].rates,
                             op_results[2].rates, op_results[3].rates))
            time_4, x_end = timed_deplete(
                chain, x[0], rates, dt, print_out, matrix_func=_cf4_f3)
            time_5, x_end = timed_deplete(
                chain, x_end, rates, dt, print_out, matrix_func=_cf4_f4)

            # Create results, write to disk
            Results.save(
                operator, x, op_results, [t, t + dt], p, i_res + i,
                time_1 + time_2 + time_3 + time_4 + time_5)

            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        op_results = [operator(x[0], power[-1])]

        # Create results, write to disk
        Results.save(operator, x, op_results, [t, t], p, i_res + len(timesteps))
