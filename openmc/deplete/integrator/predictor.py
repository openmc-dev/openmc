"""First-order predictor algorithm."""

import copy
from collections.abc import Iterable

from .cram import timed_deplete
from ..results import Results


def predictor(operator, timesteps, power=None, power_density=None,
              print_out=True):
    r"""Deplete using a first-order predictor algorithm.

    Implements the first-order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_p &= A(y_n, t_n) \\
        y_{n+1} &= \text{expm}(A_p h) y_n
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

    proc_time = None

    # Generate initial conditions
    with operator as vec:
        # Initialize time and starting index
        if operator.prev_res is None:
            t = 0.0
            i_res = 0
        else:
            t = operator.prev_res[-1].time[-1]
            i_res = len(operator.prev_res) - 1

        chain = operator.chain

        for i, (dt, p) in enumerate(zip(timesteps, power)):
            # Get beginning-of-timestep concentrations and reaction rates
            # Avoid doing first transport run if already done in previous
            # calculation
            if i > 0 or operator.prev_res is None:
                x = [copy.deepcopy(vec)]
                op_results = [operator(x[0], p)]

                # Create results, write to disk
                Results.save(operator, x, op_results, [t, t + dt], p, i_res + i, proc_time)
            else:
                # Get initial concentration
                x = [operator.prev_res[-1].data[0]]

                # Get rates
                op_results = [operator.prev_res[-1]]
                op_results[0].rates = op_results[0].rates[0]

                # Scale reaction rates by ratio of powers
                power_res = operator.prev_res[-1].power
                ratio_power = p / power_res
                op_results[0].rates *= ratio_power[0]

            # Deplete for full timestep
            proc_time, x_end = timed_deplete(
                chain, x[0], op_results[0].rates, dt, print_out)

            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        op_results = [operator(x[0], power[-1])]

        # Create results, write to disk
        Results.save(operator, x, op_results, [t, t], p, i_res + len(timesteps), proc_time)
