"""The CE/CM integrator."""

import copy
from collections.abc import Iterable

from .cram import deplete
from .save_results import save_results


def cecm(operator, timesteps, power, print_out=True):
    r"""Deplete using the CE/CM algorithm.

    Implements the second order `CE/CM Predictor-Corrector algorithm
    <https://doi.org/10.13182/NSE14-92>`_.  This algorithm is mathematically
    defined as:

    .. math::
        y' &= A(y, t) y(t)

        A_p &= A(y_n, t_n)

        y_m &= \text{expm}(A_p h/2) y_n

        A_c &= A(y_m, t_n + h/2)

        y_{n+1} &= \text{expm}(A_c h) y_n

    Parameters
    ----------
    operator : openmc.deplete.TransportOperator
        The operator object to simulate on.
    timesteps : iterable of float
        Array of timesteps in units of [s]. Note that values are not cumulative.
    power : float or iterable of float
        Power of the reactor in [W]. A single value indicates that the power is
        constant over all timesteps. An iterable indicates potentially different
        power levels for each timestep. For a 2D problem, the power can be given
        in [W/cm] as long as the "volume" assigned to a depletion material is
        actually an area in [cm^2].
    print_out : bool, optional
        Whether or not to print out time.

    """
    if not isinstance(power, Iterable):
        power = [power]*len(timesteps)

    # Generate initial conditions
    with operator as vec:
        chain = operator.chain
        t = 0.0
        for i, (dt, p) in enumerate(zip(timesteps, power)):
            # Get beginning-of-timestep reaction rates
            x = [copy.deepcopy(vec)]
            results = [operator(x[0], p)]

            # Deplete for first half of timestep
            x_middle = deplete(chain, x[0], results[0], dt/2, print_out)

            # Get middle-of-timestep reaction rates
            x.append(x_middle)
            results.append(operator(x_middle, p))

            # Deplete for full timestep using beginning-of-step materials
            x_end = deplete(chain, x[0], results[1], dt, print_out)

            # Create results, write to disk
            save_results(operator, x, results, [t, t + dt], i)

            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        results = [operator(x[0], power[-1])]

        # Create results, write to disk
        save_results(operator, x, results, [t, t], len(timesteps))
