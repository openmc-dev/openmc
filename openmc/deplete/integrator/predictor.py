"""First-order predictor algorithm."""

import copy
from collections.abc import Iterable

from .cram import deplete
from ..results import Results


def predictor(operator, timesteps, power, print_out=True):
    r"""Deplete using a first-order predictor algorithm.

    Implements the first-order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        y' &= A(y, t) y(t)

        A_p &= A(y_n, t_n)

        y_{n+1} &= \text{expm}(A_p h) y_n

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
            op_results = [operator(x[0], p)]

            # Create results, write to disk
            Results.save(operator, x, op_results, [t, t + dt], i)

            # Deplete for full timestep
            x_end = deplete(chain, x[0], op_results[0], dt, print_out)

            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        op_results = [operator(x[0], power[-1])]

        # Create results, write to disk
        Results.save(operator, x, op_results, [t, t], len(timesteps))
