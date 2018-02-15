"""The Predictor algorithm."""

import copy

from .cram import deplete
from .save_results import save_results


def predictor(operator, print_out=True):
    r"""The basic predictor integrator.

    Implements the first-order predictor algorithm. This algorithm is
    mathematically defined as:

    .. math::
        y' &= A(y, t) y(t)

        A_p &= A(y_n, t_n)

        y_{n+1} &= \text{expm}(A_p h) y_n

    Parameters
    ----------
    operator : openmc.deplete.Operator
        The operator object to simulate on.
    print_out : bool, optional
        Whether or not to print out time.

    """
    # Generate initial conditions
    with operator as vec:
        chain = operator.chain
        t = 0.0
        for i, dt in enumerate(operator.settings.dt_vec):
            # Get beginning-of-timestep reaction rates
            x = [copy.deepcopy(vec)]
            results = [operator(x[0])]

            # Create results, write to disk
            save_results(operator, x, results, [t, t + dt], i)

            # Deplete for full timestep
            x_end = deplete(chain, x[0], results[0], dt, print_out)

            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        results = [operator(x[0])]

        # Create results, write to disk
        save_results(operator, x, results, [t, t], len(operator.settings.dt_vec))
