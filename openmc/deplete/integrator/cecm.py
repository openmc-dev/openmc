"""The CE/CM integrator."""

import copy

from .cram import deplete
from .save_results import save_results


def cecm(operator, print_out=True):
    r"""The CE/CM integrator.

    Implements the second order CE/CM Predictor-Corrector algorithm [ref]_.
    This algorithm is mathematically defined as:

    .. math::
        y' &= A(y, t) y(t)

        A_p &= A(y_n, t_n)

        y_m &= \text{expm}(A_p h/2) y_n

        A_c &= A(y_m, t_n + h/2)

        y_{n+1} &= \text{expm}(A_c h) y_n

    .. [ref]
        Isotalo, Aarno. "Comparison of Neutronics-Depletion Coupling Schemes
        for Burnup Calculations-Continued Study." Nuclear Science and
        Engineering 180.3 (2015): 286-300.

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

            # Deplete for first half of timestep
            x_middle = deplete(chain, x[0], results[0], dt/2, print_out)

            # Get middle-of-timestep reaction rates
            x.append(x_middle)
            results.append(operator(x_middle))

            # Deplete for full timestep using beginning-of-step materials
            x_end = deplete(chain, x[0], results[1], dt, print_out)

            # Create results, write to disk
            save_results(operator, x, results, [t, t + dt], i)

            # Advance time, update vector
            t += dt
            vec = copy.deepcopy(x_end)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        results = [operator(x[0])]

        # Create results, write to disk
        save_results(operator, x, results, [t, t], len(operator.settings.dt_vec))
