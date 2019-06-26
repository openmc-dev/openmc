"""The CE/LI CFQ4 integrator."""

import copy
from collections.abc import Iterable

from .cram import timed_deplete
from ..results import Results


# Functions to form the special matrix for depletion
def _celi_f1(chain, rates):
    return 5/12 * chain.form_matrix(rates[0]) + \
           1/12 * chain.form_matrix(rates[1])


def _celi_f2(chain, rates):
    return 1/12 * chain.form_matrix(rates[0]) + \
           5/12 * chain.form_matrix(rates[1])


def celi(operator, timesteps, power=None, power_density=None,
         print_out=True):
    r"""Deplete using the CE/LI CFQ4 algorithm.

    Implements the CE/LI Predictor-Corrector algorithm using the `fourth order
    commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    "CE/LI" stands for constant extrapolation on predictor and linear
    interpolation on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_0 &= A(y_n, t_n) \\
        y_p &= \text{expm}(h A_0) y_n \\
        A_1 &= A(y_p, t_n + h) \\
        y_{n+1} &= \text{expm}(\frac{h}{12} A_0 + \frac{5h}{12} A1)
                   \text{expm}(\frac{5h}{12} A_0 + \frac{h}{12} A1) y_n
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

        for i, (dt, p) in enumerate(zip(timesteps, power)):
            vec, t, _ = celi_inner(operator, vec, p, i, i_res, t, dt,
                                   print_out)

        # Perform one last simulation
        x = [copy.deepcopy(vec)]
        op_results = [operator(x[0], power[-1])]

        # Create results, write to disk
        Results.save(operator, x, op_results, [t, t], p, i_res + len(timesteps))


def celi_inner(operator, vec, p, i, i_res, t, dt, print_out):
    """ The inner loop of CE/LI CFQ4.

    Parameters
    ----------
    operator : Operator
        The operator object to simulate on.
    x : list of nuclide vector
        Nuclide vector, beginning of time.
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

    Returns
    -------
    list of numpy.array
        Nuclide vector, end of time.
    float
        Next time
    OperatorResult
        Operator result from beginning of step.
    """

    chain = operator.chain

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

    # Deplete to end
    proc_time, x_new = timed_deplete(chain, x[0], op_results[0].rates, dt, print_out)
    x.append(x_new)
    op_results.append(operator(x[1], p))

    # Deplete with two matrix exponentials
    rates = list(zip(op_results[0].rates, op_results[1].rates))
    time_1, x_end = timed_deplete(chain, x[0], rates, dt, print_out,
                    matrix_func=_celi_f1)
    time_2, x_end = timed_deplete(chain, x_end, rates, dt, print_out,
                    matrix_func=_celi_f2)

    # Create results, write to disk
    Results.save(operator, x, op_results, [t, t + dt], p, i_res + i, proc_time + time_1 + time_2)

    # return updated time and vectors
    return x_end, t + dt, op_results[0]
