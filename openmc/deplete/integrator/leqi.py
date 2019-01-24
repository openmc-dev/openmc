"""The LE/QI CFQ4 integrator."""

import copy
from collections.abc import Iterable
from itertools import repeat

from .celi import celi_inner
from .cram import deplete
from ..results import Results
from ..abc import OperatorResult


# Functions to form the special matrix for depletion
def _leqi_f1(chain, inputs):
     f1 = chain.form_matrix(inputs[0])
     f2 = chain.form_matrix(inputs[1])
     dt_l, dt = inputs[2], inputs[3]
     return -dt / (12 * dt_l) * f1 + (dt + 6 * dt_l) / (12 * dt_l) * f2

def _leqi_f2(chain, inputs):
     f1 = chain.form_matrix(inputs[0])
     f2 = chain.form_matrix(inputs[1])
     dt_l, dt = inputs[2], inputs[3]
     return -5 * dt / (12 * dt_l) * f1 + (5 * dt + 6 * dt_l) / (12 * dt_l) * f2

def _leqi_f3(chain, inputs):
     f1 = chain.form_matrix(inputs[0])
     f2 = chain.form_matrix(inputs[1])
     f3 = chain.form_matrix(inputs[2])
     dt_l, dt = inputs[3], inputs[4]
     return (-dt**2 / (12 * dt_l * (dt + dt_l)) * f1 +
             (dt**2 + 6*dt*dt_l + 5*dt_l**2) / (12 * dt_l * (dt + dt_l)) * f2 +
             dt_l / (12 * (dt + dt_l)) * f3)

def _leqi_f4(chain, inputs):
     f1 = chain.form_matrix(inputs[0])
     f2 = chain.form_matrix(inputs[1])
     f3 = chain.form_matrix(inputs[2])
     dt_l, dt = inputs[3], inputs[4]
     return (-dt**2 / (12 * dt_l * (dt + dt_l)) * f1 +
             (dt**2 + 2*dt*dt_l + dt_l**2) / (12 * dt_l * (dt + dt_l)) * f2 +
             (4 * dt * dt_l + 5 * dt_l**2) / (12 * dt_l * (dt + dt_l)) * f3)

def leqi(operator, timesteps, power=None, power_density=None, print_out=True):
    r"""Deplete using the LE/QI CFQ4 algorithm.

    Implements the LE/QI Predictor-Corrector algorithm using the [fourth order
    commutator-free integrator]_.

    The LE/QI algorithm is mathematically defined as:

    .. math:
        y' = A(y, t) y(t)
        A_m1 = A(y_n-1, t_n-1)
        A_0 = A(y_n, t_n)
        A_l(t) linear extrapolation of A_m1, A_0
        Integrate to t_n+1 to get y_p
        A_c = A(y_p, y_n+1)
        A_q(t) quadratic interpolation of A_m1, A_0, A_c

    Here, A(t) is integrated using the fourth order algorithm CFQ4.

    It is initialized using the CE/LI algorithm.

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

    References
    ----------
    .. [fourth order commutator-free integrator]
       Thalhammer, Mechthild. "A fourth-order commutator-free exponential
       integrator for nonautonomous differential equations." SIAM journal on
       numerical analysis 44.2 (2006): 851-864.
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
            # LE/QI needs the last step results to start
            # Perform CE/LI CFQ4 or restore results for the first step
            if i == 0:
                if i_res <= 1:
                    dt_l = dt
                    x_new, t, op_res_last = celi_inner(operator, vec, p, i,
                                                       i_res, t, dt, print_out)
                    continue
                else:
                    dt_l = t - operator.prev_res[-2].time[0]
                    op_res_last = operator.prev_res[-2]
                    op_res_last.rates = op_res_last.rates[0]
                    x_new = operator.prev_res[-1].data[0]

            # Perform remaining LE/QI
            x = [copy.deepcopy(x_new)]
            op_results = [operator(x[0], p)]

            inputs = list(zip(op_res_last.rates, op_results[0].rates,
                              repeat(dt_l), repeat(dt)))
            x_new = deplete(chain, x[0], inputs, dt, print_out,
                            matrix_func=_leqi_f1)
            x_new = deplete(chain, x_new, inputs, dt, print_out,
                            matrix_func=_leqi_f2)
            x.append(x_new)
            op_results.append(operator(x[1], p))

            inputs = list(zip(op_res_last.rates, op_results[0].rates,
                              op_results[1].rates, repeat(dt_l), repeat(dt)))
            x_new = deplete(chain, x[0], inputs, dt, print_out,
                            matrix_func=_leqi_f3)
            x_new = deplete(chain, x_new, inputs, dt, print_out,
                            matrix_func=_leqi_f4)

            # Create results, write to disk
            Results.save(operator, x, op_results, [t, t+dt], p, i_res+i)

            # update results
            op_res_last = copy.deepcopy(op_results[0])
            t += dt
            dt_l = dt

        # Perform one last simulation
        x = [copy.deepcopy(x_new)]
        op_results = [operator(x[0], power[-1])]

        # Create results, write to disk
        Results.save(operator, x, op_results, [t, t], p, i_res + len(timesteps))
