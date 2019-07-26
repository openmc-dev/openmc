"""The SI-CE/LI CFQ4 integrator."""

import copy

from uncertainties import ufloat

from .abc import SI_Integrator
from .cram import timed_deplete
from .celi import _celi_f1, _celi_f2
from ..abc import OperatorResult
from ..results import Results


class SI_CELI_Integrator(SI_Integrator):
    r"""Deplete using the si-ce/li cfq4 algorithm.

    Implements the stochastic implicit ce/li predictor-corrector algorithm using
    the `fourth order commutator-free integrator <https://doi.org/10.1137/05063042>`_.

    detailed algorithm can be found in section 3.2 in `colin josey's thesis
    <http://hdl.handle.net/1721.1/113721>`_.
    """
    def __call__(self, bos_conc, bos_rates, dt, power, _i):
        """Perform the integration across one time step

        Parameters
        ----------
        bos_conc : numpy.ndarray
            Initial bos_concentrations for all nuclides in [atom]
        bos_rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system [W]
        _i : int
            Current depletion step index. Unused

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials
        bos_conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final bos_concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from intermediate transport
            simulations
        """
        proc_time, eos_conc = timed_deplete(
            self.chain, bos_conc, bos_rates, dt)
        inter_conc = copy.deepcopy(eos_conc)

        # Begin iteration
        for j in range(self.n_stages + 1):
            inter_res = self.operator(inter_conc, power)

            if j <= 1:
                res_bar = copy.deepcopy(inter_res)
            else:
                rates = 1/j * inter_res.rates + (1 - 1 / j) * res_bar.rates
                k = 1/j * inter_res.k + (1 - 1 / j) * res_bar.k
                res_bar = OperatorResult(k, rates)

            list_rates = list(zip(bos_rates, res_bar.rates))
            time1, inter_conc = timed_deplete(
                self.chain, bos_conc, list_rates, dt, matrix_func=_celi_f1)
            time2, inter_conc = timed_deplete(
                self.chain, inter_conc, list_rates, dt, matrix_func=_celi_f2)
            proc_time += time1 + time2

        # end iteration
        return proc_time, [eos_conc, inter_conc], [res_bar]


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
