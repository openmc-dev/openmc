"""The CE/CM integrator."""

from .abc import Integrator
from .cram import timed_deplete


class CECMIntegrator(Integrator):
    r"""Deplete using the CE/CM algorithm.

    Implements the second order `CE/CM predictor-corrector algorithm
    <https://doi.org/10.13182/NSE14-92>`_.

    "CE/CM" stands for constant extrapolation on predictor and constant
    midpoint on corrector. This algorithm is mathematically defined as:

    .. math::
        \begin{aligned}
        y' &= A(y, t) y(t) \\
        A_p &= A(y_n, t_n) \\
        y_m &= \text{expm}(A_p h/2) y_n \\
        A_c &= A(y_m, t_n + h/2) \\
        y_{n+1} &= \text{expm}(A_c h) y_n
        \end{aligned}
    """
    _num_stages = 2

    def __call__(self, conc, rates, dt, power, _i=-1):
        """Integrate using CE/CM

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system [W]
        _i : int, optional
            Current iteration count. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at each of the intermediate points with
            the final concentration as the last element
        op_results : list of openmc.deplete.OperatorResult
            Eigenvalue and reaction rates from transport simulations
        """
        # deplete across first half of inteval
        time0, x_middle = timed_deplete(self.chain, conc, rates, dt / 2)
        res_middle = self.operator(x_middle, power)

        # deplete across entire interval with BOS concentrations,
        # MOS reaction rates
        time1, x_end = timed_deplete(self.chain, conc, res_middle.rates, dt)

        return time0 + time1, [x_middle, x_end], [res_middle]
