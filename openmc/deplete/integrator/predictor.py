"""First-order predictor algorithm."""

from .abc import Integrator
from .cram import timed_deplete


class PredictorIntegrator(Integrator):
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
        Array of timesteps in units of [s]. Note that values are not
        cumulative.
    power : float or iterable of float, optional
        Power of the reactor in [W]. A single value indicates that
        the power is constant over all timesteps. An iterable
        indicates potentially different power levels for each timestep.
        For a 2D problem, the power can be given in [W/cm] as long
        as the "volume" assigned to a depletion material is actually
        an area in [cm^2]. Either ``power`` or ``power_density`` must be
        specified.
    power_density : float or iterable of float, optional
        Power density of the reactor in [W/gHM]. It is multiplied by
        initial heavy metal inventory to get total power if ``power``
        is not speficied.
    """
    _N_STAGES = 1

    def __call__(self, conc, rates, dt, power, _i=None):
        """Perform the integration across one time step

        Parameters
        ----------
        conc : numpy.ndarray
            Initial concentrations for all nuclides in [atom]
        rates : openmc.deplete.ReactionRates
            Reaction rates from operator
        dt : float
            Time in [s] for the entire depletion interval
        power : float
            Power of the system in [W]
        _i : int or None
            Iteration index. Not used

        Returns
        -------
        proc_time : float
            Time spent in CRAM routines for all materials in [s]
        conc_list : list of numpy.ndarray
            Concentrations at end of interval
        op_results : empty list
            Kept for consistency with API. No intermediate calls to
            operator with predictor

        """
        proc_time, conc_end = timed_deplete(self.chain, conc, rates, dt)
        return proc_time, [conc_end], []
