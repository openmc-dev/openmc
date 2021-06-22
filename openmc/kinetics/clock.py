import numpy as np

TIME_POINTS = ['START',
               'PREVIOUS_OUTER',
               'FORWARD_OUTER',
               'PREVIOUS_INNER',
               'FORWARD_INNER',
               'END']


class Clock:
    """Class to track the time evolution of the transient.

    Parameters
    ----------
    dt_inner : float
        Time step size for inner point kinetics calculations.
    t_outer : numpy array
        Time steps for outer Monte Carlo calculations.

    Attributes
    ----------
    dt_inner : float
        Time step size for inner point kinetics calculations.
    dt_outer : float
        Time steps size for outer Monte Carlo calculations.
    num_outer_steps : int
        Number of outer time steps.
    outer_step : int
        Counter that tracks how many outer steps have been taken.
    times : dict
        Relative time points corresponding to outer time steps.
    t_outer : numpy array
        Time steps for outer Monte Carlo calculations.

    """

    def __init__(self, dt_inner, t_outer):

        # Initialize coordinates
        self.dt_inner = dt_inner
        self.t_outer = t_outer

    def __repr__(self):

        string = 'Clock\n'
        string += '{0: <24}{1}{2}\n'.format('\tdt inner', '=\t', self.dt_inner)
        string += '{0: <24}{1}{2}\n'.format('\tt outer', '=\t', self.t_outer)

        for t in TIME_POINTS:
            string += '{0: <24}{1}{2}\n'.format('\tTime ' + t, '=\t', self.times[t])

        return string

    @property
    def t_outer(self):
        return self._t_outer

    @property
    def times(self):
        return self._times

    @property
    def outer_step(self):
        return self._outer_step

    @property
    def num_outer_steps(self):
        return len(self.t_outer) - 1

    @property
    def dt_inner(self):
        return self._dt_inner

    @dt_inner.setter
    def dt_inner(self, dt_inner):
        self._dt_inner = float(dt_inner)

    @t_outer.setter
    def t_outer(self, t_outer):
        self._t_outer = np.float64(t_outer)
        
        # Create a dictionary of clock times
        self._times = {}
        for t in TIME_POINTS:
            self._times[t] = t_outer[0]

        # Reset the end time
        self._times['END'] = t_outer[-1]
        self._outer_step = 0

    @property
    def dt_outer(self):
        return self.t_outer[self.outer_step + 1] - self.t_outer[self.outer_step]

    @outer_step.setter
    def outer_step(self, outer_step):
        self._outer_step = outer_step
