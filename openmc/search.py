from numbers import Real
from types import FunctionType
from inspect import getfullargspec

import openmc
from openmc.checkvalue import check_type, check_length


class KeffSearch(object):
    """Class to perform a search for a certain keff value given changes on an
    arbitrary scalar input.

    Parameters
    ----------
    model_builder : FunctionType
        Callable function which builds a model according to a single, passed
        parameter. This function must return an openmc.ModelContainer object.
    guess : Real
        Initial guess for the parameter modified in :param:`model_builder`.
    target_keff : Real
        keff value to search for, defaults to 1.0.

    Attributes
    ----------
    model_builder : FunctionType
        Callable function which builds a model according to a single, passed
        parameter. This function must return an openmc.ModelContainer object.
    guess : Real
        Initial guess for the parameter modified in :param:`model_builder`.
    target_keff : Real
        keff value to search for.
    guesses : List of Real
        List of guesses attempted by the search
    keffs : List of Real
        List of keffs corresponding to the guess attempted by the search
    keff_uncs : List of Real
        List of keff uncertainties corresponding to the guess attempted by the
        search

    """

    def __init__(self, model_builder, guess, target_keff=1.0):
        self.initial_guess = guess
        self.model_builder = model_builder
        self.target_keff = target_keff
        self.guesses = []
        self.keffs = []
        self.keff_uncs = []

    @property
    def model_builder(self):
        return self._model_builder

    @model_builder.setter
    def model_builder(self, model_builder):
        # Make sure model_builder is a function
        check_type('model_builder', model_builder, FunctionType)

        # Make sure model_builder has only one parameter
        argspec = getfullargspec(model_builder)
        check_length('model_builder arguments', argspec.args, 1)

        # Run the model builder function once to make sure it provides the
        # correct output types
        model = model_builder(self.initial_guess)
        check_type('model_builder return', model, openmc.ModelContainer)
        self._model_builder = model_builder

    @property
    def initial_guess(self):
        return self._initial_guess

    @initial_guess.setter
    def initial_guess(self, initial_guess):
        self._initial_guess = initial_guess

    @property
    def target_keff(self):
        return self._target_keff

    @target_keff.setter
    def target_keff(self, target_keff):
        check_type('target_keff', target_keff, Real)
        self._target_keff = target_keff

    def _search_function(self, guess):
        # Build the model
        model = self.model_builder(guess)

        # Run the model
        keff = model.execute(output=False)

        # Close the model to ensure HDF5 will allow access during the next
        # OpenMC execution
        model.close()

        # Record the history
        self.guesses.append(guess)
        self.keffs.append(keff[0])
        self.keff_uncs.append(keff[1])

        return (keff[0] - self.target_keff)

    def search(self, **kwargs):
        """Searches for the target eigenvalue with the Newton-Raphson method

        Parameters
        ----------
        **kwargs
            Keyword arguments passed to :func:`scipy.optimize.newton`.

        Returns
        zero_value : Real
            Estimated value of the variable parameter where keff is the
            targeted value
        keff : Iterable of Real
            keff calculated at the zero_value

        """

        import scipy.optimize as sopt

        zero_value = sopt.newton(self._search_function, self.initial_guess,
                                 **kwargs)

        return zero_value
