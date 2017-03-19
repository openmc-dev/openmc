from numbers import Real
from types import FunctionType

import openmc
import openmc.model
from openmc.checkvalue import check_type, check_iterable_type, check_length, \
    check_value


_SCALAR_BRACKETED_METHODS = ['brentq', 'brentq', 'ridder', 'bisect']


class KeffSearch(object):
    """Class to perform a keff search by modifying an arbitrarily parametrized
    model.

    Parameters
    ----------
    model_builder : FunctionType
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
    guess : Real, optional
        Initial guess for the parameter to be searched in
        :param:`model_builder`. One of :param:`guess` or :param`bracket` must
        be provided.
    target_keff : Real, optional
        keff value to search for, defaults to 1.0.
    bracket : None or Iterable of Real, optional
        Bracketing interval to search for the solution; if not provided,
        a generic non-bracketing method is used. If provided, the brackets
        are used. Defaults to no brackets provided. One of :param:`guess` or
        :param`bracket` must be provided. If both are provided, the bracket
        will be preferentially used.

    Attributes
    ----------
    model_builder : FunctionType
        Callable function which builds a model according to parameters passed.
        This function must return an openmc.model.Model
        object.
    initial_guess : Iterable of Real
        Initial guess for the parameter to be searched in
        :param:`model_builder`.
    target_keff : Real
        keff value to search for.
    bracket : None or Iterable of Real
        Bracketing interval to search for the solution; if not provided,
        a generic non-bracketing method is used. If provided, the brackets
        are used.
    print_iterations : bool
        Whether or not to print the guess and the resultant keff during the
        iteration process. Defaults to False.
    print_output : bool
        Whether or not to print the OpenMC output during the iterations.
        Defaults to False.
    guesses : List of Real
        List of guesses attempted by the search
    keffs : List of Real
        List of keffs corresponding to the guess attempted by the search
    keff_uncs : List of Real
        List of keff uncertainties corresponding to the guess attempted by the
        search

    """

    def __init__(self, model_builder, guess=None, bracket=None,
                 target_keff=1.0):
        self.initial_guess = guess
        self.bracket = bracket
        self.model_builder = model_builder
        self.target_keff = target_keff
        self.guesses = []
        self.keffs = []
        self.keff_uncs = []
        self.print_iterations = False
        self.print_output = False

    @property
    def model_builder(self):
        return self._model_builder

    @model_builder.setter
    def model_builder(self, model_builder):
        # Make sure model_builder is a function
        check_type('model_builder', model_builder, FunctionType)

        # Run the model builder function once to make sure it provides the
        # correct output type
        if self.bracket is not None:
            model = model_builder(self.bracket[0])
        elif self.initial_guess is not None:
            model = model_builder(self.initial_guess)
        check_type('model_builder return', model, openmc.model.Model)
        self._model_builder = model_builder

    @property
    def initial_guess(self):
        return self._initial_guess

    @initial_guess.setter
    def initial_guess(self, initial_guess):
        if initial_guess is not None:
            check_type('initial_guess', initial_guess, Real)
        self._initial_guess = initial_guess

    @property
    def bracket(self):
        return self._bracket

    @bracket.setter
    def bracket(self, bracket):
        if bracket is not None:
            check_iterable_type('bracket', bracket, Real)
            check_length('bracket', bracket, 2)
        self._bracket = bracket

    @property
    def target_keff(self):
        return self._target_keff

    @target_keff.setter
    def target_keff(self, target_keff):
        check_type('target_keff', target_keff, Real)
        self._target_keff = target_keff

    @property
    def print_iterations(self):
        return self._print_iterations

    @print_iterations.setter
    def print_iterations(self, print_iterations):
        check_type('print_iterations', print_iterations, bool)
        self._print_iterations = print_iterations

    @property
    def print_output(self):
        return self._print_output

    @print_output.setter
    def print_output(self, print_output):
        check_type('print_output', print_output, bool)
        self._print_output = print_output

    def _search_function(self, guess):
        # Build the model
        model = self.model_builder(guess)

        # Run the model
        keff = model.execute(output=self._print_output)

        # Close the model to ensure HDF5 will allow access during the next
        # OpenMC execution
        model.close()

        # Record the history
        self.guesses.append(guess)
        self.keffs.append(keff[0])
        self.keff_uncs.append(keff[1])

        if self._print_iterations:
            print(guess, '{:1.5f} +/- {:1.5f}'.format(keff[0], keff[1]))

        return (keff[0] - self.target_keff)

    def search(self, tol=None, bracketed_method='brentq', **kwargs):
        """Searches for the target eigenvalue with the Newton-Raphson method

        Parameters
        ----------
        tol : Real
            Tolerance to pass to the search method
        bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
            Solution method to use; only applies if
            :param:`bracket` is set, otherwise the Newton method is used.
            Defaults to 'brentq'.
        **kwargs
            Keyword arguments passed to :func:`scipy.optimize.newton`.

        Returns
        zero_value : Real
            Estimated value of the variable parameter where keff is the
            targeted value
        keff : Iterable of Real
            keff calculated at the zero_value

        """

        check_value('bracketed_method', bracketed_method,
                    _SCALAR_BRACKETED_METHODS)

        import scipy.optimize as sopt

        if self.bracket is not None:
            # Generate our arguments
            args = {'f': self._search_function, 'a': self.bracket[0],
                    'b': self.bracket[1]}
            if tol is not None:
                args['rtol'] = tol

            # Set the root finding method
            if bracketed_method == 'brentq':
                root_finder = sopt.brentq
            elif bracketed_method == 'brenth':
                root_finder = sopt.brenth
            elif bracketed_method == 'ridder':
                root_finder = sopt.ridder
            elif bracketed_method == 'bisect':
                root_finder = sopt.bisect

        elif self.initial_guess is not None:

            # Generate our arguments
            args = {'func': self._search_function, 'x0': self.initial_guess}
            if tol is not None:
                args['tol'] = tol

            # Set the root finding method
            root_finder = sopt.newton

        else:
            raise ValueError("One of the 'bracket' or 'initial_guess' "
                             "parameters must be set")

        # Perform the search
        zero_value = root_finder(**args, **kwargs)

        return zero_value
