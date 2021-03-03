from collections.abc import Callable
from numbers import Real

import scipy.optimize as sopt

import openmc
import openmc.model
import openmc.checkvalue as cv


_SCALAR_BRACKETED_METHODS = ['brentq', 'brenth', 'ridder', 'bisect']


def _search_keff(guess, target, model_builder, model_args, print_iterations,
                 print_output, guesses, results):
    """Function which will actually create our model, run the calculation, and
    obtain the result. This function will be passed to the root finding
    algorithm

    Parameters
    ----------
    guess : Real
        Current guess for the parameter to be searched in `model_builder`.
    target_keff : Real
        Value to search for
    model_builder : collections.Callable
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
    model_args : dict
        Keyword-based arguments to pass to the `model_builder` method.
    print_iterations : bool
        Whether or not to print the guess and the resultant keff during the
        iteration process.
    print_output : bool
        Whether or not to print the OpenMC output during the iterations.
    guesses : Iterable of Real
        Running list of guesses thus far, to be updated during the execution of
        this function.
    results : Iterable of Real
        Running list of results thus far, to be updated during the execution of
        this function.

    Returns
    -------
    float
        Value of the model for the current guess compared to the target value.

    """

    # Build the model
    model = model_builder(guess, **model_args)

    # Run the model and obtain keff
    sp_filepath = model.run(output=print_output)
    with openmc.StatePoint(sp_filepath) as sp:
        keff = sp.k_combined

    # Record the history
    guesses.append(guess)
    results.append(keff)

    if print_iterations:
        text = 'Iteration: {}; Guess of {:.2e} produced a keff of ' + \
            '{:1.5f} +/- {:1.5f}'
        print(text.format(len(guesses), guess, keff.n, keff.s))

    return keff.n - target


def search_for_keff(model_builder, initial_guess=None, target=1.0,
                    bracket=None, model_args=None, tol=None,
                    bracketed_method='bisect', print_iterations=False,
                    print_output=False, **kwargs):
    """Function to perform a keff search by modifying a model parametrized by a
    single independent variable.

    Parameters
    ----------
    model_builder : collections.Callable
        Callable function which builds a model according to a passed
        parameter. This function must return an openmc.model.Model object.
    initial_guess : Real, optional
        Initial guess for the parameter to be searched in
        `model_builder`. One of `guess` or `bracket` must be provided.
    target : Real, optional
        keff value to search for, defaults to 1.0.
    bracket : None or Iterable of Real, optional
        Bracketing interval to search for the solution; if not provided,
        a generic non-bracketing method is used. If provided, the brackets
        are used. Defaults to no brackets provided. One of `guess` or `bracket`
        must be provided. If both are provided, the bracket will be
        preferentially used.
    model_args : dict, optional
        Keyword-based arguments to pass to the `model_builder` method. Defaults
        to no arguments.
    tol : float
        Tolerance to pass to the search method
    bracketed_method : {'brentq', 'brenth', 'ridder', 'bisect'}, optional
        Solution method to use; only applies if
        `bracket` is set, otherwise the Secant method is used.
        Defaults to 'bisect'.
    print_iterations : bool
        Whether or not to print the guess and the result during the iteration
        process. Defaults to False.
    print_output : bool
        Whether or not to print the OpenMC output during the iterations.
        Defaults to False.
    **kwargs
        All remaining keyword arguments are passed to the root-finding
        method.

    Returns
    -------
    zero_value : float
        Estimated value of the variable parameter where keff is the
        targeted value
    guesses : List of Real
        List of guesses attempted by the search
    results : List of 2-tuple of Real
        List of keffs and uncertainties corresponding to the guess attempted by
        the search

    """

    if initial_guess is not None:
        cv.check_type('initial_guess', initial_guess, Real)
    if bracket is not None:
        cv.check_iterable_type('bracket', bracket, Real)
        cv.check_length('bracket', bracket, 2)
        cv.check_less_than('bracket values', bracket[0], bracket[1])
    if model_args is None:
        model_args = {}
    else:
        cv.check_type('model_args', model_args, dict)
    cv.check_type('target', target, Real)
    cv.check_type('tol', tol, Real)
    cv.check_value('bracketed_method', bracketed_method,
                   _SCALAR_BRACKETED_METHODS)
    cv.check_type('print_iterations', print_iterations, bool)
    cv.check_type('print_output', print_output, bool)
    cv.check_type('model_builder', model_builder, Callable)

    # Run the model builder function once to make sure it provides the correct
    # output type
    if bracket is not None:
        model = model_builder(bracket[0], **model_args)
    elif initial_guess is not None:
        model = model_builder(initial_guess, **model_args)
    cv.check_type('model_builder return', model, openmc.model.Model)

    # Set the iteration data storage variables
    guesses = []
    results = []

    # Set the searching function (for easy replacement should a later
    # generic function be added.
    search_function = _search_keff

    if bracket is not None:
        # Generate our arguments
        args = {'f': search_function, 'a': bracket[0], 'b': bracket[1]}
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

    elif initial_guess is not None:

        # Generate our arguments
        args = {'func': search_function, 'x0': initial_guess}
        if tol is not None:
            args['tol'] = tol

        # Set the root finding method
        root_finder = sopt.newton

    else:
        raise ValueError("Either the 'bracket' or 'initial_guess' parameters "
                         "must be set")

    # Add information to be passed to the searching function
    args['args'] = (target, model_builder, model_args, print_iterations,
                    print_output, guesses, results)

    # Create a new dictionary with the arguments from args and kwargs
    args.update(kwargs)

    # Perform the search
    zero_value = root_finder(**args)

    return zero_value, guesses, results
