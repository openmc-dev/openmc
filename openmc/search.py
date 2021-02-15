from collections.abc import Callable
from math import ceil
from numbers import Real
import openmc
import openmc.model
import openmc.checkvalue as cv


def search_for_keff(model_builder, guesses, target=1.0,
                    kwargs={}, tolerance=1e-3):

    _run_safety_checks(model_builder, guesses, target,
                       kwargs, tolerance)

    batches = model_builder(guesses[0], **kwargs).settings.batches

    print("Solving your guesses...\n")

    z_old = guesses[0]
    z_new = guesses[1]
    k_old = _calculate_keff(z_old, model_builder, kwargs, batches)
    k_new = _calculate_keff(z_new, model_builder, kwargs, batches)

    keffs = [k_old, k_new]

    print("Beginning search sequence...\n")

    while not _critical(k_new, target, tolerance):

        z_next = _make_guess(z_new, z_old, k_new, k_old, target)

        z_old = z_new
        k_old = k_new
        z_new = z_next
        delta = abs(z_new-z_old)

        print("guess = {}".format(z_new))

        k_new = _calculate_keff(z_new, model_builder,
                                kwargs, batches)

        while delta <= _sigma_z(z_new, z_old, k_new, k_old, target):

            batches = _increase_batches(batches, k_new, k_old)

            k_new = _calculate_keff(
                z_new, model_builder, kwargs, batches)

    _store_results(z_new, k_new, guesses, keffs)

    return guesses, keffs


def _calculate_keff(guess, model_builder, kwargs, batches):

    model = model_builder(guess, **kwargs)
    model.settings.batches = batches
    sp = model.run(output=False)
    k = openmc.StatePoint(sp).k_combined

    print("keff  = {}\n".format(k))

    return k


def _make_guess(z_new, z_old, k_new, k_old, target):
    DW = (k_new.n - k_old.n)/(z_new - z_old)
    z_next = z_new + (target - k_new.n)/DW
    return z_next


def _increase_batches(batches, k_new, k_old):

    # Optional
    # This function implements Morrow's reccomended increase.
    # However much of the time k_new.s is smaller than k_old.s
    # IMO this increases code complexity with little benefit.

    uncertainty_ratio = k_old.s/k_new.s

    if uncertainty_ratio > 1:
        batches *= ceil(uncertainty_ratio**2)
    else:
        batches *= 2

    print("increasing batches to {}".format(batches))

    return batches


def _critical(k, target, tolerance):
    return abs(k - target) < tolerance


def _sigma_z(z_new, z_old, k_new, k_old, target):

    LHS_numerator = (z_new - z_old) * (k_old.n - target)
    RHS_numerator = (z_new - z_old) * (target - k_new.n)
    denominator = (k_new.n - k_old.n)**2

    LHS = k_new.s**2 * (LHS_numerator/denominator)**2
    RHS = k_old.s**2 * (RHS_numerator/denominator)**2

    sigma_z = (LHS + RHS)**0.5

    return sigma_z


def _store_results(z_new, k_new, guesses, keffs):
    guesses.append(z_new)
    keffs.append(k_new)
    return None


def _run_safety_checks(model_builder, guesses, target,
                       model_kwargs, tolerance):

    print("Checking inputs are correctly formatted...\n")

    cv.check_iterable_type("guesses", guesses, Real)
    cv.check_type("target", target, Real)
    cv.check_type("tolerance", tolerance, Real)
    cv.check_type("model_builder", model_builder, Callable)
    cv.check_type("model_kwargs", model_kwargs, dict)
    cv.check_type("model_builder", model_builder(
        guesses[0], **model_kwargs), openmc.model.Model)

    return None
