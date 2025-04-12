import numpy as np
from math import sqrt, log
from scipy.stats import norm, chi2, shapiro
import matplotlib.pyplot as plt

def _calc_b1_b2(n, mean, sum_sq, sum_rd, sum_th):
    """
    Compute raw skewness, sample skewness (b1) and kurtosis (b2) based on
    the central moments of the tallies data.
    """
    m2 = (sum_sq/n -  mean**2)  
    m3 = 1.0/n * (sum_rd - 3*mean*sum_sq) + 2*mean**3
    m4 = 1.0/n * (sum_th - 4*mean*sum_rd + 6*mean**2*sum_sq) - 3*mean**4
    
    raw_skew = np.where(m2 > 0, m3 / (m2 ** (1.5)), 0.0)
    b1 = raw_skew ** 2
    b2 = np.where(m2 > 0, m4 / (m2 ** 2), 0.0)
    return raw_skew, b1, b2

def dagostino_pearson_skew_test(n, mean, sum_sq, sum_rd, sum_th, alternative = "two-sided"):
    """
    Performs the D'Agostino-Pearson test for skewness.
    """

    raw_skew, _, _ = _calc_b1_b2(n, mean, sum_sq, sum_rd, sum_th)
 
    y = raw_skew * np.sqrt(((n + 1.0) * (n + 3.0)) / (6.0 * (n - 2.0)))

    beta2 = (3.0 * (n**2 + 27.0*n - 70.0) * (n + 1.0) * (n + 3.0)) / ((n - 2.0) * (n + 5.0) * (n + 7.0) * (n + 9.0))
    W_squared = -1.0 + np.sqrt(2.0 * (beta2 - 1.0))
    delta = 1.0 / np.sqrt(np.log(np.sqrt(W_squared)))
    alpha = np.sqrt(2.0 / (W_squared - 1.0))
    
    Zb1 = np.where(
        y >= 0,
        delta * np.log((y / alpha) + np.sqrt((y / alpha) ** 2 + 1)),
         -delta * np.log((-y / alpha) + np.sqrt((y / alpha) ** 2 + 1))
    )
    
    if alternative == "two-sided":
        pval = 2.0 * (1.0 - norm.cdf(abs(Zb1)))
    elif alternative == "greater":
        pval = 1.0 - norm.cdf(Zb1)
    elif alternative == "less":
        pval = norm.cdf(Zb1)
    
    return Zb1, pval, raw_skew

def dagostino_pearson_kurt_test(n, mean, sum_sq, sum_rd, sum_th, alternative = "two-sided"):
    """
    Performs the D'Agostino-Pearson test for kurtosis.
    """

    _, _, b2 = _calc_b1_b2(n, mean, sum_sq, sum_rd, sum_th)

    mean_b2 = 3.0 * (n - 1.0) / (n + 1.0)
    var_b2  = 24.0 * n * (n - 2.0) * (n - 3.0) / ((n + 1.0) ** 2 * (n + 3.0) * (n + 5.0))
    
    x = (b2 - mean_b2) / np.sqrt(var_b2)
    
    # Compute the third moment and transformation parameter A
    moment = (6.0 * (n**2 - 5.0*n + 2.0)) / ((n + 7.0) * (n + 9.0)) * np.sqrt(6.0 * (n + 3.0) * (n + 5.0) / (n * (n - 2.0) * (n - 3.0)))
    A = 6.0 + (8.0 / moment) * ( (2.0 / moment) + np.sqrt(1.0 + 4.0 / (moment ** 2)) )
    Zb2 = (1.0 - 2.0 / (9.0 * A) -
           ((1.0 - 2.0 / A) / (1.0 + x * np.sqrt(2.0 / (A - 4.0)))) ** (1.0 / 3.0)) / np.sqrt(2.0 / (9.0 * A))
    
    if alternative == "two-sided":
        pval = 2.0 * (1.0 - norm.cdf(abs(Zb2)))
    elif alternative == "greater":
        pval = 1.0 - norm.cdf(Zb2)
    elif alternative == "less":
        pval = norm.cdf(Zb2)
    
    return Zb2, pval, b2

def dagostino_pearson_k2_test(Zb1, Zb2):
    """
    Omnibus test of normality that combines skewness and kurtosis.
    """
    K2 = Zb1**2 + Zb2**2
    pval = 1.0 - chi2.cdf(K2, 2)
    return K2, pval

def print_normality_tests_summary(n, mean, sum_sq, sum_rd, sum_th,
                                  skew_alternative = "two-sided",
                                  kurt_alternative = "two-sided"):
    """
    Computes and prints a summary of the skewness and kurtosis tests,
    including interpretations.
    """
    # Ensure inputs are NumPy arrays for consistent handling
    mean = np.asarray(mean)
    sum_sq = np.asarray(sum_sq)
    sum_rd = np.asarray(sum_rd)
    sum_th = np.asarray(sum_th)
    
    # Determine the number of filters (assumes the first dimension corresponds to filters)
    num_filters = mean.shape[0] if mean.ndim > 1 else 1

    for i in range(num_filters):
        print(f"\n--- Filter {i + 1} ---")

        # Extract data for the current filter
        current_mean = mean[i] if num_filters > 1 else mean
        current_sum_sq = sum_sq[i] if num_filters > 1 else sum_sq
        current_sum_rd = sum_rd[i] if num_filters > 1 else sum_rd
        current_sum_th = sum_th[i] if num_filters > 1 else sum_th

        # Perform skewness and kurtosis tests
        Zb1, p_skew, raw_skew = dagostino_pearson_skew_test(
            n, current_mean, current_sum_sq, current_sum_rd, current_sum_th, alternative=skew_alternative
        )
        Zb2, p_kurt, b2 = dagostino_pearson_kurt_test(
            n, current_mean, current_sum_sq, current_sum_rd, current_sum_th, alternative=kurt_alternative
        )
        K2, p_omni = dagostino_pearson_k2_test(Zb1, Zb2)

        # Print skewness results
        print("\nSkewness Test Result:")
        print(f"Raw skewness: {raw_skew}")
        print(f"Z statistic: {Zb1}, p-value: {p_skew}")
        if np.all(raw_skew > 0):
            print("Interpretation: Data is right-skewed (long tail on the right).")
        elif np.all(raw_skew < 0):
            print("Interpretation: Data is left-skewed (long tail on the left).")
        else:
            print("Interpretation: Data is symmetric.")

        # Print kurtosis results
        print("\nKurtosis Test Result:")
        print(f"Kurtosis: {b2}")
        print(f"Z statistic: {Zb2}, p-value: {p_kurt}")
        if np.all(b2 > 3):
            print("Interpretation: Distribution is leptokurtic (more peaked with heavy tails than normal).")
        elif np.all(b2 < 3):
            print("Interpretation: Distribution is platykurtic (flatter with lighter tails than normal).")
        else:
            print("Interpretation: Distribution has normal kurtosis.")

        print("\nOmnibus Test Result:")
        print(f"K2 statistic: {K2}, p-value: {p_omni}")
        if p_omni < 0.05:
            print("Interpretation: Data is not normally distributed.")
        else:
            print("Interpretation: Data is normally distributed.")
