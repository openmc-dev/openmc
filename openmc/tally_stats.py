import numpy as np
from math import sqrt, log
from scipy.stats import norm, chi2, shapiro

"""
[1] "A Suggestion for Using Powerful and Informative Tests of Normality"
Ralph B. D'Agostino, Albert Belanger, Ralph B. D'Agostino, Jr.
The American Statistician, Vol. 44, No. 4 (Nov., 1990), pp. 316-321 
"""

def _calc_b1_b2(n, mean, sum_sq, sum_rd, sum_th):
    """
    Compute sample skewness and kurtosis based on the tallies data.
    """
    m2 = (sum_sq/n -  mean**2)  
    m3 = 1.0/n * (sum_rd - 3*mean*sum_sq) + 2*mean**3
    m4 = 1.0/n * (sum_th - 4*mean*sum_rd + 6*mean**2*sum_sq) - 3*mean**4
    
    # sample skewness
    sqrt_b1 = np.where(m2 > 0, m3 / (m2 ** (1.5)), 0.0)
    #b1 = raw_skew ** 2.0

    # sample kurtosis
    b2 = np.where(m2 > 0, m4 / (m2 ** 2), 0.0)
    return sqrt_b1, b2

def skewness_test(n, mean, sum_sq, sum_rd, sum_th, alternative = "two-sided"):
    """
    Performs the D'Agostino-Pearson test for skewness [1].
    """

    sqrt_b1, _ = _calc_b1_b2(n, mean, sum_sq, sum_rd, sum_th)
 
    y = sqrt_b1 * np.sqrt(((n + 1.0) * (n + 3.0)) / (6.0 * (n - 2.0)))

    beta2 = (3.0 * (n**2 + 27.0*n - 70.0) * (n + 1.0) * (n + 3.0)) / ((n - 2.0) * (n + 5.0) * (n + 7.0) * (n + 9.0))
    W_squared = -1.0 + np.sqrt(2.0 * (beta2 - 1.0))
    delta = 1.0 / np.sqrt(np.log(np.sqrt(W_squared)))
    alpha = np.sqrt(2.0 / (W_squared - 1.0))
    
    Zb1 = np.where(y >= 0, delta * np.log((y / alpha) + np.sqrt((y / alpha) ** 2 + 1)),
         -delta * np.log((-y / alpha) + np.sqrt((y / alpha) ** 2 + 1))
    )
    
    if alternative == "two-sided":
        pval = 2.0 * (1.0 - norm.cdf(abs(Zb1)))
    elif alternative == "greater":
        pval = 1.0 - norm.cdf(Zb1)
    elif alternative == "less":
        pval = norm.cdf(Zb1)
    
    return Zb1, pval, sqrt_b1

def kurtosis_test(n, mean, sum_sq, sum_rd, sum_th, alternative = "two-sided"):
    """
    Performs the D'Agostino-Pearson test for kurtosis [1].
    """

    _, b2 = _calc_b1_b2(n, mean, sum_sq, sum_rd, sum_th)

    mean_b2 = 3.0 * (n - 1.0) / (n + 1.0)
    var_b2  = 24.0 * n * (n - 2.0) * (n - 3.0) / ((n + 1.0) ** 2 * (n + 3.0) * (n + 5.0))
    
    x = (b2 - mean_b2) / np.sqrt(var_b2)
    
    # Compute the third moment and transformation parameter A
    moment = ((6.0 * (n**2 - 5.0*n + 2.0)) / ((n + 7.0) * (n + 9.0))) * np.sqrt((6.0 * (n + 3.0) * (n + 5.0)) / ((n * (n - 2.0) * (n - 3.0))))
    A = 6.0 + (8.0 / moment) * ( (2.0 / moment) + np.sqrt(1.0 + 4.0 / (moment ** 2)) )
    Zb2 = (1.0 - 2.0 / (9.0 * A) - ((1.0 - 2.0 / A) / (1.0 + x * np.sqrt(2.0 / (A - 4.0)))) ** (1.0 / 3.0)) / np.sqrt(2.0 / (9.0 * A))
    
    if alternative == "two-sided":
        pval = 2.0 * (1.0 - norm.cdf(abs(Zb2)))
    elif alternative == "greater":
        pval = 1.0 - norm.cdf(Zb2)
    elif alternative == "less":
        pval = norm.cdf(Zb2)
    
    return Zb2, pval, b2

def k2_test(Zb1, Zb2):
    """
    Omnibus test of normality that combines skewness and kurtosis [1].
    """
    K2 = Zb1**2 + Zb2**2
    pval = 1.0 - chi2.cdf(K2, 2)
    return K2, pval

def print_normality_tests_summary(n, mean, sum_sq, sum_rd, sum_th, sig_level,
                                  skew_alternative = "two-sided",
                                  kurt_alternative = "two-sided"):
    """
    Summary of the skewness and kurtosis tests.
    """
    # Ensure inputs are NumPy arrays for consistent handling
    mean = np.asarray(mean)
    sum_sq = np.asarray(sum_sq)
    sum_rd = np.asarray(sum_rd)
    sum_th = np.asarray(sum_th)
    
    # Determine the number of filters (assumes the first dimension corresponds to filters)
    num_filters = mean.shape[0] if mean.ndim > 1 else 1

    for i in range(num_filters):
        
        # Extract data for the current filter
        current_mean = mean[i] if num_filters > 1 else mean
        current_sum_sq = sum_sq[i] if num_filters > 1 else sum_sq
        current_sum_rd = sum_rd[i] if num_filters > 1 else sum_rd
        current_sum_th = sum_th[i] if num_filters > 1 else sum_th

        # Perform skewness and kurtosis tests
        Zb1, p_skew, sqrt_b1 = skewness_test(
            n, current_mean, current_sum_sq, current_sum_rd, current_sum_th, alternative=skew_alternative
        )
        Zb2, p_kurt, b2 = kurtosis_test(
            n, current_mean, current_sum_sq, current_sum_rd, current_sum_th, alternative=kurt_alternative
        )
        K2, p_omni = k2_test(Zb1, Zb2)

        # Print skewness results
        print("\nSkewness Test Result:")
        print(f"Raw skewness: {sqrt_b1}")
        print(f"Z statistic: {Zb1}, p-value: {p_skew}")
        if p_skew < sig_level:
            if np.all(sqrt_b1 > 0):
                print("Data is right-skewed")
            elif np.all(sqrt_b1 < 0):
                print("Data is left-skewed")
        else:
            print("Data is symmetric.")

        # Print kurtosis results
        print("\nKurtosis Test Result:")
        print(f"Kurtosis: {b2}")
        print(f"Z statistic: {Zb2}, p-value: {p_kurt}")
        if p_kurt < sig_level:
            if np.all(b2 > (3*(n-1)/(n+1)) ):
                print("Distribution is leptokurtic")
            elif np.all(b2 < (3*(n-1)/(n+1))):
                print("Distribution is platykurtic")
        else:
            print("Distribution has normal kurtosis.")

        print("\nOmnibus Test Result:")
        print(f"K2 statistic: {K2}, p-value: {p_omni}")
        if p_omni < sig_level:
            print("Data is not normally distributed.")
        else:
            print("Data is normally distributed.")