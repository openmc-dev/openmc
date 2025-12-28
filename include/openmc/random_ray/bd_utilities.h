#ifndef OPENMC_RANDOM_RAY_BD_UTILITIES_H
#define OPENMC_RANDOM_RAY_BD_UTILITIES_H

#include <algorithm>
#include <deque>
#include <map>

#include "openmc/error.h"
#include "openmc/vector.h"

namespace openmc {

//----------------------------------------------------------------------------
// Helper Variables
// Coefficients come from Table 3 in Fornberg (1988)
// DOI: 10.1090/S0025-5718-1988-0935077-0
// Note that the signs are flipped compared to the citation, as the author was
// formulating weights for a forward difference
const std::map<int, vector<double>> bd_coefficients_first_order_ = {
  {1, {1.0, -1.0}}, {2, {1.5, -2.0, 0.5}},
  {3, {1.833333333333333, -3.0, 1.5, -0.333333333333333}},
  {4, {2.083333333333333, -4.0, 3.0, -1.333333333333333, 0.25}},
  {5, {2.283333333333333, -5.0, 5.0, -3.333333333333333, 1.25, -0.2}},
  {6, {2.45, -6.0, 7.5, -6.666666666666667, 3.75, -1.2, 0.166666666666667}}};

// Coefficients come from Table 3 in Fornberg (1988)
// DOI: 10.1090/S0025-5718-1988-0935077-0
const std::map<int, vector<double>> bd_coefficients_second_order_ = {
  {1, {1.0, -2.0, 1.0}}, {2, {2.0, -5.0, 4, -1}},
  {3, {2.916666666666667, -8.666666666666667, 9.5, -4.666666666666667,
        0.916666666666667}},
  {4, {3.75, -12.833333333333333, 17.833333333333333, -13.0, 5.083333333333333,
        -0.833333333333333}},
  {5, {4.511111111111111, -17.4, 29.25, -28.222222222222222, 16.5, -5.4,
        0.761111111111111}},
  {6, {5.211111111111111, -22.3, 43.95, -52.722222222222222, 41.0, -20.1,
        5.661111111111111, -0.7}}};

// Take RHS derivative to solve for the current timestep
template<typename T>
T rhs_backwards_difference(
  std::deque<T>& bd_vector, int bd_order, double dt, int derivative_order = 1)
{
  vector<double> bd_coeffs;
  int n_bd_terms;
  double time_factor;
  if (derivative_order == 1) {
    bd_coeffs = bd_coefficients_first_order_.at(bd_order);
    time_factor = 1 / dt;
    n_bd_terms = bd_order;
  } else if (derivative_order == 2) {
    bd_coeffs = bd_coefficients_second_order_.at(bd_order);
    n_bd_terms = bd_order + 1;
    time_factor = 1 / (dt * dt);
  } else {
    fatal_error("Only first or second order bd derivatives are allowed.");
  }
  T rhs_bd = 0.0;
  for (int i = 0; i < n_bd_terms; i++)
    rhs_bd += bd_coeffs[i + 1] * bd_vector[i];
  rhs_bd *= time_factor;
  return rhs_bd;
}

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_BD_UTILITIES_H
