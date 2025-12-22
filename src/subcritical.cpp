#include "openmc/subcritical.h"
#include "openmc/simulation.h"
#include <cmath>
#include <utility>

namespace openmc{
    
void convert_to_subcritical_k(double& k, double& k_std) {
    double k_sub = k / (k + 1);
    double k_sub_std = k_sub * sqrt( pow(k_std / k, 2) + pow(k_std / (k + 1), 2) );
    k = k_sub;
    k_std = k_sub_std;
}

double convert_to_subcritical_k(double k) {
    // Used to convert keff estimators in fixed source mode, which are multiplicities, specifically = M-1
    // into the corresponding estimators for subcritical k
    double k_sub = (k/simulation::total_weight) / (k/simulation::total_weight + 1) * simulation::total_weight;
    return k_sub;
}

} // namespace openmc