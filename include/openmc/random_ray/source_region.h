#ifndef OPENMC_RANDOM_RAY_SOURCE_REGION_H
#define OPENMC_RANDOM_RAY_SOURCE_REGION_H

#include <vector>

#include "openmc/openmp_interface.h"
#include "openmc/position.h"

namespace openmc {

namespace random_ray {

// Scalars
extern int64_t n_source_elements;
extern int64_t n_source_regions; // number of source regions (a.k.a. cells)

// 1D arrays representing values for each OpenMC "Cell"
extern std::vector<int64_t> source_region_offsets;

// Cell-wise Data
extern std::vector<OpenMPMutex> lock;
extern std::vector<int> material;
extern std::vector<int> position_recorded;
extern std::vector<Position> position;
extern std::vector<float> scalar_flux_new;
extern std::vector<float> scalar_flux_old;
extern std::vector<float> source;
extern std::vector<double> volume;
extern std::vector<double> volume_t;
extern std::vector<int> was_hit;

} // namespace random_ray

void initialize_source_regions();

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_SOURCE_REGION_H
