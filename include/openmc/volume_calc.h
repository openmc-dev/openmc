#ifndef OPENMC_VOLUME_CALC_H
#define OPENMC_VOLUME_CALC_H

#include <algorithm> // for find
#include <cstdint>
#include <string>
#include <vector>

#include "openmc/array.h"
#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/tallies/trigger.h"
#include "openmc/vector.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"
#include <gsl/gsl-lite.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace openmc {

//==============================================================================
// Volume calculation class
//==============================================================================

class VolumeCalculation {

public:
  // Aliases, types
  struct Result {
    array<double, 2> volume;    //!< Mean/standard deviation of volume
    vector<int> nuclides;       //!< Index of nuclides
    vector<double> atoms;       //!< Number of atoms for each nuclide
    vector<double> uncertainty; //!< Uncertainty on number of atoms
    int iterations; //!< Number of iterations needed to obtain the results
  };                // Results for a single domain

  // Constructors
  VolumeCalculation(pugi::xml_node node);

  VolumeCalculation() = default;

  // Methods

  //! \brief Stochastically determine the volume of a set of domains along with
  //! the
  //!   average number densities of nuclides within the domain
  //
  //! \return Vector of results for each user-specified domain
  vector<Result> execute() const;

  //! \brief Write volume calculation results to HDF5 file
  //
  //! \param[in] filename  Path to HDF5 file to write
  //! \param[in] results   Vector of results for each domain
  void to_hdf5(
    const std::string& filename, const vector<Result>& results) const;

  // Tally filter and map types
  enum class TallyDomain { UNIVERSE, MATERIAL, CELL };

  // Data members
  TallyDomain domain_type_; //!< Type of domain (cell, material, etc.)
  size_t n_samples_;        //!< Number of samples to use
  double threshold_ {-1.0}; //!< Error threshold for domain volumes
  TriggerMetric trigger_type_ {
    TriggerMetric::not_active}; //!< Trigger metric for the volume calculation
  Position lower_left_;         //!< Lower-left position of bounding box
  Position upper_right_;        //!< Upper-right position of bounding box
  vector<int> domain_ids_;      //!< IDs of domains to find volumes of

private:
  //! \brief Check whether a material has already been hit for a given domain.
  //! If not, add new entries to the vectors
  //
  //! \param[in] i_material Index in global materials vector
  //! \param[in,out] indices Vector of material indices
  //! \param[in,out] hits Number of hits corresponding to each material
  void check_hit(
    int i_material, vector<uint64_t>& indices, vector<uint64_t>& hits) const;
};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
extern vector<VolumeCalculation> volume_calcs;
}

//==============================================================================
// Non-member functions
//==============================================================================

//! Reduce vector of indices and hits from each thread to a single copy
//
//! \param[in] local_indices Indices specific to each thread
//! \param[in] local_hits Hit count specific to each thread
//! \param[out] indices Reduced vector of indices
//! \param[out] hits Reduced vector of hits
template<typename T, typename T2>
void reduce_indices_hits(const vector<T>& local_indices,
  const vector<T2>& local_hits, vector<T>& indices, vector<T2>& hits)
{
  const int n_threads = num_threads();

#pragma omp for ordered schedule(static, 1)
  for (int i = 0; i < n_threads; ++i) {
#pragma omp ordered
    for (int j = 0; j < local_indices.size(); ++j) {
      // Check if this material has been added to the master list and if
      // so, accumulate the number of hits
      auto it = std::find(indices.begin(), indices.end(), local_indices[j]);
      if (it == indices.end()) {
        indices.push_back(local_indices[j]);
        hits.push_back(local_hits[j]);
      } else {
        hits[it - indices.begin()] += local_hits[j];
      }
    }
  }
}

void free_memory_volume();

} // namespace openmc

#endif // OPENMC_VOLUME_CALC_H
