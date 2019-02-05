#ifndef OPENMC_TALLIES_TALLY_H
#define OPENMC_TALLIES_TALLY_H

#include "openmc/constants.h"
#include "openmc/tallies/trigger.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <memory> // for unique_ptr
#include <string>
#include <vector>

namespace openmc {

//==============================================================================
//! A user-specified flux-weighted (or current) measurement.
//==============================================================================

class Tally {
public:
  Tally() {}

  //----------------------------------------------------------------------------
  // Methods for getting and setting filter/stride data.

  const std::vector<int32_t>& filters() const {return filters_;}

  int32_t filters(int i) const {return filters_[i];}

  void set_filters(const int32_t filter_indices[], int n);

  int32_t strides(int i) const {return strides_[i];}

  int32_t n_filter_bins() const {return n_filter_bins_;}

  //----------------------------------------------------------------------------
  // Other methods.

  void init_triggers(pugi::xml_node node, int i_tally);

  //----------------------------------------------------------------------------
  // Major public data members.

  int id_; //!< user-defined identifier

  int type_ {TALLY_VOLUME}; //!< volume, surface current

  //! Event type that contributes to this tally
  int estimator_ {ESTIMATOR_TRACKLENGTH};

  //! Whether this tally is currently being updated
  bool active_ {false};

  //! Index of each nuclide to be tallied.  -1 indicates total material.
  std::vector<int> nuclides_;

  //! True if this tally has a bin for every nuclide in the problem
  bool all_nuclides_ {false};

  //----------------------------------------------------------------------------
  // Miscellaneous public members.

  // We need to have quick access to some filters.  The following gives indices
  // for various filters that could be in the tally or C_NONE if they are not
  // present.
  int energyout_filter_ {C_NONE};
  int delayedgroup_filter_ {C_NONE};

  std::vector<Trigger> triggers_;

  int deriv_ {C_NONE}; //!< Index of a TallyDerivative object for diff tallies.

private:
  //----------------------------------------------------------------------------
  // Private data.

  std::vector<int32_t> filters_; //!< Filter indices in global filters array

  //! Index strides assigned to each filter to support 1D indexing.
  std::vector<int32_t> strides_;

  int32_t n_filter_bins_ {0};
};

//==============================================================================
// Global variable declarations
//==============================================================================

extern "C" double total_weight;

namespace model {
  extern std::vector<std::unique_ptr<Tally>> tallies;

  extern std::vector<int> active_tallies;
  extern std::vector<int> active_analog_tallies;
  extern std::vector<int> active_tracklength_tallies;
  extern std::vector<int> active_collision_tallies;
  extern std::vector<int> active_meshsurf_tallies;
  extern std::vector<int> active_surface_tallies;
}

// Threadprivate variables
extern "C" double global_tally_absorption;
extern "C" double global_tally_collision;
extern "C" double global_tally_tracklength;
extern "C" double global_tally_leakage;
#pragma omp threadprivate(global_tally_absorption, global_tally_collision, \
  global_tally_tracklength, global_tally_leakage)

//==============================================================================
// Non-member functions
//==============================================================================

// Alias for the type returned by xt::adapt(...). N is the dimension of the
// multidimensional array
template <std::size_t N>
using adaptor_type = xt::xtensor_adaptor<xt::xbuffer_adaptor<double*&, xt::no_ownership>, N>;

//! Get the global tallies as a multidimensional array
//! \return Global tallies array
adaptor_type<2> global_tallies();

//! Get tally results as a multidimensional array
//! \param idx Index in tallies array
//! \return Tally results array
adaptor_type<3> tally_results(int idx);

#ifdef OPENMC_MPI
//! Collect all tally results onto master process
extern "C" void reduce_tally_results();
#endif

extern "C" void free_memory_tally_c();

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_H
