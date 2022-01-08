#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H

#include <cstdint>
#include <unordered_map>

#include <hdf5.h>
#include <pugixml.hpp>

#include "openmc/constants.h"
#include "openmc/memory.h"
#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

//! Free memory associated with weight windows
void free_memory_weight_windows();

//==============================================================================
// Global variables
//==============================================================================

class WeightWindows;

namespace variance_reduction {

extern std::unordered_map<int32_t, int32_t> ww_map;
extern vector<unique_ptr<WeightWindows>> weight_windows;

} // namespace variance_reduction

//==============================================================================
//! Individual weight window information
//==============================================================================

struct WeightWindow {
  double lower_weight {-1}; // -1 indicates not valid state
  double upper_weight {1};
  double survival_weight {0.5};
  double weight_cutoff {WEIGHT_CUTOFF};
  int max_split {1};

  //! Whether the weight window is in a valid state
  operator bool() const { return lower_weight >= 0.0; }
};

//==============================================================================
//! Weight window settings
//==============================================================================

class WeightWindows {
public:
  // Constructors
  WeightWindows();
  WeightWindows(pugi::xml_node node);

  // Methods

  //! Set the weight window ID
  void set_id(int32_t id = -1);

  // NOTE: This is unused for now but may be used in the future
  //! Write weight window settings to an HDF5 file
  //! \param[in] group  HDF5 group to write to
  void to_hdf5(hid_t group) const;

  //! Return a weight window at the specified index
  //! \param[in] p  Particle to get weight window for
  WeightWindow get_weight_window(const Particle& p) const;

  // Accessors
  int32_t id() const { return id_; }
  const Mesh& mesh() const { return *model::meshes[mesh_idx_]; }

private:
  // Data members
  int32_t id_;                     //!< Unique ID
  ParticleType particle_type_;     //!< Particle type to apply weight windows to
  vector<double> energy_bins_;     //!< Energy bins [eV]
  vector<double> lower_ww_;        //!< Lower weight window bounds
  vector<double> upper_ww_;        //!< Upper weight window bounds
  double survival_ratio_ {3.0};    //!< Survival weight ratio
  double weight_cutoff_ {1.0e-38}; //!< Weight cutoff
  int max_split_ {10};             //!< Maximum value for particle splitting
  int32_t mesh_idx_;               //!< index in meshes vector
};

} // namespace openmc
#endif // OPENMC_WEIGHT_WINDOWS_H
