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
//new
#include "openmc/tallies/tally.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

constexpr double DEFAULT_WEIGHT_CUTOFF {1.0e-38}; // default low weight cutoff

//==============================================================================
// Non-member functions
//==============================================================================

//! Apply weight windows to a particle
//! \param[in] p  Particle to apply weight windows to
void apply_weight_windows(Particle& p);

//! Free memory associated with weight windows
void free_memory_weight_windows();

//==============================================================================
// Global variables
//==============================================================================

class WeightWindows;

namespace variance_reduction {

extern std::unordered_map<int32_t, int32_t> ww_map;
extern vector<unique_ptr<WeightWindows>> weight_windows;
// new in GVR
extern bool global_on;             //!< are GVR are enabled?
extern vector<double> iteration;   //!< iteration during WW generation
extern int32_t tally_idx;          //!< index in tallies vector
extern int source_space;       //!< S vaule in the PS-GVR method
extern int n_statistics;   //!< n_statistics vaule in the PS-GVR method
extern int64_t nps_backup;
// new in asynchronous
extern int max_nps_per_task;       //!< maximum number of particles per task in the asynchronous scheduling function
//new in LVR
extern bool local_on;                     //!< are LVR are enabled?
extern vector<double> target_lower_left;  //!< lower left point of the target 
extern vector<double> target_upper_right; //!< upper right point of the target
extern vector<double> total_weight;       //!< total weight for each mesh bin 
extern vector<double> total_score;        //!< total score for each mesh bin
extern double score_in_target;            //!< total score in the target region
//

} // namespace variance_reduction

//==============================================================================
//! Individual weight window information
//==============================================================================

struct WeightWindow {
  double lower_weight {-1}; // -1 indicates invalid state
  double upper_weight {1};
  double max_lb_ratio {1};
  double survival_weight {0.5};
  double weight_cutoff {DEFAULT_WEIGHT_CUTOFF};
  int max_split {1};

  //! Whether the weight window is in a valid state
  bool is_valid() const { return lower_weight >= 0.0; }

  //! Adjust the weight window by a constant factor
  void scale(double factor)
  {
    lower_weight *= factor;
    upper_weight *= factor;
  }
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

  //! Retrieve the weight window for a particle
  //! \param[in] p  Particle to get weight window for
  WeightWindow get_weight_window(const Particle& p) const;

  // Accessors
  int32_t id() const { return id_; }
  const Mesh& mesh() const { return *model::meshes[mesh_idx_]; }
  
  // new in GVR
  void calculate_WW();

  // new in LVR
  int get_energy_bin(double E);
  int get_energy_size() {return energy_bounds_.size()-1;}

private:
  // Data members
  int32_t id_;                     //!< Unique ID
  ParticleType particle_type_;     //!< Particle type to apply weight windows to
  vector<double> energy_bounds_;   //!< Energy boundaries [eV]
  vector<double> lower_ww_;        //!< Lower weight window bounds
  vector<double> upper_ww_;        //!< Upper weight window bounds
  double survival_ratio_ {3.0};    //!< Survival weight ratio
  double max_lb_ratio_ {1.0}; //!< Maximum lower bound to particle weight ratio
  double weight_cutoff_ {DEFAULT_WEIGHT_CUTOFF}; //!< Weight cutoff
  int max_split_ {10};             //!< Maximum value for particle splitting
  int32_t mesh_idx_;               //!< index in meshes vector
};

} // namespace openmc
#endif // OPENMC_WEIGHT_WINDOWS_H
