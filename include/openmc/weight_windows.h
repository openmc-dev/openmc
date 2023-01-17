#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H

#include <cstdint>
#include <unordered_map>

#include <gsl/gsl-lite.hpp>
#include <pugixml.hpp>

#include "openmc/constants.h"
#include "openmc/memory.h"
#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/tallies/tally.h"
#include "openmc/vector.h"

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
  WeightWindows(int32_t id = -1);
  WeightWindows(pugi::xml_node node);
  ~WeightWindows();
  static WeightWindows* create(int32_t id = -1);
  static WeightWindows* from_hdf5(
    hid_t wws_group, const std::string& group_name);
  // Methods

  //! Ready the weight window class for use
  void set_defaults();

  //! Set the weight window ID
  void set_id(int32_t id = -1);

  void set_energy_bounds(gsl::span<const double> bounds);

  void set_mesh(const std::unique_ptr<Mesh>& mesh);

  void set_mesh(const Mesh* mesh);

  void set_mesh(int32_t mesh_idx);

  // ! Update weight window boundaries using tally results
  void update_weight_windows(const std::unique_ptr<Tally>& tally,
                             const std::string& score,
                             const std::string& value,
                             const std::string& method);

  void export_to_hdf5(const std::string& filename = "weight_windows.h5") const;

  // NOTE: This is unused for now but may be used in the future
  //! Write weight window settings to an HDF5 file
  //! \param[in] group  HDF5 group to write to
  void to_hdf5(hid_t group) const;

  //! Retrieve the weight window for a particle
  //! \param[in] p  Particle to get weight window for
  WeightWindow get_weight_window(const Particle& p) const;

  double bounds_size() const;

  const vector<double>& energy_bounds() const { return energy_bounds_; }

  void set_weight_windows(
    gsl::span<const double> lower_bounds, gsl::span<const double> upper_bounds);

  void set_weight_windows(
    gsl::span<const double> lower_bounds, double bound_ratio);

  void set_particle_type(ParticleType p_type);

  // Accessors
  int32_t id() const { return id_; }
  int32_t& id() { return id_; }

  vector<double>& energy_bounds() { return energy_bounds_; }

  const std::unique_ptr<Mesh>& mesh() const { return model::meshes[mesh_idx_]; }

  const xt::xarray<double>& lower_bounds() const { return lower_ww_; }
  xt::xarray<double>& lower_bounds() { return lower_ww_; }

  const xt::xarray<double>& upper_bounds() const { return upper_ww_; }
  xt::xarray<double>& upper_bounds() { return upper_ww_; }

  ParticleType particle_type() const { return particle_type_; }

private:
  // Data members
  int32_t id_;       //!< Unique ID
  gsl::index index_; //!< Index into weight windows vector
  ParticleType particle_type_ {
    ParticleType::neutron};      //!< Particle type to apply weight windows to
  vector<double> energy_bounds_; //!< Energy boundaries [eV]
  xt::xarray<double> lower_ww_;  //!< Lower weight window bounds
  xt::xarray<double> upper_ww_;  //!< Upper weight window bounds
  double survival_ratio_ {3.0};  //!< Survival weight ratio
  double max_lb_ratio_ {1.0}; //!< Maximum lower bound to particle weight ratio
  double weight_cutoff_ {DEFAULT_WEIGHT_CUTOFF}; //!< Weight cutoff
  int max_split_ {10}; //!< Maximum value for particle splitting
  int32_t mesh_idx_;   //!< Index in meshes vector
};

} // namespace openmc
#endif // OPENMC_WEIGHT_WINDOWS_H
