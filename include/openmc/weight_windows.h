#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H

#include <array>
#include <cstdint>
#include <unordered_map>

#include <hdf5.h>
#include <pugixml.hpp>

#include "openmc/constants.h"
#include "openmc/memory.h"
#include "openmc/mesh.h"
#include "openmc/particle_type.h"
#include "openmc/span.h"
#include "openmc/tallies/tally.h"
#include "openmc/tensor.h"
#include "openmc/vector.h"

namespace openmc {

enum class WeightWindowUpdateMethod { MAGIC, FW_CADIS };

//==============================================================================
// Constants
//==============================================================================

constexpr double DEFAULT_WEIGHT_CUTOFF {1.0e-38}; // default low weight cutoff

//==============================================================================
// Global variables
//==============================================================================

class WeightWindows;
class WeightWindowsGenerator;
class SourceBias;

namespace variance_reduction {

extern std::unordered_map<int32_t, int32_t> ww_map;
extern vector<unique_ptr<WeightWindows>> weight_windows;
extern vector<unique_ptr<WeightWindowsGenerator>> weight_windows_generators;
extern vector<unique_ptr<SourceBias>> source_biases;

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
  int max_split {10};

  //! Whether the weight window is in a valid state. A non-positive lower
  //! bound indicates that no weight window information exists at this
  //! location (generators mark such cells with -1, and a lower bound of zero
  //! conventionally turns the weight window game off in a cell, as in MCNP
  //! wwinp files), in which case no weight window game is played.
  bool is_valid() const { return lower_weight > 0.0; }

  //! Adjust the weight window by a constant factor
  void scale(double factor)
  {
    lower_weight *= factor;
    upper_weight *= factor;
    survival_weight *= factor;
  }
};

//==============================================================================
//! Weight window settings
//==============================================================================

class WeightWindows {
public:
  //----------------------------------------------------------------------------
  // Constructors
  WeightWindows(int32_t id = -1);
  WeightWindows(pugi::xml_node node);
  ~WeightWindows();
  static WeightWindows* create(int32_t id = -1);
  static WeightWindows* from_hdf5(
    hid_t wws_group, const std::string& group_name);

  //----------------------------------------------------------------------------
  // Methods
private:
  template<class T>
  void check_bounds(const T& lower, const T& upper) const;

  template<class T>
  void check_bounds(const T& lower) const;

  void check_tally_update_compatibility(const Tally* tally);

public:
  //! Set the weight window ID
  void set_id(int32_t id = -1);

  void set_energy_bounds(span<const double> bounds);

  void set_mesh(const std::unique_ptr<Mesh>& mesh);

  void set_mesh(const Mesh* mesh);

  void set_mesh(int32_t mesh_idx);

  //! Ready the weight window class for use
  void set_defaults();

  //! Ensure the weight window lower bounds are properly allocated
  void allocate_ww_bounds();

  //! Update weight window boundaries using tally results
  //! \param[in] tally Pointer to the tally whose results will be used to
  //! update weight windows \param[in] value String representing the type of
  //! value to use for weight window generation (one of "mean" or "rel_err")
  //! \param[in] threshold Relative error threshold. Results over this
  //! threshold will be ignored \param[in] ratio Ratio of upper to lower
  //! weight window bounds
  void update_weights(const Tally* tally, const std::string& value = "mean",
    double threshold = 1.0, double ratio = 5.0,
    WeightWindowUpdateMethod method = WeightWindowUpdateMethod::MAGIC);

  // NOTE: This is unused for now but may be used in the future
  //! Write weight window settings to an HDF5 file
  //! \param[in] group  HDF5 group to write to
  void to_hdf5(hid_t group) const;

  //! Retrieve the weight window for a particle
  //! \param[in] p  Particle to get weight window for
  std::pair<bool, WeightWindow> get_weight_window(const Particle& p) const;

  std::array<int, 2> bounds_size() const;

  const vector<double>& energy_bounds() const { return energy_bounds_; }

  void set_bounds(const tensor::Tensor<double>& lower_ww_bounds,
    const tensor::Tensor<double>& upper_bounds);

  void set_bounds(const tensor::Tensor<double>& lower_bounds, double ratio);

  void set_bounds(
    span<const double> lower_bounds, span<const double> upper_bounds);

  void set_bounds(span<const double> lower_bounds, double ratio);

  void set_particle_type(ParticleType p_type);

  double survival_ratio() const { return survival_ratio_; }

  double& survival_ratio() { return survival_ratio_; }

  double max_lower_bound_ratio() const { return max_lb_ratio_; }

  double& max_lower_bound_ratio() { return max_lb_ratio_; }

  int max_split() const { return max_split_; }

  int& max_split() { return max_split_; }

  double weight_cutoff() const { return weight_cutoff_; }

  double& weight_cutoff() { return weight_cutoff_; }

  //----------------------------------------------------------------------------
  // Accessors
  int32_t id() const { return id_; }
  int32_t& id() { return id_; }

  int32_t index() const { return index_; }

  vector<double>& energy_bounds() { return energy_bounds_; }

  const std::unique_ptr<Mesh>& mesh() const { return model::meshes[mesh_idx_]; }

  const tensor::Tensor<double>& lower_ww_bounds() const { return lower_ww_; }
  tensor::Tensor<double>& lower_ww_bounds() { return lower_ww_; }

  const tensor::Tensor<double>& upper_ww_bounds() const { return upper_ww_; }
  tensor::Tensor<double>& upper_ww_bounds() { return upper_ww_; }

  ParticleType particle_type() const { return particle_type_; }

private:
  //----------------------------------------------------------------------------
  // Data members
  int32_t id_;                   //!< Unique ID
  int64_t index_;                //!< Index into weight windows vector
  ParticleType particle_type_;   //!< Particle type to apply weight windows to
  vector<double> energy_bounds_; //!< Energy boundaries [eV]
  tensor::Tensor<double> lower_ww_; //!< Lower weight window bounds (shape:
                                    //!< energy_bins, mesh_bins (k, j, i))
  tensor::Tensor<double>
    upper_ww_; //!< Upper weight window bounds (shape: energy_bins, mesh_bins)
  double survival_ratio_ {3.0}; //!< Survival weight ratio
  double max_lb_ratio_ {1.0}; //!< Maximum lower bound to particle weight ratio
  double weight_cutoff_ {DEFAULT_WEIGHT_CUTOFF}; //!< Weight cutoff
  int max_split_ {10};    //!< Maximum value for particle splitting
  int32_t mesh_idx_ {-1}; //!< Index in meshes vector
};

class SourceBias {
public:
  //----------------------------------------------------------------------------
  // Constructors

  //! \param[in] spatial_mesh_idx Index into model::meshes for the spatial
  //!   mesh over which the source is biased
  //! \param[in] angle_mesh_idx Index into model::meshes for the angular mesh
  //!   used to bias emission direction, or C_NONE for isotropic emission
  //! \param[in] energy_bounds Energy group boundaries [eV]
  //! \param[in] ww_id ID of the WeightWindows object this source bias is
  //!   generated alongside (used to name this object's HDF5 group)
  //! \param[in] forward_source_mesh Path to unbiased source strength mesh file
  SourceBias(int32_t spatial_mesh_idx, int32_t angle_mesh_idx,
    vector<double> energy_bounds, int32_t ww_id,
    std::string forward_source_mesh);

  //----------------------------------------------------------------------------
  // Methods

  //! \param[in] tally Tally with a "flux" score, a mesh filter matching
  //!   spatial_mesh_idx_, an energy filter matching energy_bounds_ (if
  //!   energy_bounds_ is non-empty), and a mesh-angular filter matching
  //!   angle_mesh_idx_ (if angle_mesh_idx_ != C_NONE)
  void update(const Tally* tally);

  //! Write SourceBias data to an HDF5 group
  void to_hdf5(hid_t group) const;

  //----------------------------------------------------------------------------
  // Accessors

  int32_t spatial_mesh_idx() const { return spatial_mesh_idx_; }
  int32_t angle_mesh_idx() const { return angle_mesh_idx_; }
  std::string forward_source_mesh() const { return forward_source_mesh_; }
  //! Estimate the unbiased external source's strength as a function of
  //! spatial mesh element, angular mesh element, and energy group, assuming
  //! isotropic emission when no angular mesh is present.
  //!
  //! If the source strengths are already written to a file
  //! such as "forward_source_mesh.h5", that data is used directly in place of
  //! sampling model::external_sources. Its spatial mesh, angular mesh, and
  //! energy group structure are checked against spatial_mesh_, angle_mesh_,
  //! and energy_bounds_ before use.
  //!
  //! Otherwise, group probabilities are read exactly from each source's
  //! Discrete energy spectrum, while the spatial and angular dimensions are
  //! estimated by binning samples into (spatial_mesh_, angle_mesh_).
  //!
  //! \param[in] n_samples_per_source Number of trial (position, direction)
  //!   pairs to draw per external source. Unused if forward source mesh file
  //!   is present.
  void compute_unbiased_strength(int64_t n_samples_per_source = 1000000);

  const tensor::Tensor<double>& unbiased_strength() const
  {
    return unbiased_strength_;
  }

private:
  //! Read a precomputed unbiased source strength distribution from a
  //! forward source mesh file
  void load_forward_source_mesh(const std::string& path);

  //----------------------------------------------------------------------------
  // Data members

  int32_t spatial_mesh_idx_;        //!< Index into model::meshes
  int32_t angle_mesh_idx_;          //!< Index into model::meshes, or C_NONE
  vector<double> energy_bounds_;    //!< Energy group boundaries [eV]
  int32_t ww_id_;                   //!< ID of the associated WeightWindows
  std::string forward_source_mesh_; //!< Path to unbiased source strengths file

  //! Un-normalized mean flux. Shape: (spatial_bins, angle_bins,
  //! energy_bins). angle_bins == 1 when angle_mesh_idx_ == C_NONE.
  tensor::Tensor<double> flux_;

  //! Estimated relative strength of the unbiased external source as a function
  //! of (spatial mesh element, angular mesh element, energy group), assuming
  //! isotropic emission when angle_mesh_idx_ == C_NONE.
  tensor::Tensor<double> unbiased_strength_;
};

class WeightWindowsGenerator {
public:
  // Constructors
  WeightWindowsGenerator(pugi::xml_node node);

  // Methods
  void update() const;

  //! Create the tally used for weight window generation
  void create_tally();

  // Data members
  int32_t tally_idx_; //!< Index of the tally used to update the weight windows
  int32_t sb_tally_idx_ {C_NONE}; //!< Index of tally used to update source bias
  int32_t sb_idx_ {C_NONE}; //!< Index into variance_reduction::source_biases
  int32_t ww_idx_; //!< Index of the weight windows object being generated
  WeightWindowUpdateMethod method_; //!< Method used to update weight window.
  int32_t max_realizations_;        //!< Maximum number of tally realizations
  int32_t update_interval_;         //!< Determines how often updates occur
  bool on_the_fly_; //!< Whether or not to keep tally results between batches or
                    //!< realizations

  // MAGIC update parameters
  std::string tally_value_ {
    "mean"};               //<! Tally value to use (one of {"mean", "rel_err"})
  double threshold_ {1.0}; //<! Relative error threshold for values used to
                           // update weight windows
  double ratio_ {5.0};     //<! ratio of lower to upper weight window bounds

  // Local FW-CADIS target tallies
  std::vector<size_t> targets_;
  // FW-CADIS source biasing
  bool source_biasing_;
  int32_t angle_mesh_idx_ {C_NONE}; //<! Index in mesh map of optional angular
                                    // quadrature for source biasing
  std::string forward_source_mesh_ {
    "forward_source_mesh.h5"}; //<! Path to unbiased source strengths mesh file
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Apply weight windows to a particle
//! \param[in] p  Particle to apply weight windows to
void apply_weight_windows(Particle& p);

//! Apply weight window to a particle
//! \param[in] p  Particle to apply weight window to
//! \param[in] weight_window WeightWindow to apply
void apply_weight_window(Particle& p, WeightWindow weight_window);

//! Free memory associated with weight windows
void free_memory_weight_windows();

//! Search weight window that apply to a particle
//! \param[in]  p  Particle to search weight window for
std::pair<bool, WeightWindow> search_weight_window(const Particle& p);

//! Finalize variance reduction objects after all inputs have been read
void finalize_variance_reduction();

} // namespace openmc
#endif // OPENMC_WEIGHT_WINDOWS_H
