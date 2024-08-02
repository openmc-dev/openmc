#ifndef OPENMC_TALLIES_FILTER_TIMED_MESH_H
#define OPENMC_TALLIES_FILTER_TIMED_MESH_H

#include <cstdint>

#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Indexes the time-location of particle events to a time gridn and a regular 
//  mesh. For tracklength tallies, it will produce multiple valid bins and the
//  bin weight will correspond to the fraction of the track length that lies in
//! that bin.
//==============================================================================

class TimedMeshFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~TimedMeshFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "timedmesh"; }
  FilterType type() const override { return FilterType::TIMED_MESH; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int32_t mesh() const { return mesh_; }

  void set_mesh(int32_t mesh);

  void set_translation(const Position& translation);

  void set_translation(const double translation[3]);

  const Position& translation() const { return translation_; }

  bool translated() const { return translated_; }
  
  const vector<double>& time_grid() const { return time_grid_; }
  
  void set_time_grid(gsl::span<const double> time_grid);
  
  void reset_bins();

protected:
  //----------------------------------------------------------------------------
  // Data members

  int32_t mesh_;            //!< Index of the mesh
  int mesh_n_bins_;
  bool translated_ {false}; //!< Whether or not the filter is translated
  Position translation_ {0.0, 0.0, 0.0}; //!< Filter translation
  vector<double> time_grid_ {0.0, INFTY};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_TIMED_MESH_H
