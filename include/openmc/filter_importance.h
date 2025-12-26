#ifndef OPENMC_TALLIES_FILTER_IMPORTANCE_H
#define OPENMC_TALLIES_FILTER_IMPORTANCE_H

#include <cstdint>

#include <gsl/gsl-lite.hpp>
#include "openmc/vector.h"
#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Multiplies each tally by the importance corresponding to the mesh index 
//==============================================================================

class ImportanceFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~ImportanceFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override {return "importance";}
  FilterType type() const override { return FilterType::IMPORTANCE; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  const std::vector<double>& importance() const { return importance_; }
  void set_importance(gsl::span<const double> importance);

  virtual int32_t mesh() const {return mesh_;}

  virtual void set_mesh(int32_t mesh);

protected:
  //----------------------------------------------------------------------------
  // Data members

  std::vector<double> importance_;

  int32_t mesh_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_IMPORTANCE_H