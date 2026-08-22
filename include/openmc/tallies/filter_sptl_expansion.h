#ifndef OPENMC_TALLIES_FILTER_SPTL_EXPANSION_H
#define OPENMC_TALLIES_FILTER_SPTL_EXPANSION_H

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Abstract filter for functional expansions of the particle's normalized
//! position along a Cartesian axis
//==============================================================================

class SpatialExpansionFilter : public Filter {
public:
  enum class Axis { x, y, z };

  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~SpatialExpansionFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  void from_xml(pugi::xml_node node) override;

  void to_statepoint(hid_t filter_group) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int order() const { return order_; }

  //! Set the expansion order and the corresponding number of bins
  virtual void set_order(int order) = 0;

  Axis axis() const { return axis_; }
  void set_axis(Axis axis);

  double min() const { return min_; }
  double max() const { return max_; }
  void set_minmax(double min, double max);

protected:
  //----------------------------------------------------------------------------
  // Methods

  //! Get the particle's coordinate along the expansion axis
  double position(const Particle& p) const;

  //! Get the name of the expansion axis ("x", "y", or "z")
  const char* axis_label() const;

  //----------------------------------------------------------------------------
  // Data members

  int order_;

  //! The Cartesian coordinate axis that the expansion is applied to.
  Axis axis_;

  //! The minimum coordinate along the reference axis that the expansion covers.
  double min_;

  //! The maximum coordinate along the reference axis that the expansion covers.
  double max_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_SPTL_EXPANSION_H
