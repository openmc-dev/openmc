#ifndef OPENMC_TALLIES_FILTER_ZERNIKE_H
#define OPENMC_TALLIES_FILTER_ZERNIKE_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Gives Zernike polynomial moments of a particle's position
//==============================================================================

class ZernikeFilter : public Filter
{
public:
  std::string type() const override {return "zernike";}

  ~ZernikeFilter() = default;

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  int order() const {return order_;}

  virtual void set_order(int order);

  //! Cartesian x coordinate for the origin of this expansion.
  double x_;

  //! Cartesian y coordinate for the origin of this expansion.
  double y_;

  //! Maximum radius from the origin covered by this expansion.
  double r_;

protected:
  int order_;
};

//==============================================================================
//! Gives even order radial Zernike polynomial moments of a particle's position
//==============================================================================

class ZernikeRadialFilter : public ZernikeFilter
{
public:
  std::string type() const override {return "zernikeradial";}

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;

  void set_order(int order) override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ZERNIKE_H
