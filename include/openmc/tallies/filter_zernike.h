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

  void get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  virtual void calc_n_bins() {n_bins_ = ((order_+1) * (order_+2)) / 2;}

  int order_;
  double x_, y_, r_;
};

//==============================================================================
//! Gives even order radial Zernike polynomial moments of a particle's position
//==============================================================================

class ZernikeRadialFilter : public ZernikeFilter
{
public:
  std::string type() const override {return "zernikeradial";}

  void get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;

  void calc_n_bins() override {n_bins_ = order_ / 2 + 1;}
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ZERNIKE_H
