#ifndef OPENMC_TALLY_FILTER_POLAR_H
#define OPENMC_TALLY_FILTER_POLAR_H

#include <sstream>
#include <cmath>
#include <vector>

#include "openmc/error.h"
#include "openmc/search.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

class PolarFilter : public TallyFilter
{
public:
  virtual std::string type() const override {return "polar";}

  virtual ~PolarFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    auto bins = get_node_array<double>(node, "bins");

    if (bins.size() > 1) {
      bins_ = bins;

    } else {
      // Allow a user to input a lone number which will mean that you subdivide
      // [0,pi] evenly with the input being the number of bins

      int n_angle = bins[0];

      if (n_angle <= 1) fatal_error("Number of bins for polar filter must "
                                    "be greater than 1.");

      double d_angle = PI / n_angle;
      bins_.resize(n_angle + 1);
      for (int i = 0; i < n_angle; i++) bins_[i] = i * d_angle;
      bins_[n_angle] = PI;
    }

    n_bins_ = bins_.size() - 1;
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    double theta;
    if (estimator == ESTIMATOR_TRACKLENGTH) {
      theta = std::acos(p->coord[0].uvw[2]);
    } else {
      theta = std::acos(p->last_uvw[2]);
    }

    if (theta >= bins_[0] && theta <= bins_.back()) {
      auto bin = lower_bound_index(bins_.begin(), bins_.end(), theta) + 1;
      match.bins.push_back(bin);
      match.weights.push_back(1);
    }
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    TallyFilter::to_statepoint(filter_group);
    write_dataset(filter_group, "bins", bins_);
  }

  virtual std::string
  text_label(int bin) const override
  {
    std::stringstream out;
    out << "Polar Angle [" << bins_[bin-1] << ", " << bins_[bin] << ")";
    return out.str();
  }

protected:
  std::vector<double> bins_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_POLAR_H
