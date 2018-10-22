#ifndef OPENMC_TALLIES_FILTER_SPH_HAR_H
#define OPENMC_TALLIES_FILTER_SPH_HAR_H

#include <string>

#include "openmc/error.h"
#include "openmc/tallies/filter.h"


namespace openmc {

//TODO: those integer values are not needed when Fortran interop is removed
enum class SphericalHarmonicsCosine {
  scatter = 1, particle = 2
};

//==============================================================================
//! Gives spherical harmonics expansion moments of a tally score
//==============================================================================

class SphericalHarmonicsFilter : public Filter
{
public:
  std::string type() const override {return "sphericalharmonics";}

  ~SphericalHarmonicsFilter() = default;

  void
  from_xml(pugi::xml_node node) override
  {
    order_ = std::stoi(get_node_value(node, "order"));
    n_bins_ = (order_ + 1) * (order_ + 1);

    if (check_for_node(node, "cosine")) {
      auto cos = get_node_value(node, "cosine", true);
      if (cos == "scatter") {
        cosine_ = SphericalHarmonicsCosine::scatter;
      } else if (cos == "particle") {
        cosine_ = SphericalHarmonicsCosine::particle;
      } else {
        std::stringstream err_msg;
        err_msg << "Unrecognized cosine type, \"" << cos
                << "\" in spherical harmonics filter";
        fatal_error(err_msg);
      }
    }
  }

  void
  get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override
  {
    // Determine cosine term for scatter expansion if necessary
    double wgt[order_ + 1];
    if (cosine_ == SphericalHarmonicsCosine::scatter) {
      calc_pn_c(order_, p->mu, wgt);
    } else {
      for (int i = 0; i < order_ + 1; i++) wgt[i] = 1;
    }

    // Find the Rn,m values
    double rn[n_bins_];
    calc_rn_c(order_, p->last_uvw, rn);

    int j = 0;
    for (int n = 0; n < order_ + 1; n++) {
      // Calculate n-th order spherical harmonics for (u,v,w)
      int num_nm = 2*n + 1;

      // Append the matching (bin,weight) for each moment
      for (int i = 0; i < num_nm; i++) {
        match.weights_.push_back(wgt[n] * rn[j]);
        match.bins_.push_back(++j);
      }
    }
  }

  void
  to_statepoint(hid_t filter_group) const override
  {
    Filter::to_statepoint(filter_group);
    write_dataset(filter_group, "order", order_);
    if (cosine_ == SphericalHarmonicsCosine::scatter) {
      write_dataset(filter_group, "cosine", "scatter");
    } else {
      write_dataset(filter_group, "cosine", "particle");
    }
  }

  std::string
  text_label(int bin) const override
  {
    std::stringstream out;
    for (int n = 0; n < order_ + 1; n++) {
      if (bin <= (n + 1) * (n + 1)) {
        int m = (bin - n*n - 1) - n;
        out << "Spherical harmonic expansion, Y" << n << "," << m;
        return out.str();
      }
    }
  }

  int order_;
  SphericalHarmonicsCosine cosine_ {SphericalHarmonicsCosine::particle};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_SPH_HAR_H
