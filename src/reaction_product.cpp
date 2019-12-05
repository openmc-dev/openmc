#include "openmc/reaction_product.h"

#include <memory> // for unique_ptr
#include <string> // for string

#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/secondary_correlated.h"
#include "openmc/secondary_kalbach.h"
#include "openmc/secondary_nbody.h"
#include "openmc/secondary_uncorrelated.h"

namespace openmc {

//==============================================================================
// ReactionProduct implementation
//==============================================================================

ReactionProduct::ReactionProduct(hid_t group)
{
  // Read particle type
  std::string temp;
  read_attribute(group, "particle", temp);
  if (temp == "neutron") {
    particle_ = Particle::Type::neutron;
  } else if (temp == "photon") {
    particle_ = Particle::Type::photon;
  }

  // Read emission mode and decay rate
  read_attribute(group, "emission_mode", temp);
  if (temp == "prompt") {
    emission_mode_ = EmissionMode::prompt;
  } else if (temp == "delayed") {
    emission_mode_ = EmissionMode::delayed;
  } else if (temp == "total") {
    emission_mode_ = EmissionMode::total;
  }

  // Read decay rate for delayed emission
  if (emission_mode_ == EmissionMode::delayed)
    read_attribute(group, "decay_rate", decay_rate_);

  // Read secondary particle yield
  yield_ = read_function(group, "yield");

  int n;
  read_attribute(group, "n_distribution", n);

  for (int i = 0; i < n; ++i) {
    std::string s {"distribution_"};
    s.append(std::to_string(i));
    hid_t dgroup = open_group(group, s.c_str());

    // Read applicability
    if (n > 1) {
      hid_t app = open_dataset(dgroup, "applicability");
      applicability_.emplace_back(app);
      close_dataset(app);
    }

    // Determine distribution type and read data
    read_attribute(dgroup, "type", temp);
    if (temp == "uncorrelated") {
      distribution_.push_back(std::make_unique<UncorrelatedAngleEnergy>(dgroup));
    } else if (temp == "correlated") {
      distribution_.push_back(std::make_unique<CorrelatedAngleEnergy>(dgroup));
    } else if (temp == "nbody") {
      distribution_.push_back(std::make_unique<NBodyPhaseSpace>(dgroup));
    } else if (temp == "kalbach-mann") {
      distribution_.push_back(std::make_unique<KalbachMann>(dgroup));
    }

    close_group(dgroup);
  }
}

void ReactionProduct::sample(double E_in, double& E_out, double& mu,
  uint64_t* seed) const
{
  auto n = applicability_.size();
  if (n > 1) {
    double prob = 0.0;
    double c = prn(seed);
    for (int i = 0; i < n; ++i) {
      // Determine probability that i-th energy distribution is sampled
      prob += applicability_[i](E_in);

      // If i-th distribution is sampled, sample energy from the distribution
      if (c <= prob) {
        distribution_[i]->sample(E_in, E_out, mu, seed);
        break;
      }
    }
  } else {
    // If only one distribution is present, go ahead and sample it
    distribution_[0]->sample(E_in, E_out, mu, seed);
  }
}

}
