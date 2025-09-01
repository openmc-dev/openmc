#include "openmc/reaction_product.h"

#include <string> // for string

#include <fmt/core.h>

#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/random_lcg.h"
#include "openmc/secondary_correlated.h"
#include "openmc/secondary_kalbach.h"
#include "openmc/secondary_nbody.h"
#include "openmc/secondary_thermal.h"
#include "openmc/secondary_uncorrelated.h"
#include "openmc/tallies/tally_scoring.h"

namespace openmc {

//==============================================================================
// ReactionProduct implementation
//==============================================================================

ReactionProduct::ReactionProduct(hid_t group)
{
  // Read particle type
  std::string temp;
  read_attribute(group, "particle", temp);
  particle_ = str_to_particle_type(temp);

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
  if (emission_mode_ == EmissionMode::delayed) {
    if (attribute_exists(group, "decay_rate")) {
      read_attribute(group, "decay_rate", decay_rate_);
    } else if (particle_ == ParticleType::neutron) {
      warning(fmt::format("Decay rate doesn't exist for delayed neutron "
                          "emission ({}).",
        object_name(group)));
    }
  }

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
      distribution_.push_back(make_unique<UncorrelatedAngleEnergy>(dgroup));
    } else if (temp == "correlated") {
      distribution_.push_back(make_unique<CorrelatedAngleEnergy>(dgroup));
    } else if (temp == "nbody") {
      distribution_.push_back(make_unique<NBodyPhaseSpace>(dgroup));
    } else if (temp == "kalbach-mann") {
      distribution_.push_back(make_unique<KalbachMann>(dgroup));
    }

    close_group(dgroup);
  }
}

ReactionProduct::ReactionProduct(const ChainNuclide::Product& product)
{
  particle_ = ParticleType::photon;
  emission_mode_ = EmissionMode::delayed;

  // Get chain nuclide object for radionuclide
  parent_nuclide_ = data::chain_nuclide_map.at(product.name);
  const auto& chain_nuc = data::chain_nuclides[parent_nuclide_].get();

  // Determine decay constant in [s^-1]
  decay_rate_ = chain_nuc->decay_constant();

  // Determine number of photons per decay and set yield
  double photon_per_sec = chain_nuc->photon_energy()->integral();
  double photon_per_decay = photon_per_sec / decay_rate_;
  vector<double> coef = {product.branching_ratio * photon_per_decay};
  yield_ = make_unique<Polynomial>(coef);

  // Set decay photon angle-energy distribution
  distribution_.push_back(
    make_unique<DecayPhotonAngleEnergy>(chain_nuc->photon_energy()));
}

void ReactionProduct::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
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

void ReactionProduct::get_pdf(int i_tally, double E_in, double& E_out,
  uint64_t* seed, double mu_cm) const
{


  int distribution_index;
  auto n = applicability_.size();
  if (n > 1) {
    double prob = 0.0;
    double c = prn(seed);
    for (int i = 0; i < n; ++i) {
      // Determine probability that i-th energy distribution is sampled
      prob += applicability_[i](E_in);

      // If i-th distribution is sampled, sample energy from the distribution
      if (c <= prob) {
        // distribution_[i]->sample(E_in, E_out, mu, seed);
        distribution_index = i;
        break;
      }
    }
  } else {
    // If only one distribution is present, go ahead and sample it
    // distribution_[0]->sample(E_in, E_out, mu, seed);
    distribution_index = 0;
  }
  // now extract pdf

  AngleEnergy* angleEnergyPtr = distribution_[distribution_index].get();

}

} // namespace openmc
