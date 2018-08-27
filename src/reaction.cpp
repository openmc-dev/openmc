#include "openmc/reaction.h"

#include <string>
#include <utility> // for move

#include "openmc/hdf5_interface.h"
#include "openmc/endf.h"
#include "openmc/random_lcg.h"
#include "openmc/secondary_uncorrelated.h"

namespace openmc {

Reaction::Reaction(hid_t group, const std::vector<int>& temperatures)
{
  read_attribute(group, "Q_value", q_value_);
  read_attribute(group, "mt", mt_);
  int cm;
  read_attribute(group, "center_of_mass", cm);
  scatter_in_cm_ = (cm == 1);

  // Read cross section and threshold_idx data
  for (auto t : temperatures) {
    // Get group corresponding to temperature
    std::string temp_str {std::to_string(t) + "K"};
    hid_t temp_group = open_group(group, temp_str.c_str());
    hid_t dset = open_dataset(temp_group, "xs");

    // Get threshold index
    TemperatureXS xs;
    read_attribute(dset, "threshold_idx", xs.threshold);

    // Read cross section values
    read_dataset(dset, xs.value);
    close_dataset(dset);
    close_group(temp_group);

    // create new entry in xs vector
    xs_.push_back(std::move(xs));
  }

  // Read products
  for (const auto& name : group_names(group)) {
    if (name.rfind("product_", 0) == 0) {
      hid_t pgroup = open_group(group, name.c_str());
      products_.emplace_back(pgroup);
      close_group(pgroup);
    }
  }

  // <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<
  // Before the secondary distribution refactor, when the angle/energy
  // distribution was uncorrelated, no angle was actually sampled. With
  // the refactor, an angle is always sampled for an uncorrelated
  // distribution even when no angle distribution exists in the ACE file
  // (isotropic is assumed). To preserve the RNG stream, we explicitly
  // mark fission reactions so that we avoid the angle sampling.
  if (is_fission(mt_)) {
    for (auto& p : products_) {
      if (p.particle_ == ParticleType::neutron) {
        for (auto& d : p.distribution_) {
          auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
          if (d_) d_->fission() = true;
        }
      }
    }
  }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

Reaction* reaction_from_hdf5(hid_t group, int* temperatures, int n)
{
  std::vector<int> temps {temperatures, temperatures + n};
  return new Reaction{group, temps};
}

void reaction_delete(Reaction* rx) { delete rx; }

int reaction_mt(Reaction* rx) { return rx->mt_; }

double reaction_q_value(Reaction* rx) { return rx->q_value_; }

bool reaction_scatter_in_cm(Reaction* rx) { return rx->scatter_in_cm_; }

double reaction_product_decay_rate(Reaction* rx, int product)
{
  return rx->products_[product - 1].decay_rate_;
}

int reaction_product_emission_mode(Reaction* rx, int product)
{
  switch (rx->products_[product - 1].emission_mode_) {
  case ReactionProduct::EmissionMode::prompt:
    return 1;
  case ReactionProduct::EmissionMode::delayed:
    return 2;
  case ReactionProduct::EmissionMode::total:
    return 3;
  }
}

int reaction_product_particle(Reaction* rx, int product)
{
  switch (rx->products_[product - 1].particle_) {
  case ParticleType::neutron:
    return 1;
  case ParticleType::photon:
    return 2;
  case ParticleType::electron:
    return 3;
  case ParticleType::positron:
    return 4;
  }
}

void reaction_product_sample(Reaction* rx, int product, double E_in, double* E_out, double* mu)
{
  rx->products_[product - 1].sample(E_in, *E_out, *mu);
}

double reaction_product_yield(Reaction* rx, int product, double E)
{
  return (*rx->products_[product - 1].yield_)(E);
}

int reaction_products_size(Reaction* rx) { return rx->products_.size(); }

double reaction_xs(Reaction* rx, int temperature, int energy)
{
  return rx->xs_[temperature - 1].value[energy - 1];
}

double reaction_sample_elastic_mu(Reaction* rx, double E)
{
  // Get elastic scattering distribution
  auto& d = rx->products_[0].distribution_[0];

  // Check if it is an uncorrelated angle-energy distribution
  auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
  if (d_) {
    return d_->angle().sample(E);
  } else {
    return 2.0*prn() - 1.0;
  }

}

int reaction_xs_size(Reaction* rx, int temperature)
{
  return rx->xs_[temperature - 1].value.size();
}

int reaction_xs_threshold(Reaction* rx, int temperature)
{
  return rx->xs_[temperature - 1].threshold;
}

}
