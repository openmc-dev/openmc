//! \file r2s_source.cpp
//! \brief Decay photon source generation for R2S calculations

#include "openmc/r2s_source.h"

#include <algorithm> // for sort
#include <numeric>   // for accumulate
#include <string>

#include "openmc/capi.h"
#include "openmc/chain.h"
#include "openmc/distribution.h"
#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/error.h"
#include "openmc/memory.h"
#include "openmc/openmp_interface.h"
#include "openmc/particle_type.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// DecayPhotonMixture — non-owning mixture of decay photon distributions
//==============================================================================

//! Energy distribution formed by mixing multiple decay photon spectra.
//!
//! Unlike the general Mixture distribution, this class holds non-owning
//! pointers to the component distributions (which live in
//! data::chain_nuclides). Each component is weighted by the activity
//! (atoms * decay_constant) of the corresponding nuclide.

class DecayPhotonMixture : public Distribution {
public:
  //! Construct from non-owning distribution pointers and weights
  //!
  //! \param dists  Non-owning pointers to component distributions
  //! \param weights  Activity-based weights for each component. The Mixture
  //!   probability for component i is weights[i] * dists[i]->integral()
  //!   (the integral encodes photons-per-decay).
  DecayPhotonMixture(vector<const Distribution*> dists, vector<double> weights)
    : dists_(std::move(dists))
  {
    vector<double> probs;
    probs.reserve(dists_.size());
    for (size_t i = 0; i < dists_.size(); ++i) {
      probs.push_back(weights[i] * dists_[i]->integral());
    }
    integral_ = std::accumulate(probs.begin(), probs.end(), 0.0);
    di_.assign(probs);
  }

  std::pair<double, double> sample(uint64_t* seed) const override
  {
    size_t idx = di_.sample(seed);
    return dists_[idx]->sample(seed);
  }

  double integral() const override { return integral_; }

protected:
  double sample_unbiased(uint64_t* seed) const override
  {
    size_t idx = di_.sample(seed);
    return dists_[idx]->sample(seed).first;
  }

private:
  vector<const Distribution*> dists_; //!< Non-owning component distributions
  DiscreteIndex di_; //!< Discrete index for component selection
  double integral_;  //!< Total photon emission rate
};

//==============================================================================
// create_decay_photon_sources implementation
//==============================================================================

void create_decay_photon_sources(int n_regions, const int32_t* domain_ids,
  Source::DomainType domain_type, const double* lower_left,
  const double* upper_right, int n_nuclides, const char** nuclide_names,
  const double* atom_densities, const double* volumes)
{
  // Pre-resolve nuclide names to chain nuclide indices. A value of -1
  // indicates the nuclide is not in the chain (or has no decay photon data).
  vector<int> chain_indices(n_nuclides, -1);
  for (int j = 0; j < n_nuclides; ++j) {
    auto it = data::chain_nuclide_map.find(nuclide_names[j]);
    if (it != data::chain_nuclide_map.end()) {
      chain_indices[j] = it->second;
    }
  }

  // Clear existing sources
  model::external_sources.clear();

  // Thread-local storage for sources built in parallel
  vector<vector<unique_ptr<Source>>> thread_sources;

  int n_threads = num_threads();
  thread_sources.resize(n_threads);

#pragma omp parallel
  {
    int tid = thread_num();

#pragma omp for schedule(dynamic, 64)
    for (int i = 0; i < n_regions; ++i) {
      // Collect distributions and weights for nuclides present in this region
      vector<const Distribution*> dists;
      vector<double> weights;

      for (int j = 0; j < n_nuclides; ++j) {
        double atom_dens = atom_densities[i * n_nuclides + j];
        if (atom_dens <= 0.0)
          continue;

        int chain_idx = chain_indices[j];
        if (chain_idx < 0)
          continue;

        const auto& chain_nuc = data::chain_nuclides[chain_idx];
        const Distribution* photon_dist = chain_nuc->photon_energy();
        if (!photon_dist)
          continue;

        // Weight = number of atoms * decay_constant
        // atom_dens is in [atom/b-cm], volumes in [cm^3]
        // atoms = atom_dens * 1e24 * volume
        double atoms = atom_dens * 1.0e24 * volumes[i];
        double activity = atoms * chain_nuc->decay_constant();

        dists.push_back(photon_dist);
        weights.push_back(activity);
      }

      // Skip regions with no photon-emitting nuclides
      if (dists.empty())
        continue;

      // Build combined energy distribution
      auto energy =
        make_unique<DecayPhotonMixture>(std::move(dists), std::move(weights));
      double strength = energy->integral();
      if (strength <= 0.0)
        continue;

      // Build spatial distribution from bounding box
      Position ll {
        lower_left[i * 3], lower_left[i * 3 + 1], lower_left[i * 3 + 2]};
      Position ur {
        upper_right[i * 3], upper_right[i * 3 + 1], upper_right[i * 3 + 2]};

      // Build time distribution (delta at t=0)
      double t0[] {0.0};
      double tp[] {1.0};

      // Create IndependentSource
      auto source = make_unique<IndependentSource>(
        UPtrSpace {new SpatialBox(ll, ur)}, UPtrAngle {new Isotropic()},
        std::move(energy), UPtrDist {new Discrete(t0, tp, 1)});
      source->set_particle(ParticleType::photon());
      source->set_strength(strength);
      source->set_domain(domain_type, {domain_ids[i]});

      thread_sources[tid].push_back(std::move(source));
    }
  } // end parallel

  // Merge thread-local sources into model::external_sources
  for (auto& tvec : thread_sources) {
    for (auto& src : tvec) {
      model::external_sources.push_back(std::move(src));
    }
  }

  // Enable photon transport since we are creating photon sources
  settings::photon_transport = true;

  // Rebuild probability mass function for sampling
  vector<double> source_strengths;
  source_strengths.reserve(model::external_sources.size());
  for (auto& s : model::external_sources) {
    source_strengths.push_back(s->strength());
  }
  model::external_sources_probability.assign(source_strengths);
}

} // namespace openmc

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_create_decay_photon_sources(int n_regions,
  const int32_t* domain_ids, const char* domain_type, const double* lower_left,
  const double* upper_right, int n_nuclides, const char** nuclide_names,
  const double* atom_densities, const double* volumes)
{
  using DT = openmc::Source::DomainType;

  // Validate domain type
  std::string dt_str(domain_type);
  DT dt;
  if (dt_str == "material") {
    dt = DT::MATERIAL;
  } else if (dt_str == "cell") {
    dt = DT::CELL;
  } else {
    openmc::set_errmsg(
      "Invalid domain_type '" + dt_str + "'. Must be 'material' or 'cell'.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  try {
    openmc::create_decay_photon_sources(n_regions, domain_ids, dt, lower_left,
      upper_right, n_nuclides, nuclide_names, atom_densities, volumes);
  } catch (const std::exception& e) {
    openmc::set_errmsg(e.what());
    return OPENMC_E_DATA;
  }

  return 0;
}
