//! \file r2s_source.cpp
//! \brief Decay photon source generation for R2S calculations

#include "openmc/r2s_source.h"

#include <numeric> // for accumulate
#include <string>

#include "openmc/chain.h"
#include "openmc/error.h"
#include "openmc/memory.h"
#include "openmc/vector.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// DecayPhotonMixture implementation
//==============================================================================

DecayPhotonMixture::DecayPhotonMixture(
  vector<const Distribution*> dists, vector<double> weights)
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

DecayPhotonMixture::DecayPhotonMixture(pugi::xml_node node)
{
  // Read the region volume [cm^3] needed for absolute emission rate
  if (!check_for_node(node, "volume"))
    fatal_error("DecayPhotonMixture: 'volume' attribute is required.");
  double volume = std::stod(get_node_value(node, "volume"));

  // Read nuclide names and atom densities from XML
  vector<const Distribution*> dists;
  vector<double> weights;

  for (auto nuclide_node : node.children("nuclide")) {
    std::string name = get_node_value(nuclide_node, "name");
    double density = std::stod(get_node_value(nuclide_node, "density"));

    // Look up nuclide in the depletion chain
    auto it = data::chain_nuclide_map.find(name);
    if (it == data::chain_nuclide_map.end())
      continue;

    const auto& chain_nuc = data::chain_nuclides[it->second];
    const Distribution* photon_dist = chain_nuc->photon_energy();
    if (!photon_dist)
      continue;

    // Weight = atom_count * decay_constant
    // atom_count = atom_density [atom/b-cm] * 1e24 [b-cm/cm^3] * volume [cm^3]
    double atom_count = density * 1.0e24 * volume;
    double weight = atom_count * chain_nuc->decay_constant();
    if (weight <= 0.0)
      continue;

    dists.push_back(photon_dist);
    weights.push_back(weight);
  }

  dists_ = std::move(dists);

  // Build alias table from weighted probabilities
  vector<double> probs;
  probs.reserve(dists_.size());
  for (size_t i = 0; i < dists_.size(); ++i) {
    probs.push_back(weights[i] * dists_[i]->integral());
  }
  integral_ = std::accumulate(probs.begin(), probs.end(), 0.0);
  di_.assign(probs);
}

std::pair<double, double> DecayPhotonMixture::sample(uint64_t* seed) const
{
  size_t idx = di_.sample(seed);
  return dists_[idx]->sample(seed);
}

double DecayPhotonMixture::integral() const
{
  return integral_;
}

double DecayPhotonMixture::sample_unbiased(uint64_t* seed) const
{
  size_t idx = di_.sample(seed);
  return dists_[idx]->sample(seed).first;
}

} // namespace openmc
