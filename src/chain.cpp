//! \file chain.cpp
//! \brief Depletion chain and associated information

#include "openmc/chain.h"

#include <cstdlib> // for getenv
#include <memory>  // for make_unique
#include <string>  // for stod

#include <fmt/core.h>
#include <pugixml.hpp>

#include "openmc/distribution.h" // for distribution_from_xml
#include "openmc/error.h"
#include "openmc/reaction.h"
#include "openmc/xml_interface.h" // for get_node_value

namespace openmc {

//==============================================================================
// ChainNuclide implementation
//==============================================================================

ChainNuclide::ChainNuclide(pugi::xml_node node)
{
  name_ = get_node_value(node, "name");
  if (check_for_node(node, "half_life")) {
    half_life_ = std::stod(get_node_value(node, "half_life"));
  }
  if (check_for_node(node, "decay_energy")) {
    decay_energy_ = std::stod(get_node_value(node, "decay_energy"));
  }

  // Read reactions to store MT -> product map
  for (pugi::xml_node reaction_node : node.children("reaction")) {
    std::string rx_name = get_node_value(reaction_node, "type");
    if (!reaction_node.attribute("target"))
      continue;
    std::string rx_target = get_node_value(reaction_node, "target");
    double branching_ratio = 1.0;
    if (reaction_node.attribute("branching_ratio")) {
      branching_ratio =
        std::stod(get_node_value(reaction_node, "branching_ratio"));
    }
    int mt = reaction_type(rx_name);
    reaction_products_[mt].push_back({rx_target, branching_ratio});
  }

  for (pugi::xml_node source_node : node.children("source")) {
    auto particle = get_node_value(source_node, "particle");
    if (particle == "photon") {
      photon_energy_ = distribution_from_xml(source_node);
      break;
    }
  }

  // Set entry in mapping
  data::chain_nuclide_map[name_] = data::chain_nuclides.size();
}

ChainNuclide::~ChainNuclide()
{
  data::chain_nuclide_map.erase(name_);
}

//==============================================================================
// DecayPhotonAngleEnergy implementation
//==============================================================================

void DecayPhotonAngleEnergy::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  E_out = photon_energy_->sample(seed);
  mu = Uniform(-1., 1.).sample(seed);
}

double DecayPhotonAngleEnergy::sample_energy_and_pdf(
  double E_in, double mu, double& E_out, uint64_t* seed) const
{
  E_out = photon_energy_->sample(seed);
  return 0.5;
}

//==============================================================================
// Global variables
//==============================================================================

namespace data {

std::unordered_map<std::string, int> chain_nuclide_map;
vector<unique_ptr<ChainNuclide>> chain_nuclides;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

void read_chain_file_xml()
{
  char* chain_file_path = std::getenv("OPENMC_CHAIN_FILE");
  if (!chain_file_path) {
    return;
  }

  write_message(5, "Reading chain file: {}...", chain_file_path);

  pugi::xml_document doc;
  auto result = doc.load_file(chain_file_path);
  if (!result) {
    fatal_error(
      fmt::format("Error processing chain file: {}", chain_file_path));
  }

  // Get root element
  pugi::xml_node root = doc.document_element();

  for (auto node : root.children("nuclide")) {
    data::chain_nuclides.push_back(std::make_unique<ChainNuclide>(node));
  }
}

} // namespace openmc
