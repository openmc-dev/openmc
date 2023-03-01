//! \file chain.h
//! \brief Depletion chain and associated information

#ifndef OPENMC_CHAIN_H
#define OPENMC_CHAIN_H

#include <string>
#include <unordered_map>

#include "pugixml.hpp"

#include "openmc/distribution.h" // for UPtrDist
#include "openmc/memory.h"       // for unique_ptr
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Data for a nuclide in the depletion chain
//==============================================================================

class ChainNuclide {
public:
  // Constructors, destructors
  ChainNuclide(pugi::xml_node node);
  ~ChainNuclide();

  const Distribution* photon_energy() const { return photon_energy_.get(); }

private:
  // Data members
  std::string name_;
  double half_life_ {0.0};
  double decay_energy_ {0.0};
  // TODO: Decay modes
  // TODO: Reactions?
  UPtrDist photon_energy_;
};

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern std::unordered_map<std::string, int> chain_nuclide_map;
extern vector<unique_ptr<ChainNuclide>> chain_nuclides;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

void read_chain_file_xml();

} // namespace openmc

#endif // OPENMC_CHAIN_H
