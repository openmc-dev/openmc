//! \file chain.h
//! \brief Depletion chain and associated information

#ifndef OPENMC_CHAIN_H
#define OPENMC_CHAIN_H

#include <cmath>
#include <string>
#include <unordered_map>

#include "pugixml.hpp"

#include "openmc/angle_energy.h" // for AngleEnergy
#include "openmc/distribution.h" // for UPtrDist
#include "openmc/endf.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Fission Yield Data
//==============================================================================

class FissionYields {
public:
  // Constructors, destructors
  FissionYields(pugi::xml_node node);

  // Data members
  std::unordered_map<std::string, Tabulated1D> yields_;
};

//==============================================================================
// Data for a nuclide in the depletion chain
//==============================================================================

class ChainNuclide {
public:
  // Types
  struct Product {
    std::string name;       //!< Reaction product name
    double branching_ratio; //!< Branching ratio
  };

  // Constructors, destructors
  ChainNuclide(pugi::xml_node node);
  ~ChainNuclide();

  //! Compute the decay constant for the nuclide
  //! \return Decay constant in [1/s]
  double decay_constant() const { return std::log(2.0) / half_life_; }
  FissionYields* fission_yields();

  const Distribution* photon_energy() const { return photon_energy_.get(); }
  const std::unordered_map<int, vector<Product>>& reaction_products() const
  {
    return reaction_products_;
  }

private:
  // Data members
  std::string name_;                        //!< Name of nuclide
  double half_life_ {0.0};                  //!< Half-life in [s]
  double decay_energy_ {0.0};               //!< Decay energy in [eV]
  FissionYields* fission_yields_ = nullptr; //!< Fission yields
  std::string fission_yields_parent_ {""};  //!< Fission yields parent
  std::unordered_map<int, vector<Product>>
    reaction_products_;    //!< Map of MT to reaction products
  UPtrDist photon_energy_; //!< Decay photon energy distribution
};

//==============================================================================
// Angle-energy distribution for decay photon
//==============================================================================

class DecayPhotonAngleEnergy : public AngleEnergy {
public:
  explicit DecayPhotonAngleEnergy(const Distribution* dist)
    : photon_energy_(dist)
  {}

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;

private:
  const Distribution* photon_energy_;
};

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern std::unordered_map<std::string, int> chain_nuclide_map;
extern vector<unique_ptr<ChainNuclide>> chain_nuclides;
extern vector<unique_ptr<FissionYields>> fission_yields;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

void read_chain_file_xml();

} // namespace openmc

#endif // OPENMC_CHAIN_H
