#include "openmc/photon.h"

#include "openmc/hdf5_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {

std::vector<PhotonInteraction> elements;

} // namespace data

//==============================================================================
// PhotonInteraction implementation
//==============================================================================

PhotonInteraction::PhotonInteraction(hid_t group)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  read_attribute(group, "Z", Z_);
}

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" void photon_from_hdf5_c(hid_t group)
{
  data::elements.emplace_back(group);
}

} // namespace openmc
