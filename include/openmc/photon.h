#ifndef OPENMC_PHOTON_H
#define OPENMC_PHOTON_H

#include <hdf5.h>

#include <string>
#include <vector>

namespace openmc {

//==============================================================================
//! Photon interaction data for a single element
//==============================================================================

class PhotonInteraction {
public:
  // Constructors
  PhotonInteraction(hid_t group);

  // Data members
  std::string name_; //! Name of element, e.g. "Zr"
  int Z_; //! Atomic number
};

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern std::vector<PhotonInteraction> elements;

} // namespace data

} // namespace openmc

#endif // OPENMC_PHOTON_H
