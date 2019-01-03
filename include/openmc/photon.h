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
//! Cached microscopic photon cross sections for a particular element at the
//! current energy
//==============================================================================

struct ElementMicroXS {
  int index_grid; //!< index on element energy grid
  double last_E {0.0}; //!< last evaluated energy in [eV]
  double interp_factor; //!< interpolation factor on energy grid
  double total; //!< microscopic total photon xs
  double coherent; //!< microscopic coherent xs
  double incoherent; //!< microscopic incoherent xs
  double photoelectric; //!< microscopic photoelectric xs
  double pair_production; //!< microscopic pair production xs
};

//==============================================================================
// Global variables
//==============================================================================

namespace data {
extern std::vector<PhotonInteraction> elements;
} // namespace data

namespace simulation {
extern ElementMicroXS* micro_photon_xs;
#pragma omp threadprivate(micro_photon_xs)
} // namespace simulation

} // namespace openmc

#endif // OPENMC_PHOTON_H
