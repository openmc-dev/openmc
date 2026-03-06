//! \file r2s_source.h
//! \brief Decay photon source generation for R2S calculations

#ifndef OPENMC_R2S_SOURCE_H
#define OPENMC_R2S_SOURCE_H

#include <cstddef> // for size_t
#include <cstdint> // for int32_t

#include "openmc/source.h"

namespace openmc {

//! Create decay photon IndependentSource objects from activated materials
//!
//! For each spatial region, this function combines the decay photon energy
//! spectra from all nuclides present (weighted by atom count and decay
//! constant) and creates an IndependentSource with a SpatialBox spatial
//! distribution and isotropic angular distribution. The resulting sources
//! are stored in model::external_sources with the probability mass function
//! rebuilt accordingly.
//!
//! \param n_regions Number of spatial regions
//! \param domain_ids Material or cell IDs for source constraints [n_regions]
//! \param domain_type MATERIAL or CELL
//! \param lower_left Flattened lower-left bbox corners [n_regions * 3]
//! \param upper_right Flattened upper-right bbox corners [n_regions * 3]
//! \param n_nuclides Number of nuclides in the atom density matrix
//! \param nuclide_names Array of nuclide name strings [n_nuclides]
//! \param atom_densities Row-major atom densities in [atom/b-cm],
//!   shape [n_regions * n_nuclides]
//! \param volumes Volume of each region in [cm^3], shape [n_regions]
void create_decay_photon_sources(int n_regions, const int32_t* domain_ids,
  Source::DomainType domain_type, const double* lower_left,
  const double* upper_right, int n_nuclides, const char** nuclide_names,
  const double* atom_densities, const double* volumes);

} // namespace openmc

#endif // OPENMC_R2S_SOURCE_H
