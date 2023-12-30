#ifndef OPENMC_ATOMIC_DATA_H
#define OPENMC_ATOMIC_DATA_H

#include <map>
#include <string>

#include "openmc/constants.h"

struct atomic_data {
  double mass_excess;
  double binding_energy;
  double mass;
};

const int line_in_header = 36;

namespace openmc {

//============================================================================
//! Stores all atomics masses (and mass excess and binding energy)
//============================================================================

/*
 *
 *
 *
 */
class AtomicData {
public:
  //----------------------------------------------------------------------------
  // Constructors
  AtomicData(std::string data_file);

private:
  //==========================================================================
  // Data members (accessor methods are below)
  std::map<std::string, atomic_data> atomic_mass_data;

public:
  //==========================================================================
  // Methods and accessors
  double get_atomic_mass(std::string nuclide) const;
  double get_atomic_mass_excess(std::string nuclide) const;
  double get_atomic_binding_energy(std::string nuclide) const;
  atomic_data get_atomic_data(std::string nuclide) const
}
} // namespace openmc

#endif // OPENMC_ATOMIC_DATA_H
