#ifndef OPENMC_ATOMIC_DATA_H
#define OPENMC_ATOMIC_DATA_H

#include <map>
#include <string>

#include "openmc/constants.h"

struct ZAI_data {
  double mass_excess;
  double binding_energy;
  double mass;
};

const int line_in_header = 36;

namespace openmc {

//============================================================================
//! Defines
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
  std::map<std::string, ZAI_data> atomic_mass_data;

public:
  //==========================================================================
  // Methods and accessors
}
} // namespace openmc

#endif // OPENMC_ATOMIC_DATA_H
