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

const int line_in_header = 37;

namespace openmc {
namespace data {
//============================================================================
//! Stores all atomic masses (and mass excess and binding energy)
//============================================================================

/* Reads the Atomic mass data comes from the `Atomic Mass Evaluation 2020
    <https://doi.org/10.1088/1674-1137/abddaf>`_.
 */
class AtomicData {
public:
  //----------------------------------------------------------------------------
  // Constructors
  AtomicData(std::string data_file = "mass_1.mas20.txt");

  //==========================================================================
  // Methods and accessors
  double get_atomic_mass(std::string nuclide) const;
  double get_atomic_mass_excess(std::string nuclide) const;
  double get_atomic_binding_energy(std::string nuclide) const;
  atomic_data get_atomic_data(std::string nuclide) const;

private:
  //==========================================================================
  // Data members (accessor methods are above)
  std::map<std::string, atomic_data> atomic_mass_data;
};
} // namespace data
} // namespace openmc
#endif // OPENMC_ATOMIC_DATA_H
