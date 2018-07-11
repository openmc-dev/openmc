#ifndef OPENMC_REACTION_H
#define OPENMC_REACTION_H

#include <vector>

#include "hdf5.h"
#include "reaction_product.h"

namespace openmc {

class Reaction {
public:
  explicit Reaction(hid_t group, const std::vector<int>& temperatures);

  struct TemperatureXS {
    int threshold;
    std::vector<double> value;
  };

  int mt_;             //!< ENDF MT value
  double q_value_;     //!< Reaction Q value in [eV]
  bool scatter_in_cm_; //!< scattering system in center-of-mass?
  std::vector<TemperatureXS> xs_;
  std::vector<ReactionProduct> products_;
};

}

#endif // OPENMC_REACTION_H
