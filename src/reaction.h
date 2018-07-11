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

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Reaction* reaction_from_hdf5(hid_t group, int* temperatures, int n);
  void reaction_delete(Reaction* rx);
  int reaction_mt(Reaction* rx);
  double reaction_q_value(Reaction* rx);
  bool reaction_scatter_in_cm(Reaction* rx);
  double reaction_product_decay_rate(Reaction* rx, int product);
  int reaction_product_emission_mode(Reaction* rx, int product);
  int reaction_product_particle(Reaction* rx, int product);
  void reaction_product_sample(Reaction* rx, int product, double E_in,
                               double* E_out, double* mu);
  int reaction_products_size(Reaction* rx);
  double reaction_product_yield(Reaction* rx, int product, double E);
  double reaction_sample_elastic_mu(Reaction* rx, double E);
  double reaction_xs(Reaction* xs, int temperature, int energy);
  int reaction_xs_size(Reaction* xs, int temperature);
  int reaction_xs_threshold(Reaction* xs, int temperature);
}

}

#endif // OPENMC_REACTION_H
