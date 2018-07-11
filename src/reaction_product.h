#ifndef OPENMC_REACTION_PRODUCT_H
#define OPENMC_REACTION_PRODUCT_H

#include <memory> // for unique_ptr
#include <vector> // for vector

#include "hdf5.h"
#include "angle_energy.h"
#include "endf.h"
#include "particle.h"

namespace openmc {

class ReactionProduct {
public:
  explicit ReactionProduct(hid_t group);
  void sample(double E_in, double& E_out, double& mu) const;

  enum class EmissionMode {
    prompt,  // Prompt emission of secondary particle
    total,   // Delayed emission of secondary particle
    delayed  // Yield represents total emission (prompt + delayed)
  };

  using Secondary = std::unique_ptr<AngleEnergy>;

  ParticleType particle_;
  EmissionMode emission_mode_;
  double decay_rate_;
  std::unique_ptr<Function1D> yield_;
  std::vector<Tabulated1D> applicability_;
  std::vector<Secondary> distribution_;
};

}

#endif // OPENMC_REACTION_PRODUCT_H
