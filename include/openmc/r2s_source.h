//! \file r2s_source.h
//! \brief Decay photon source generation for R2S calculations

#ifndef OPENMC_R2S_SOURCE_H
#define OPENMC_R2S_SOURCE_H

#include <cstddef> // for size_t
#include <cstdint> // for int32_t

#include "openmc/distribution.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// DecayPhotonMixture — non-owning mixture of decay photon distributions
//==============================================================================

//! Energy distribution formed by mixing multiple decay photon spectra.
//!
//! Unlike the general Mixture distribution, this class holds non-owning
//! pointers to the component distributions (which live in
//! data::chain_nuclides). Each component is weighted by the activity
//! (atoms * decay_constant) of the corresponding nuclide.

class DecayPhotonMixture : public Distribution {
public:
  //! Construct from non-owning distribution pointers and weights
  //!
  //! \param dists  Non-owning pointers to component distributions
  //! \param weights  Activity-based weights for each component. The Mixture
  //!   probability for component i is weights[i] * dists[i]->integral()
  //!   (the integral encodes photons-per-decay).
  DecayPhotonMixture(vector<const Distribution*> dists, vector<double> weights);

  //! Construct from an XML node containing nuclide names and atom densities.
  //!
  //! Reads child ``<nuclide>`` elements with ``name`` and ``density``
  //! attributes, resolves them against the loaded depletion chain, and
  //! constructs the mixed distribution.
  explicit DecayPhotonMixture(pugi::xml_node node);

  std::pair<double, double> sample(uint64_t* seed) const override;
  double integral() const override;

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  vector<const Distribution*> dists_; //!< Non-owning component distributions
  DiscreteIndex di_; //!< Discrete index for component selection
  double integral_;  //!< Total photon emission rate
};

} // namespace openmc

#endif // OPENMC_R2S_SOURCE_H
