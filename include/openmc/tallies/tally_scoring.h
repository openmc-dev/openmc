#ifndef OPENMC_TALLIES_TALLY_SCORING_H
#define OPENMC_TALLIES_TALLY_SCORING_H

#include "openmc/particle.h"
#include "openmc/tallies/tally.h"

namespace openmc {

//==============================================================================
//! An iterator over all combinations of a tally's matching filter bins.
//
//! This iterator handles two distinct tasks.  First, it maps the N-dimensional
//! space created by the indices of N filters onto a 1D sequence.  In other
//! words, it provides a single number that uniquely identifies a combination of
//! bins for many filters.  Second, it handles the task of finding each valid
//! combination of filter bins given that each filter can have 1 or 2 or many
//! bins that are valid for the current tally event.
//==============================================================================

class FilterBinIter
{
public:

  //! Construct an iterator over bins that match a given particle's state.
  FilterBinIter(const Tally& tally, Particle* p);

  //! Construct an iterator over all filter bin combinations.
  //
  //! \param end if true, the returned iterator indicates the end of a loop.
  FilterBinIter(const Tally& tally, bool end);

  bool operator==(const FilterBinIter& other) const
  {return index_ == other.index_;}

  bool operator!=(const FilterBinIter& other) const
  {return !(*this == other);}

  FilterBinIter& operator++();

  int index_ {1};
  double weight_ {1.};

private:
  void compute_index_weight();

  const Tally& tally_;
};

//==============================================================================
// Non-member functions
//==============================================================================

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_SCORING_H
