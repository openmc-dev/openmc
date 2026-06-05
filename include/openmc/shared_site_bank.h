#ifndef OPENMC_SHARED_SITE_BANK_H
#define OPENMC_SHARED_SITE_BANK_H

//! \file shared_site_bank.h
//! \brief Shared site bank data structure for storing SourceSite objects

#include "openmc/particle.h"
#include "openmc/shared_array.h"
#include "openmc/vector.h"

#include <cstdint>

namespace openmc {

//==============================================================================
//! SharedSiteBank - A thread-safe bank for storing SourceSite objects
//!
//! This class wraps a SharedArray<SourceSite> and provides additional
//! functionality for sorting sites deterministically and managing the bank
//! during Monte Carlo simulations. It is used for fission banks and secondary
//! particle banks where multiple threads need to append sites concurrently.
//==============================================================================

class SharedSiteBank {
public:
  //==========================================================================
  // Constructors

  //! Default constructor
  SharedSiteBank() = default;

  //==========================================================================
  // Methods

  //! Reserve space for a specified number of sites
  //! \param capacity Number of sites to allocate space for
  void reserve(int64_t capacity) { sites_.reserve(capacity); }

  //! Clear all sites from the bank
  void clear() { sites_.clear(); }

  //! Resize the bank to contain a specified number of sites
  //! \param size New size of the bank
  void resize(int64_t size) { sites_.resize(size); }

  //! Thread-safe append of a site to the bank
  //! \param site The SourceSite to append
  //! \return Index where site was stored, or -1 if bank is full
  int64_t thread_safe_append(const SourceSite& site)
  {
    return sites_.thread_safe_append(site);
  }

  //! Sort the bank deterministically using parent/progeny IDs
  //!
  //! Performs an O(n) sort on the bank by leveraging the parent_id and
  //! progeny_id fields of banked sites. Uses progeny_per_particle to
  //! determine starting indices for each parent.
  //!
  //! \param progeny_per_particle Vector tracking progeny count per parent
  void sort(vector<int64_t>& progeny_per_particle);

  //! Load balance sites across MPI ranks
  //!
  //! Distributes the sites currently in this bank across MPI ranks such that
  //! each rank ends up with approximately the same number of sites. This is
  //! used by both fission and secondary banks after particles have been
  //! accumulated locally on each rank.
  //!
  //! \param dest_bank Destination bank where load-balanced sites are written
  //! \param work_index Starting index for each rank's work assignment
  //! \param work_per_rank Number of particles per rank
  void balance_across_ranks(vector<SourceSite>& dest_bank,
    const vector<int64_t>& work_index, int64_t work_per_rank);

  //! Sample sites using uniform combing and write to destination
  //!
  //! Uses the uniform combing algorithm to sample exactly n_samples sites
  //! from the bank and writes them to dest_bank. This is used specifically
  //! for fission banks in eigenvalue calculations to select source sites
  //! for the next generation.
  //!
  //! \param dest_bank Destination bank where sampled sites are written
  //! \param n_samples Number of sites to sample
  void sample(vector<SourceSite>& dest_bank, int64_t n_samples);

  //! Sample sites from bank and load balance across MPI ranks
  //!
  //! Convenience method for fission banks that combines sampling via
  //! uniform combing with load balancing across ranks. This performs:
  //! 1. Uniform combing to sample exactly n_particles sites
  //! 2. Load balancing to distribute sites evenly across MPI ranks
  //!
  //! \param dest_bank Destination bank where sampled sites are written
  //! \param n_particles Number of particles to sample
  //! \param work_index Starting index for each rank's work assignment
  //! \param work_per_rank Number of particles per rank
  void sample_and_synchronize(vector<SourceSite>& dest_bank,
    int64_t n_particles, const vector<int64_t>& work_index,
    int64_t work_per_rank);

  //==========================================================================
  // Accessors

  //! Get number of sites currently in bank
  int64_t size() { return sites_.size(); }

  //! Get capacity of bank
  int64_t capacity() { return sites_.capacity(); }

  //! Get pointer to underlying data
  SourceSite* data() { return sites_.data(); }
  const SourceSite* data() const { return sites_.data(); }

  //! Check if bank is full
  bool full() const { return sites_.full(); }

  //! Array access operators
  SourceSite& operator[](int64_t i) { return sites_[i]; }
  const SourceSite& operator[](int64_t i) const { return sites_[i]; }

  //! Iterator support
  SourceSite* begin() { return sites_.begin(); }
  const SourceSite* cbegin() const { return sites_.cbegin(); }
  SourceSite* end() { return sites_.end(); }
  const SourceSite* cend() const { return sites_.cend(); }

private:
  //==========================================================================
  // Data members

  SharedArray<SourceSite> sites_; //!< Underlying storage for sites
};

} // namespace openmc

#endif // OPENMC_SHARED_SITE_BANK_H
