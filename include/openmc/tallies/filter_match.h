#ifndef OPENMC_TALLIES_FILTERMATCH_H
#define OPENMC_TALLIES_FILTERMATCH_H

#include <cassert>
#include <utility>

using std::pair;

namespace openmc {

/*
 * Defines the interface for storing instances
 * of matches against a Filter class's bins. For example,
 * if using a mesh filter and track length estimator,
 * each crossing of a mesh segment corresponds to a single
 * entry in the vector stored by FilterMatch. FilterMatch
 * is stored in a vector on each Particle class with indices
 * corresponding to each filter in the problem.
 *
 * For the vast majority of filters, only a single match
 * gets stored. Another example of a filter with multiple
 * possible matches are the functional expansions like
 * LegendreFilter.
 *
 * Because many, many types of filter only have a single
 * match, this class optimizes for that case. For example,
 * computationally expensive depletion calculations tend
 * to only employ filters which only match against one
 * bin. As such, a cache-friendliness optimization is that
 * we store the first bin/weight pair locally, then the
 * rest in a vector. Most of the time, the indirection
 * into the vector need not be used.
 *
 * In 2/2/2024, this sped up a depletion simulation
 * on a laptop by XX% and the improvement is similar
 * elsewhere.
 */
class FilterMatch {
public:
  // If using this interface, we assume this filter has
  // only single matches. Therefore we store the match
  // info locally in this object to avoid indirection.
  void set(int i)
  {
    w_ = 1.0;
    i_ = i;
  }
  void set(int i, double w)
  {
    w_ = w;
    i_ = i;
  }

  pair<int, double> get()
  {
    assert(!is_vector_);
    assert(i_ >= 0); // ensure it was set
    return {i_, w_};
  }

  // FilterBinIter uses this and increments it through
  // the tallying process.
  int& i_bin() { return i_bin_; }

  // Determine whether this match has been set already
  bool bins_present() { return is_vector_ || i_ >= 0; }

  /* When using this interface, it's assumed
   * that many possible matches are going to get
   * pushed back and the cache optimization technique
   * is abandoned. We therefore treat the starting
   * match as a zero weight entry, and then return
   * the vector storage for mesh to push back into.
   */
  vector<pair<int, double>>& vector_pairs()
  {
    is_vector_ = true;
    return bins_;
  }

  void reset()
  {
    if (is_vector_) {
      bins_.clear();
      is_vector_ = false;
    }

    w_ = 0.0;
    i_ = -1;
    i_bin_ = 0;
  }

  // These return the bin/weight at the current i_bin_ value
  int bin()
  {
    assert(!is_vector_ ? i_bin_ == 0 : true);
    if (is_vector_) {
      return bins_[i_bin_].first;
    } else {
      return i_;
    }
  }
  double weight()
  {
    assert(!is_vector_ ? i_bin_ == 0 : true);
    if (is_vector_) {
      return bins_[i_bin_].second;
    } else {
      return w_;
    }
  }

  int size()
  {
    if (is_vector_)
      return bins_.size();
    else
      return 1;
  }

private:
  // For scalar filter matches, these local variables
  // are used
  double w_ {0.0};
  int i_ {-1};

  int i_bin_ {0};

  // For vector filter matches, the vector
  // below gets used. It is best to avoid using
  // it in order to avoid on-the-fly allocations.
  vector<pair<int, double>> bins_;
  bool is_vector_ {false};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTERMATCH_H
