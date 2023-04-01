//! \file search.h
//! Search algorithms

#ifndef OPENMC_SEARCH_H
#define OPENMC_SEARCH_H

#include <algorithm> // for lower_bound, upper_bound

namespace openmc {

//! Perform binary search

template<class It, class T>
typename std::iterator_traits<It>::difference_type lower_bound_index(
  It first, It last, const T& value)
{
  if (*first == value)
    return 0;
  return std::lower_bound(first, last, value) - first - 1;
}

template<class It, class T>
typename std::iterator_traits<It>::difference_type upper_bound_index(
  It first, It last, const T& value)
{
  return std::upper_bound(first, last, value) - first - 1;
}

} // namespace openmc

#endif // OPENMC_SEARCH_H
