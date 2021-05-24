//! \file search.h
//! Search algorithms

#ifndef OPENMC_SEARCH_H
#define OPENMC_SEARCH_H

#include <algorithm> // for lower_bound, upper_bound

namespace openmc {

//! Perform binary search

template<class It, class T>
typename std::iterator_traits<It>::difference_type
lower_bound_index(It first, It last, const T& value)
{
  if (*first == value) return 0;
  It index = std::lower_bound(first, last, value) - 1;
  return (index == last) ? -1 : index - first;
}

template<class It, class T>
typename std::iterator_traits<It>::difference_type
upper_bound_index(It first, It last, const T& value)
{
  It index = std::upper_bound(first, last, value) - 1;
  return (index == last) ? -1 : index - first;
}

} // namespace openmc

#endif // OPENMC_SEARCH_H
