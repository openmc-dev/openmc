//! \file search.h
//! Search algorithms

#ifndef OPENMC_SEARCH_H
#define OPENMC_SEARCH_H

#include <algorithm> // for lower_bound, upper_bound

namespace openmc {

//! Perform binary search

template<class It, class T>
typename std::iterator_traits<It>::difference_type HD inline
  lower_bound_index(It first, It last, const T& value)
{
  if (value > *(last - 1))
    return last - 1 - first;
  It orig_first = first;
  while (last > first + 2) {
    // unsigned interpolation = static_cast<unsigned>(
    //   (value - *first) / (*(last - 1) - *first) * (last - first - 1));
    unsigned interpolation = static_cast<unsigned>(last - first - 1) / 2;
    if (!interpolation)
      interpolation++;
    else if (interpolation == last - first)
      interpolation--;
    It mid = first + interpolation;
    if (mid == first)
      mid++;
    else if (mid == last)
      mid--;
    if (*mid == value)
      return mid - orig_first - 1;
    else if (*mid < value)
      first = mid;
    else
      last = mid + 1;
  }
  return first - orig_first;
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
