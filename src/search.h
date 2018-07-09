#ifndef OPENMC_SEARCH_H
#define OPENMC_SEARCH_H

#include <algorithm> // for lower_bound

namespace openmc {

template<class It, class T>
typename std::iterator_traits<It>::difference_type
lower_bound_index(It first, It last, const T& value)
{
  It index = std::lower_bound(first, last, value) - 1;
  return (index == last) ? -1 : index - first;
}

}

#endif // OPENMC_SEARCH_H