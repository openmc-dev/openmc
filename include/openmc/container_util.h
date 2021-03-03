#ifndef OPENMC_CONTAINER_UTIL_H
#define OPENMC_CONTAINER_UTIL_H

#include <algorithm> // for find
#include <iterator> // for begin, end

namespace openmc {

template<class C, class T>
inline bool contains(const C& v, const T& x)
{
  return std::end(v) != std::find(std::begin(v), std::end(v), x);
}

}

#endif // OPENMC_CONTAINER_UTIL_H
