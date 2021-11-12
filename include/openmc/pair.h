#pragma once

namespace openmc {

#ifdef __CUDACC__
// This one will work on the GPU since we have the =default definition,
// which will produce both host and device code in CUDA.
template<typename First, typename Second>
struct pair {
  pair() = default;
  pair(pair const&) = default;
  pair(pair&&) = default;
  pair& operator=(pair const&) = default;
  pair& operator=(pair&&) = default;
  ~pair() = default;

  First first;
  Second second;
};
#else
#include <utility>
template<typename First, typename Second>
using pair = std::pair<First, Second>;
#endif

} // namespace openmc
