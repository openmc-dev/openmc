#include "openmc/static_map.h"

namespace openmc {

//==============================================================================
// Static map implementation
//==============================================================================

  /*
 `openmc::static_map<openmc::CellInstance, long, openmc::CellInstanceHash>::calc_num_buckets(unsigned long) const'
 `openmc::static_map<int, int, openmc::integer_hash<int> >::calc_num_buckets(unsigned long) const'
 */
//! Find the nearest prime above a value, used for hash table bucket sizing
/*
template<typename Key, typename T, typename Hash>
auto static_map<Key, T, Hash>::calc_num_buckets(size_type count) const
  -> size_type
{
  // Table of largest prime <= 2^n is pulled from https://oeis.org/A014234
  static const size_type primes[] = {3, 7, 13, 31, 61, 127, 251, 509, 1021,
    2039, 4093, 8191, 16381, 32749, 65521, 131071, 262139, 524287, 1048573,
    2097143, 4194301, 8388593, 16777213, 33554393, 67108859, 134217689,
    268435399, 536870909, 1073741789, 2147483647};

  // Find the first prime number larger than count in the table
  auto iter = std::lower_bound(std::begin(primes), std::end(primes), count);
  return *iter;
}
*/

} // namespace openmc
