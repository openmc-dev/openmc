#ifndef OPENMC_GEOMETRY_H
#define OPENMC_GEOMETRY_H

#include <cstdint>
#include <vector>

namespace openmc {

extern "C" int openmc_root_universe;

//TODO: free this memory
extern std::vector<int64_t> overlap_check_count;

} // namespace openmc

#endif // OPENMC_GEOMETRY_H
