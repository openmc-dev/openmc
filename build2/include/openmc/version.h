#ifndef OPENMC_VERSION_H
#define OPENMC_VERSION_H

#include "openmc/array.h"

namespace openmc {

// OpenMC major, minor, and release numbers
// clang-format off
constexpr int VERSION_MAJOR {0};
constexpr int VERSION_MINOR {15};
constexpr int VERSION_RELEASE {1};
constexpr bool VERSION_DEV {true};
constexpr const char* VERSION_COMMIT_COUNT = "277";
constexpr const char* VERSION_COMMIT_HASH = "b8d2cde33fbc53917f4a86596ae3ed6cec251e12";
constexpr std::array<int, 3> VERSION {VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE};
// clang-format on

} // namespace openmc

#endif // OPENMC_VERSION_H
