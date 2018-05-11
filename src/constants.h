#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>


namespace openmc{

extern "C" double FP_COINCIDENT;
extern "C" double FP_PRECISION;
constexpr double INFTY {std::numeric_limits<double>::max()};
constexpr int C_NONE {-1};

} // namespace openmc

#endif // CONSTANTS_H
