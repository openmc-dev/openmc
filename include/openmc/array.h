#ifndef OPENMC_ARRAY_H
#define OPENMC_ARRAY_H

/*
 * See notes in include/openmc/vector.h
 *
 * In an implementation of OpenMC that uses an accelerator, we may remove the
 * use of array below and replace it with a custom
 * implementation behaving as expected on the device.
 */

#include <array>

namespace openmc {
using std::array;
} // namespace openmc

#endif // OPENMC_ARRAY_H
