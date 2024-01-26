#ifndef OPENMC_MEMORY_H
#define OPENMC_MEMORY_H

/*
 * See notes in include/openmc/vector.h
 *
 * In an implementation of OpenMC that uses an accelerator, we may remove the
 * use of unique_ptr, etc. below and replace it with a custom
 * implementation behaving as expected on the device.
 */

#include <memory>

namespace openmc {
using std::make_unique;
using std::unique_ptr;
} // namespace openmc

#endif // OPENMC_MEMORY_H
