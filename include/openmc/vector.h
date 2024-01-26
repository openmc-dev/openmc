#ifndef OPENMC_VECTOR_H
#define OPENMC_VECTOR_H

/*
 * In an implementation of OpenMC that offloads computations to an accelerator,
 * we may need to provide replacements for standard library containers and
 * algorithms that have no native implementations on the device of interest.
 * Because some developers are currently in the process of creating such code,
 * introducing the below typedef lessens the amount of rebase conflicts that
 * happen as they rebase their code on OpenMC's develop branch.
 *
 * In an implementation of OpenMC that uses such an accelerator, we may remove
 * the use of vector below and replace it with a custom implementation
 * behaving as expected on the device.
 */

#include <vector>

namespace openmc {
using std::vector;
}

#endif // OPENMC_VECTOR_H
