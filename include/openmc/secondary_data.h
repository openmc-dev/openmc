//! Native secondary-distribution resolution and immutable inspection snapshots.
#ifndef OPENMC_SECONDARY_DATA_H
#define OPENMC_SECONDARY_DATA_H

#include <cstdint>
#include <string>

#include "openmc/capi.h"
#include "openmc/vector.h"

namespace openmc {
class AngleEnergy;
class AngleDistribution;
class Nuclide;
class Reaction;
class ReactionProduct;

Nuclide& secondary_nuclide(int index);
Reaction& secondary_reaction(int index, int mt);
ReactionProduct& secondary_product(int index, int mt, int product);
AngleEnergy& secondary_distribution(int index, int mt, int product, int dist);

void clear_nuclear_data_backups();

//! Atomically replace PDFs on existing native angular grids/interpolation.
void replace_secondary_angle(AngleDistribution& angle,
  const vector<double>& incident, const vector<int>& offsets,
  const vector<int>& interpolation, const vector<double>& mu,
  const vector<double>& pdf);

//! A node owns its data; container children are indices in a snapshot.
struct SecondaryNode {
  std::string name;
  int type {OPENMC_SECONDARY_NULL};
  std::string text;
  vector<int64_t> integers;
  vector<double> reals;
  vector<size_t> children;
};

struct SecondarySnapshot {
  vector<SecondaryNode> nodes;
};

//! Inspect a native law, also used by loader-based unit tests.
SecondarySnapshot inspect_secondary(const AngleEnergy& distribution);
} // namespace openmc
#endif
