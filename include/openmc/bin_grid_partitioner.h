#ifndef OPENMC_BIN_GRID_PARTITIONER_H
#define OPENMC_BIN_GRID_PARTITIONER_H

#include "hdf5_interface.h"
#include "partitioner_utils.h"
#include "universe.h"
#include "zplane_partitioner.h"
#include <vector>

namespace openmc {

// The bin grid partitioner divides the universe into a set of bins. It then
// searches for cells in those bins and stores the list of cells it finds for
// each bin.
class BinGridPartitioner : public UniversePartitioner {
public:
  explicit BinGridPartitioner(
    const Universe& univ, const AABB& bounds, int32_t grid_res = 32);
  explicit BinGridPartitioner(
    const Universe& univ, const AABB& bounds, hid_t file);
  virtual ~BinGridPartitioner() override;

  virtual void export_to_hdf5(const std::string& file_path) const override;

  //! Return the list of cells that could contain the given coordinates.
  virtual const std::vector<int32_t>& get_cells(
    Position r, Direction u) const override;
  virtual const std::vector<int32_t>& get_cells_fallback(
    Position r, Direction u) const override;

private:
  // bounds of the bin grid partitioner
  AABB bounds_;

  // Resolution of the bin grid on each axis
  int32_t grid_res_;

  // Dimensions of each bin
  Position bin_dim_;

  // The actual array of bins. Each bin is a vector of ints that denotes the
  // cells within that bin
  std::vector<std::vector<int>> bin_grid_;

  // Use a z-plane partitioner as a fallback
  ZPlanePartitioner fallback_;

  // This maps a position to a specific bin in bin_grid_
  int convert_to_index(const Position& r) const;
};

}; // namespace openmc

#endif