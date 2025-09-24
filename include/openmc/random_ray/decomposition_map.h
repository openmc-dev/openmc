#ifndef OPENMC_DECOMPOSITION_MAP_H
#define OPENMC_DECOMPOSITION_MAP_H

#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {

// // Container that contains the data to calculate the centroid 
// // of the Voronoi cell corresponding to each rank
// struct VoronoiCellData {
//   Position position;
//   int cell_count;
// };

class DecompositionMap;

namespace mpi {
    extern DecompositionMap decomp_map;
} // namespace mpi

class DecompositionMap {
public:
  //----------------------------------------------------------------------------
  // Constructors
  DecompositionMap();

  //----------------------------------------------------------------------------
  // Methods

  // Methods to find rank centres that divide spatial domain up into equal 
  // Voronoi volumes
  void generate_rank_centers(); //TODO: put in constructor
  void calculate_grid_points(int grid_points_per_dimension);
  void initialize_points();
  void calculate_voronoi(vector<Position>& position_sum_per_rank, vector<int>& num_points_per_rank);
  Position calculate_centroids(const Position position_sum, const int num_points, int rank);

  // Methods to create and update subdomain list
  void create(FlatSourceDomain* domain);
  void update();

  // Method to check if source region key is in subdomain of calling rank
  bool is_SRK_in_domain(int sr_key);
  // bool is_SRK_in_domain(SourceRegionKey sr_key);
  // int n_source_regions();

  //----------------------------------------------------------------------------
  // Static data members


  //----------------------------------------------------------------------------
  // Public data members

  // Map that relates a SourceRegionKey to the index of the MPI rank that
  // contains that source region in its subdomain.
  // std::unordered_map<SourceRegionKey, int, SourceRegionKey::HashFunctor>
  //   subdomain_map_;

  std::unordered_map<int, int> subdomain_map_;
  vector<Position> rank_centers_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  // vector<int> my_subdomain_list_;
  SpatialBox* spatial_box_ = nullptr; // Add this member variable
  // std::unordered_map<int, vector<Position>> voronoi_map_;
  vector<Position> grid_points_;
  int grid_points_per_rank_{10};
      // vector<Position> position_sum_per_rank(mpi::n_procs, Position(0,0,0));
      // vector<int> num_points_per_rank(mpi::n_procs, 0);

}; // class DecompositionMap

} // namespace openmc

#endif // OPENMC_DECOMPOSITION_MAP_H
