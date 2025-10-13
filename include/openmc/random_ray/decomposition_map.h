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
  void initialize();
  void generate_rank_centers(); //TODO: put in constructor
  // void calculate_grid_points(int grid_points_per_dimension);
  void calculate_grid_points(int grid_points_total);
  void initialize_points();
  void calculate_voronoi(vector<Position>& position_sum_per_rank, vector<int>& num_points_per_rank);
  Position calculate_centroids(const Position position_sum, const int num_points, int rank);

  // Methods to create and update subdomain list
  void update(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions);
  void exchange_sr_info(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions);
  bool any_discovered_source_regions(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions);
  // void transfer_sr_data(int sender, int receiver, SourceRegionKey sr_key, ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
  //   discovered_source_regions);
  void transfer_sr_data(int sender, int receiver, SourceRegion& sr_send, SourceRegion& sr_recv);
  void send_sr_data(int receiver, SourceRegion& sr_send);
  void receive_sr_data(int sender, SourceRegion& sr_recv);


  // Method to check if source region key is in subdomain of calling rank
  // bool is_SRK_in_domain(SourceRegionKey sr_key, Position r);
  // bool is_SRK_in_domain(SourceRegionKey sr_key, Position r,  
  //   ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
  //   discovered_source_regions);
  int find_owner(SourceRegionKey sr_key, Position r,  
    ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions, int64_t ID);
  // bool record_new_SRK(SourceRegionKey sr_key, Position r);
  // bool is_closest_rank(Position r);
  int find_closest_rank(Position r);
  // bool is_SRK_in_domain(SourceRegionKey sr_key);
  // int n_source_regions();

  // void calculate_rank_load(uint64_t total_geometric_intersections);

  // Calculates the load per rank based on the total number of hits all source regions
  // experience
  void calculate_rank_load(FlatSourceDomain* domain);
  // void calculate_rank_load(uint64_t total_geometric_intersections);

  void balance_load(FlatSourceDomain* domain);
  void update_load(FlatSourceDomain* domain);
  void redistribute_source_regions(FlatSourceDomain* domain);
  bool load_balanced() const { return max_imbalance_ < imbalance_tolerance_; }

  

  //----------------------------------------------------------------------------
  // Static data members


  //----------------------------------------------------------------------------
  // Public data members

  // Map that relates a SourceRegionKey to the index of the MPI rank that
  // contains that source region in its subdomain.
  std::unordered_map<SourceRegionKey, int, SourceRegionKey::HashFunctor>
    subdomain_map_;

  // Map that contains all newly discovered source region keys in a batch
  // that have not yet assigned a rank with certainty.
  // std::unordered_map<SourceRegionKey, int, SourceRegionKey::HashFunctor>
  //   discovered_regions_map_; //TODO: could exisiting discovered_source_regions container be re-used?

  // std::unordered_map<int, int> subdomain_map_;
  vector<Position> rank_centers_;

  vector<double> rank_load_;
  // std::unordered_set<SourceRegionKey, SourceRegionKey::HashFunctor> outside_source_regions_;
  // std::unordered_map<SourceRegionKey, int, SourceRegionKey::HashFunctor> outside_source_regions_;
  double target_load_;
  vector<double> rank_weights_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  // vector<int> my_subdomain_list_;
  SpatialBox* spatial_box_ = nullptr; // Add this member variable
  // std::unordered_map<int, vector<Position>> voronoi_map_;
  vector<Position> grid_points_;
  // int grid_points_per_rank_{2};
  int grid_points_per_rank_{100};
  // int grid_points_per_rank_{1000};
  int negroups_;
  uint64_t n_hits_sum_;

  double max_domain_length_;
  bool is_linear_;
  double max_imbalance_ = 0.0;
  double imbalance_tolerance_ = 0.01; // 1% imbalance tolerance
      // vector<Position> position_sum_per_rank(mpi::n_procs, Position(0,0,0));
      // vector<int> num_points_per_rank(mpi::n_procs, 0);

}; // class DecompositionMap

} // namespace openmc

#endif // OPENMC_DECOMPOSITION_MAP_H
