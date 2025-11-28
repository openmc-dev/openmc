#ifndef OPENMC_DECOMPOSITION_MAP_H
#define OPENMC_DECOMPOSITION_MAP_H

#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {

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
  void generate_rank_centers(); //TODO: put in constructor?
  void calculate_grid_points(int grid_points_total);
  void initialize_voronoi_centers();
  void calculate_voronoi(vector<Position>& position_sum_per_rank, vector<int>& num_points_per_rank);
  Position calculate_centroids(const Position position_sum, const int num_points, int rank);

  // Methods to create and update subdomain list and exchange source region data
  void update(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions);
  void exchange_sr_info(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions);
  bool any_discovered_source_regions(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions);
  void send_sr_data(int receiver, SourceRegion& sr_send);
  void receive_sr_data(int sender, SourceRegion& sr_recv);

  // Methods for balancing the load between ranks
  void balance_load(FlatSourceDomain* domain);
  void update_load(FlatSourceDomain* domain, bool check_all_ranks, vector<double>& rank_load_combined, vector<double>& load_ratio);
  void redistribute_source_regions(FlatSourceDomain* domain);

  // Methods to find owner of source region
  int find_owner(SourceRegionKey sr_key, Position r,  
    ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions);
  int find_closest_rank(Position r, bool test_all_ranks);

  // Method to calculate the load per rank based on the total number of hits in all source regions
  // of a rank
  void calculate_rank_load(FlatSourceDomain* domain, double batch_transport_time);  
  double calculate_load_ratio(int rank);
  
  //----------------------------------------------------------------------------
  // Public data members

  // Map that relates a SourceRegionKey to the index of the MPI rank that
  // contains that source region in its subdomain.
  std::unordered_map<SourceRegionKey, int, SourceRegionKey::HashFunctor>
    subdomain_map_;

  // Neighbors of each rank's Voronoi cell
  std::unordered_set<int> my_neighbors;

  // Data to estimate rank loads
  vector<uint64_t> num_base_source_region_RT_tot_; // total number of base source region ray trace operations per base source region
  vector<uint64_t> num_mesh_bin_RT_tot_; // total number of mesh bin ray trace operations per base source region
  vector<uint64_t> num_base_source_region_RT_batch_; // number of base source region ray trace operations per base source region for current batch
  vector<uint64_t> num_mesh_bin_RT_batch_; // number of mesh bin ray trace operations per base source region for current batch
  vector<double> ray_tracing_cost_;
  vector<double> volume_base_sr_;
  vector<double> measured_rank_load_fractions_;

  // Load optimization
  vector<double> rank_weights_;
  double target_load_;
  double max_load_imbalance_measured_ = 0.0;
  int cnt_unconverged_optimizations_total_ = 0;
  int cnt_optimizations_total_ = 0;

private:
  //----------------------------------------------------------------------------
  // Private data members
  SpatialBox* spatial_box_ = nullptr; 

  // Voronoi cell calculation
  vector<Position> grid_points_; //TODO: Should this be only defined in generate_rank_centers?
  int grid_points_per_rank_{125}; // default 5x5x5 grid points per rank
  vector<Position> rank_centers_; // Centers of each rank's Voronoi cell

  // Load calculation
  vector<double> estimated_rank_load_fractions_;
  vector<double> estimated_rank_load_totals_;
  double estimated_load_sum_;
  vector<vector<double>> load_history_; // stores values of load per rank for last few batches
  int load_history_size_ = 10;
  int history_idx = 0;

  // coefficients for load calculation
  double C1_ = 1.0;
  double C2_ = 0.1; 
  double C3_ = 0.1; 

  // Load optimization
  double imbalance_tolerance_ = 0.01; // 1% imbalance tolerance
  double optimization_history_factor_ = 1.0; 
  int cnt_unconverged_optimizations_ = 0;

  // Miscellaneous
  uint64_t n_base_sr_ = 0;
  int negroups_;
  double max_domain_length_;
  bool is_linear_;

}; // class DecompositionMap

} // namespace openmc

#endif // OPENMC_DECOMPOSITION_MAP_H
