#include "openmc/random_ray/decomposition_map.h"

#include "openmc/message_passing.h"
#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/random_lcg.h"
#include "openmc/mgxs_interface.h"
#include "openmc/timer.h"

#include "openmc/simulation.h"

namespace openmc {

namespace mpi {
    DecompositionMap decomp_map;
}

// Constructor
DecompositionMap::DecompositionMap() {}

void DecompositionMap::initialize(){
  negroups_ = data::mg.num_energy_groups_;
  rank_load_.resize(mpi::n_procs, 0.0);
  target_load_ = 1.0/mpi::n_procs;
  rank_weights_.resize(mpi::n_procs, 1.0);

  spatial_box_ = dynamic_cast<SpatialBox*>(
    dynamic_cast<IndependentSource*>(RandomRay::ray_source_.get())->space());

  double x_length = spatial_box_->upper_right().x - spatial_box_->lower_left().x;
  double y_length = spatial_box_->upper_right().y - spatial_box_->lower_left().y;
  double z_length = spatial_box_->upper_right().z - spatial_box_->lower_left().z;

  max_domain_length_ = sqrt(x_length*x_length + y_length*y_length + z_length*z_length);

  is_linear_ = RandomRay::source_shape_ != RandomRaySourceShape::FLAT; //TODO: SHould that be decided gloabbly or on sr by sr base?
}

void DecompositionMap::generate_rank_centers(){

  // Calculate grid points that are used for Voronoi cells
  int grid_points_total = grid_points_per_rank_ * mpi::n_procs;
  printf("Calculating %d grid points for Voronoi tessellation...\n", grid_points_total);
  calculate_grid_points(grid_points_total);

  // Initialize points with random positions
  initialize_points();

  double err = 1.0;
  double precision = 1e-3; // 0.001 cm
  int it = 0;  
  int max_iterations = 100;

  // Lloyd's algorithm
  // https://en.wikipedia.org/wiki/Lloyd%27s_algorithm
  while (err > precision && it < max_iterations)
  {
    // Reset error to determine maximum movement below
    err = 0.0;

    vector<Position> position_sum_per_rank(mpi::n_procs, Position(0,0,0));
    vector<int> num_points_per_rank(mpi::n_procs, 0);

    // Compute Voronoi cells by summing up all position values of mesh grid points 
    // that are closest to a Voronoi center
    calculate_voronoi(position_sum_per_rank, num_points_per_rank);

    for (int rank = 0; rank < mpi::n_procs; rank++) {

      if (mpi::master){
        //TODO: temporary
        printf("  Rank %d: Point (%f, %f, %f)\n", rank, rank_centers_[rank].x, rank_centers_[rank].y, rank_centers_[rank].z); 
      }

      // Calculate centroid of the cell
      Position centroid = calculate_centroids(position_sum_per_rank[rank], num_points_per_rank[rank], rank);

      // Calculate movement for convergence check
      double movement = (centroid - rank_centers_[rank]).norm();

      // Move rank center to centroid
      rank_centers_[rank] = centroid;

      // record maximum movement
      if (movement > err) {
        err = movement;
      }
    }

    it++;
  }

  if (mpi::master) {
    if (it == max_iterations) {
      warning("Lloyd's algorithm did not converge within the maximum number of iterations.");
    } else {
      printf("Lloyd's algorithm converged in %d iterations.\n", it);
    }
    printf("The following Voronoi centres are being used:\n");
    for (int rank = 0; rank < mpi::n_procs; rank++) {
      printf("  Rank %d: Point (%f, %f, %f)\n", rank, rank_centers_[rank].x, rank_centers_[rank].y, rank_centers_[rank].z);
      }
  }
  
}

void DecompositionMap::calculate_grid_points(int grid_points_total){

    // Calculate length along each dimension    
    vector <double> domain_length(3);
    domain_length[0] = spatial_box_->upper_right().x - spatial_box_->lower_left().x;
    domain_length[1] = spatial_box_->upper_right().y - spatial_box_->lower_left().y;
    domain_length[2] = spatial_box_->upper_right().z - spatial_box_->lower_left().z;

    double volume = domain_length[0] * domain_length[1] * domain_length[2];

    //TODO: temporary
    printf("Spatial box volume: %f \n", volume);

    // For each dimension, determine grid points along that direction based on aspect ratio: 
    // domain_length / volume^(1/3) = grid_points_dimension / grid_points_total^(1/3).
    // Check if any dimension is so distorted that it would only receive minimun of 1 grid 
    // point and flag that direction to correct total number of grid points.
    int excluded_dimension = -1;
    vector<int> grid_points_per_dimension(3);
    for (int i = 0; i < 3; i++){
      double grid_points_estimate = ((domain_length[i] / cbrt(volume)) * cbrt(grid_points_total));

      if (grid_points_estimate > 1){
        grid_points_per_dimension[i] = round(grid_points_estimate);
      } else {
        excluded_dimension = i;
        grid_points_per_dimension[i] = 1;
      }
    }

    // If problem is 2D, exclude z direction
    if (RandomRay::geom_dim_ == RandomRayGeomDim::TWO_DIM){
      excluded_dimension = 2;
      grid_points_per_dimension[2] = 1;
    }

    // If one dimension is excluded, recalculate grid points in other two dimensions based on area.
    if (excluded_dimension != -1){
      double area = 1.0;
      
      for (int i = 0; i < 3; i++){
        if (i == excluded_dimension) continue;
        area *= domain_length[i];
      }

      for (int i = 0; i < 3; i++){
        if (i == excluded_dimension) continue;
        grid_points_per_dimension[i] = round((domain_length[i] / sqrt(area)) * sqrt(grid_points_total));
      }
    } 

    //TODO: temporary
    printf("Grid points in x, y, z: %d, %d, %d \n", grid_points_per_dimension[0], grid_points_per_dimension[1], grid_points_per_dimension[2]);

    // Adjust grid points in each dimension to match actual total number of grid points to the number of grid points requested
    double new_total = grid_points_per_dimension[0] * grid_points_per_dimension[1] * grid_points_per_dimension[2];
    double adjustment;
    if (excluded_dimension != -1) {
      // When one dimension is fixed at 2, use square root for the other two dimensions
      adjustment = sqrt(grid_points_total / new_total);
    } else {
      // When all dimensions are used, use cubic root
      adjustment = cbrt(grid_points_total / new_total);
    }

    // Multiply adjustment factor
    for (int i = 0; i < 3; i++){
      if (i == excluded_dimension) continue;
      grid_points_per_dimension[i] = round(grid_points_per_dimension[i] * adjustment);
    }
  
    //TODO: temporary
    printf("Grid points in x, y, z: %d, %d, %d \n", grid_points_per_dimension[0], grid_points_per_dimension[1], grid_points_per_dimension[2]);

    // Calculate spacing between grid points in each dimension
    vector <double> delta_value(3, 0.0);

    for (int i = 0; i < 3; i++){
      if (grid_points_per_dimension[i] > 1){
        delta_value[i] = domain_length[i] / (grid_points_per_dimension[i] - 1);
      }
    }

    double x = spatial_box_->lower_left().x + domain_length[0] * 0.5;
    double y = spatial_box_->lower_left().y + domain_length[1] * 0.5;
    double z = spatial_box_->lower_left().z + domain_length[2] * 0.5;

    // Generate all grid points
    for (int i = 0; i < grid_points_per_dimension[0]; i++) {
        if (grid_points_per_dimension[0] > 1) {        
          x = spatial_box_->lower_left().x + i * delta_value[0];
        }
        for (int j = 0; j < grid_points_per_dimension[1]; j++) {
            if (grid_points_per_dimension[1] > 1) {        
              y = spatial_box_->lower_left().y + j * delta_value[1];
            }
            for (int k = 0; k < grid_points_per_dimension[2]; k++) {
                if (grid_points_per_dimension[2] > 1) {        
                  z = spatial_box_->lower_left().z + k * delta_value[2];
                }
                grid_points_.push_back({x, y, z});
            }
        }
    }
}

// Places random points in the spatial domain. 
// Each point corresponds to the initial center of a rank.
void DecompositionMap::initialize_points(){
  rank_centers_.resize(mpi::n_procs);

  uint64_t seed = openmc_get_seed();

  // Sample random positions to start with
  for (int rank = 0; rank < mpi::n_procs; rank++){

    double x = prn(&seed);
    double y = prn(&seed);
    double z = 0.0;
    if (RandomRay::geom_dim_ == RandomRayGeomDim::THREE_DIM){
      z = prn(&seed);
    } else{
      z = 0.5; // Mid-plane in z direction
    }

    Position xi {x, y, z};

    // make a small shift in position to avoid geometry floating point issues //TODO: necessary? Adopted from halton sampling
    Position shift {FP_COINCIDENT, FP_COINCIDENT, FP_COINCIDENT};
    rank_centers_[rank] = (spatial_box_->lower_left() + shift) +
            xi * ((spatial_box_->upper_right() - shift) - (spatial_box_->lower_left() + shift));
  }
}

// Determine the distance of each mesh grid point to all rank centers.
// Sum up the positions of all mesh grid points that are closest to a given rank center
// for computation of centroid later. Record number of grid points per rank center.
void DecompositionMap::calculate_voronoi(vector<Position>& position_sum_per_rank,
    vector<int>& num_points_per_rank){

      // Assign each point to the closest rank center
      #pragma omp parallel for schedule(static)         
        for (int p = 0; p < grid_points_.size(); p++) {
            Position point = grid_points_[p];
            int closest_rank = C_NONE;
            double min_distance = INFTY;
            
            // Find closest rank center
            for (int rank = 0; rank < mpi::n_procs; rank++) {
                double dist = (point - rank_centers_[rank]).norm();
                // Power Voronoi diagram uses sqaured distances
                dist *= dist;

                if (dist < min_distance) {
                    min_distance = dist;
                    closest_rank = rank;
                }
            }

            if (mpi::master && closest_rank == C_NONE) {
              fatal_error("Could not find closest rank for Voronoi cell point " + std::to_string(p) + ".");
            }
            
            // Accumulate point coordinates for the closest rank
            #pragma omp atomic
            position_sum_per_rank[closest_rank].x += point.x;
            #pragma omp atomic
            position_sum_per_rank[closest_rank].y += point.y;
            #pragma omp atomic
            position_sum_per_rank[closest_rank].z += point.z;

            // Record number of mesh grid points for closest rank
            #pragma omp atomic
            num_points_per_rank[closest_rank]++;
        }

  }

Position DecompositionMap::calculate_centroids(const Position position_sum, const int num_points, int rank){

  // check if any points have been recorded in rank //TODO: Message not precise enough, specify what is meant by "mesh".
  if (num_points == 0){
    fatal_error("Rank " + std::to_string(rank) + " has no Voronoi cell points. Mesh is too coarse.");
  }
  
  Position centroid = position_sum; 

  // Divide by number of points
  double n = static_cast<double>(num_points);
  centroid.x /= n;
  centroid.y /= n;
  centroid.z /= n;
  
  return centroid;
}

// Update subdomain list for each decomposition map. //TODO: This function seems obsolete and should be included into main code
void DecompositionMap::update(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

    // Check if any new regions discovered and if so, exchange discovered cell data between ranks
    if (any_discovered_source_regions(discovered_source_regions)){
      simulation::time_source_region_exchange.start();
      exchange_sr_info(discovered_source_regions);
      simulation::time_source_region_exchange.stop();
    }

}

bool DecompositionMap::any_discovered_source_regions(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

  int flag = 0;

  if(discovered_source_regions.begin() != discovered_source_regions.end()) {
    flag = 1;
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, mpi::intracomm);
  
  return flag > 0;
}

void DecompositionMap::exchange_sr_info(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

  // Communicate maps
  for (int rank = 0; rank < mpi::n_procs; rank++) {

    // Send size
    uint64_t bcast_size = 0;
    if (rank == mpi::rank) {
      for (const auto& pair : discovered_source_regions) { 
        // Only broadcast source regions that have non-zero volume, i.e. regions discovered during active phase of ray
        if (pair.second.scalars_.volume_ > 0.0) { //TODO: can this check be avoided?
          bcast_size++;
        }
      }
    }

    MPI_Bcast(&bcast_size, 1, MPI_UINT64_T, rank, mpi::intracomm); //TODO: MPI_UINT64_T moght not be available?

    if (bcast_size > 0) {

      vector<int64_t> local_base_ids(bcast_size);
      vector<int64_t> local_mesh_bins(bcast_size);
        
      if(rank == mpi::rank) {
        // fill in vectors to be sent
        int i = 0;
        for (const auto& pair : discovered_source_regions) {
          if (pair.second.scalars_.volume_ > 0.0) {
            SourceRegionKey sr_key = pair.first;
            local_base_ids[i] = sr_key.base_source_region_id;
            local_mesh_bins[i] = sr_key.mesh_bin;
            i++;
          }
        }
      }

      // Broadcast all data
      MPI_Bcast(local_base_ids.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);
      MPI_Bcast(local_mesh_bins.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);
      
      // Update subdomain map
      for (int j = 0; j < bcast_size; j++) {

        SourceRegionKey sr_key(local_base_ids[j], local_mesh_bins[j]);

        // check if already in map or not, i.e. has someone else discovered that region already?
        if (subdomain_map_.find(sr_key) == subdomain_map_.end()){
          subdomain_map_[sr_key] = rank;
        } else{
          int resident_rank = subdomain_map_[sr_key]; // current owner
          int challenger_rank = rank; // current broadcasting rank
          double resident_rank_load = rank_load_[resident_rank];
          double challenger_rank_load = rank_load_[challenger_rank];

          // If load of challenger rank is lower, assign source region to that rank,
          // otherwise resident keeps it
          int sender;
          int receiver;
          if (challenger_rank_load < resident_rank_load){
            subdomain_map_[sr_key] = challenger_rank;
            sender = resident_rank;
            receiver = challenger_rank;
          } else{
            // resident keeps it
            sender = challenger_rank;
            receiver = resident_rank;
          }

          // Broadcast n_hits such that each rank can update load balance
          int bcast_n_hits = 0;
          if (mpi::rank == sender) {
            SourceRegion& contested_sr = discovered_source_regions[sr_key];
            bcast_n_hits = contested_sr.scalars_.n_hits_;
          }
          MPI_Bcast(&bcast_n_hits, 1, MPI_INT, sender, mpi::intracomm);

          double load_change = bcast_n_hits / static_cast<double>(n_hits_sum_);
          rank_load_[sender] -= load_change;
          rank_load_[receiver] += load_change;

          // Communciate source region data and merge on receiver side
          if (mpi::rank == sender){

            SourceRegion& sr = discovered_source_regions[sr_key];
            send_sr_data(receiver, sr);

            // clear old source region data from discovered regions map
            discovered_source_regions.erase(sr_key);
          }
          if (mpi::rank == receiver){
            SourceRegion& sr = discovered_source_regions[sr_key];
            SourceRegion sr_recv(negroups_, is_linear_);
            receive_sr_data(sender, sr_recv);

            sr.merge(sr_recv, is_linear_);
          }
        }
      }
    }
  }

}

void DecompositionMap::send_sr_data(int receiver, SourceRegion& sr_send){

  int num_scalar_messages = 1;
  //! NOTE: update if new vector fields are added to SourceExchangeVectors struct
  int num_vector_messages = 4;
  if (is_linear_){  
    num_vector_messages += 4;
  }
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    num_vector_messages += 1;
  }
  int num_requests = num_scalar_messages + num_vector_messages;

  vector<MPI_Request> requests(num_requests);
  int req_idx = 0;

  // Send scalar data to receiver
  MPI_Isend(&sr_send.scalars_, sizeof(ScalarSourceRegionFields), MPI_BYTE, receiver, 1, mpi::intracomm, &requests[req_idx]);
  req_idx ++;

  // Send vector data to receiver
  // Tags hardcoded to avoid confusion if new fields are not added sequentially
  MPI_Isend(sr_send.scalar_flux_old_.data(), negroups_, MPI_DOUBLE, receiver, 2, mpi::intracomm, &requests[req_idx]); 
  req_idx ++;

  MPI_Isend(sr_send.scalar_flux_new_.data(), negroups_, MPI_DOUBLE, receiver, 3, mpi::intracomm, &requests[req_idx]); 
  req_idx ++;

  MPI_Isend(sr_send.source_.data(), negroups_, MPI_FLOAT, receiver, 4, mpi::intracomm, &requests[req_idx]); 
  req_idx ++;

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
      MPI_Isend(sr_send.external_source_.data(), negroups_, MPI_FLOAT, receiver, 5, mpi::intracomm, &requests[req_idx]); 
      req_idx ++;
  }

  MPI_Isend(sr_send.scalar_flux_final_.data(), negroups_, MPI_DOUBLE, receiver, 6, mpi::intracomm, &requests[req_idx]); 
  req_idx ++;

  if (is_linear_){
    MPI_Isend(sr_send.source_gradients_.data(), 3*negroups_, MPI_DOUBLE, receiver, 7, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;

    MPI_Isend(sr_send.flux_moments_old_.data(), 3*negroups_, MPI_DOUBLE, receiver, 8, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;

    MPI_Isend(sr_send.flux_moments_new_.data(), 3*negroups_, MPI_DOUBLE, receiver, 9, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;

    MPI_Isend(sr_send.flux_moments_t_.data(), 3*negroups_, MPI_DOUBLE, receiver, 10, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;
  }

  if (req_idx != num_requests){
    fatal_error(fmt::format("Number of MPI requests does not match number of messages sent."
      "Check if num_vector_messages is up to date in DecompositionMap::transfer_sr_data()."));
  }

  // Wait for all communication to complete
  MPI_Waitall(num_requests, requests.data(), MPI_STATUSES_IGNORE);

}

void DecompositionMap::receive_sr_data(int sender, SourceRegion& sr_recv){

  // Receive scalar data from sender
  MPI_Recv(&sr_recv.scalars_, sizeof(ScalarSourceRegionFields), MPI_BYTE, sender, 1, mpi::intracomm, MPI_STATUS_IGNORE);

  // Receive vector data from sender
  MPI_Recv(sr_recv.scalar_flux_old_.data(), negroups_, MPI_DOUBLE, sender, 2, mpi::intracomm, MPI_STATUS_IGNORE); 
  MPI_Recv(sr_recv.scalar_flux_new_.data(), negroups_, MPI_DOUBLE, sender, 3, mpi::intracomm, MPI_STATUS_IGNORE);     
  MPI_Recv(sr_recv.source_.data(), negroups_, MPI_FLOAT, sender, 4, mpi::intracomm, MPI_STATUS_IGNORE); 

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    MPI_Recv(sr_recv.external_source_.data(), negroups_, MPI_FLOAT, sender, 5, mpi::intracomm, MPI_STATUS_IGNORE); 
  }

  MPI_Recv(sr_recv.scalar_flux_final_.data(), negroups_, MPI_DOUBLE, sender, 6, mpi::intracomm, MPI_STATUS_IGNORE); 

  if (is_linear_){
    MPI_Recv(sr_recv.source_gradients_.data(), 3*negroups_, MPI_DOUBLE, sender, 7, mpi::intracomm, MPI_STATUS_IGNORE); 
    MPI_Recv(sr_recv.flux_moments_old_.data(), 3*negroups_, MPI_DOUBLE, sender, 8, mpi::intracomm, MPI_STATUS_IGNORE); 
    MPI_Recv(sr_recv.flux_moments_new_.data(), 3*negroups_, MPI_DOUBLE, sender, 9, mpi::intracomm, MPI_STATUS_IGNORE); 
    MPI_Recv(sr_recv.flux_moments_t_.data(), 3*negroups_, MPI_DOUBLE, sender, 10, mpi::intracomm, MPI_STATUS_IGNORE); 
  }

}

int DecompositionMap::find_owner(SourceRegionKey sr_key, Position r, 
    ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){
      
  // Check if SRK is in subdomain map
  auto it = subdomain_map_.find(sr_key);
  if (it != subdomain_map_.end()){
    return it->second;
  } 

  // Check if already recorded in newly discovered source regions
  discovered_source_regions.lock(sr_key);
  bool sr_key_discovered = discovered_source_regions.contains(sr_key);
  discovered_source_regions.unlock(sr_key);
  if (sr_key_discovered){
    return mpi::rank;
  }

  // If not found in either map, check which rank would own source 
  // region beased on location
  int closest_rank = find_closest_rank(r, true); //TODO: Maybe make condition if all ranks need to be check dependend on whether new neighbours were identified!
  return closest_rank;
}

int DecompositionMap::find_closest_rank(Position r, bool test_all_ranks) {
  // Determine which rank the position belongs to
  int closest_rank = C_NONE;
  double min_distance = INFTY;

  vector<int> test_ranks;
  // test_ranks = mpi::decomp_map.my_neighbors;


  if (test_all_ranks){
    test_ranks.resize(mpi::n_procs);
    std::iota(test_ranks.begin(), test_ranks.end(), 0); // fill with 0, 1, ..., n_procs-1
  } else {
    // convert unordered set to vector
    test_ranks=vector<int>(mpi::decomp_map.my_neighbors.begin(), mpi::decomp_map.my_neighbors.end());
    test_ranks.push_back(mpi::rank); // Always include own rank
  }
  
  // Find closest rank center
  // for (int rank = 0; rank < mpi::n_procs; rank++) {
  for (int rank : test_ranks) {
      double dist = (r - rank_centers_[rank]).norm();
      // Distance function corresponding to weighted power Voronoi diagram.
      dist = dist*dist - rank_weights_[rank];
      if (dist < min_distance) {
          min_distance = dist;
          closest_rank = rank;
      }
  }

  if (mpi::master && closest_rank == C_NONE) {
    fatal_error("Could not find closest rank for new source region at position (" 
      + std::to_string(r.x) + ", " + std::to_string(r.y) + ", " + std::to_string(r.z) + ").");
  }

  return closest_rank;
}

void DecompositionMap::calculate_rank_load(FlatSourceDomain* domain){ //TODO: This uses accumulated value, not batch-wise values!

  // Reset rank load
  std::fill(rank_load_.begin(), rank_load_.end(), 0.0);
  int64_t n_hits_rank = 0;
  
  // Add number of hits for each source region by going through all exisiting source regions.
  #pragma omp parallel for reduction(+ : n_hits_rank)
    for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
      n_hits_rank += domain->source_regions_.n_hits(sr);
    }

  // Add newly discovered source region hits.
  for (const auto & [sr_key, sr] : domain->discovered_source_regions_) {
    n_hits_rank += sr.scalars_.n_hits_;
  }

  //TODO: Temporary print-outs
  // if (mpi::master) {
  //   printf("Load distribution: ");
  // }

  // Calculate rank load based on number of source region hits
  uint64_t n_hits_total = 0;
  MPI_Allreduce(&n_hits_rank, &n_hits_total, 1, MPI_INT64_T, MPI_SUM, mpi::intracomm);
  double load = (n_hits_rank/static_cast<double>(n_hits_total));
  n_hits_sum_ = n_hits_total;
  
  //TODO: Temporary print-outs
  // if(mpi::master){
  //   printf("Total hits: %ld \n", n_hits_total);
  // }

  // Communicate load across ranks
  for (int rank = 0; rank < mpi::n_procs; rank++) {
    double bcast_load = 0.0;
    if (rank == mpi::rank) {
      bcast_load = load;
    }
    MPI_Bcast(&bcast_load, 1, MPI_DOUBLE, rank, mpi::intracomm);
    rank_load_[rank] = bcast_load;
    //TODO: Temporary print-outs
    // if (mpi::master) {
    //   printf("RANK %d: %.2f%% ", rank, rank_load_[rank]*100.0);
    // }
  }
  //TODO: Temporary print-outs
  // if (mpi::master) {
  //   printf("\n");
  // }

  // Calculate imbalance
  // vector<double> imbalance_per_rank(mpi::n_procs, 0.0);
  // for (int rank = 0; rank < mpi::n_procs; rank++) {
  //   imbalance_per_rank[rank] = std::abs(rank_load_[rank] - target_load_) / target_load_;
  // }
  // max_imbalance_ = *std::max_element(imbalance_per_rank.begin(), imbalance_per_rank.end());
  double max_load = *std::max_element(rank_load_.begin(), rank_load_.end());
  double avg_load = std::accumulate(rank_load_.begin(), rank_load_.end(), 0.0)/mpi::n_procs;
  max_imbalance_ = (max_load - avg_load) / avg_load;
}

void DecompositionMap::balance_load(FlatSourceDomain* domain){

  //TODO: The optimisation seems really messy
  cnt_optimizations_total_ ++;

  int max_iterations = 200;
  double max_imbalance = max_imbalance_;

  int it_outer = 0;
  double adaptation_factor = 1; //0.5 //0.1;
  double min_adaptation_factor = 0.01;
  double max_adaptation_factor = 2;
  double prev_imbalance = max_imbalance;
  double avg_rank_distance = max_domain_length_/cbrt(mpi::n_procs); // rough estimate of average distance between ranks //TODO: cubic root might not be adaquate for problems that are not box like
  double weight_scale = avg_rank_distance * avg_rank_distance * optimization_history_factor_; //TODO: maybe adjust dynamically dependent on wether previous batch has converged
  // double weight_scale = std::accumulate(rank_weights_.begin(), rank_weights_.end(), 0.0) / mpi::n_procs; // scale weights based on average load
  double beta = 0.6; // momentum damping
  bool check_all_ranks = true;

  // History tracking
  vector<double> imbalance_history;
  vector<vector<double>> weight_history;
  imbalance_history.push_back(max_imbalance);
  weight_history.push_back(rank_weights_);

  vector<double> weight_change(mpi::n_procs, 0.0);
  // vector<double> old_weights = rank_weights_;
  // double max_imbalance_old = max_imbalance_;


  // vector<double> imbalance_per_rank(mpi::n_procs, 0.0);

  std::unordered_map<int, vector<int>> sr_send; // contains regions to be send

  while (max_imbalance > imbalance_tolerance_ && it_outer < max_iterations){

    //TODO: temporary
    // if (mpi::master)  printf("MPI load balancing iteration %d, max. imbalance: %.2f%% \n", it_outer, max_imbalance*100.0);

    // Optimise rank load rank by rank
    // for (int rank = 0; rank < mpi::n_procs; rank++) {
    //   double correction = adaptation_factor * ((rank_load_[rank] - target_load_) / target_load_); // * weight_scale;
    //   rank_weights_[rank] -= correction;
    // }

    for (int rank = 0; rank < mpi::n_procs; rank++) {
      // weight_scale = std::accumulate(rank_weights_.begin(), rank_weights_.end(), 0.0) / mpi::n_procs;
      // imbalance_per_rank[rank] = (rank_load_[rank] - target_load_) / target_load_;
      // double corr = imbalance_per_rank[rank] * weight_scale;
      double corr = ((rank_load_[rank] - target_load_) / target_load_) * weight_scale;
      weight_change[rank] = beta * weight_change[rank] + (1.0 - beta) * corr; // keep some inertia from previous changes to prevent oscillations
      // double damping = std::min(1.0, max_imbalance / imbalance_tolerance_); // dampening factor based on how far we are from convergence
      double damping = std::clamp(max_imbalance / prev_imbalance, 0.1, 1.0); // dampening factor based on whether we are getting closer to convergence or not, prevents big jumps, if too big a change, more dampening is applied
      rank_weights_[rank] -= adaptation_factor * damping * weight_change[rank];
    }

    // if (it_outer > 0 && simulation::current_batch > 1){
    if (simulation::current_batch > 1){
      check_all_ranks = false;
    }

    update_load(domain, check_all_ranks);
    it_outer ++;

    double max_load = *std::max_element(rank_load_.begin(), rank_load_.end());
    double avg_load = std::accumulate(rank_load_.begin(), rank_load_.end(), 0.0) / mpi::n_procs;
    max_imbalance = (max_load - avg_load) / avg_load;

    // max_imbalance = 0.0;
    // for (int rank = 0; rank < mpi::n_procs; rank++) {
    //   double imbalance = std::abs(imbalance_per_rank[rank]);
    //   if (imbalance > max_imbalance) {
    //     max_imbalance = imbalance;
    //   }
    // }

    // Store imbalance history
    imbalance_history.push_back(max_imbalance);
    weight_history.push_back(rank_weights_);

    // Adaptive factor
    if (max_imbalance > prev_imbalance)
        adaptation_factor = std::max(adaptation_factor * 0.5, min_adaptation_factor); // don’t go below min
    else
        adaptation_factor = std::min(adaptation_factor * 1.05, max_adaptation_factor);     // stable - accelerate slightly

    prev_imbalance = max_imbalance;

  }

  // Check convergence and adjust history optimization factor depending on failure mode to enable better convergence in the following batch
  if (it_outer == max_iterations){
    if(mpi::master) {
      warning("MPI load balancing has not converged after "
              + std::to_string(max_iterations) + " iterations."); 
    }

    cnt_unconverged_optimizations_ ++;
    cnt_unconverged_optimizations_total_ ++;

    // Check if oscillating or simply slow convergence
    bool is_oscillating = false;
    int direction_changes = 0;
    
    for (int i = 1; i < imbalance_history.size(); i++) {
      // Calculate change
      double change = imbalance_history[i] - imbalance_history[i-1];
      
      // Count direction changes
      if (i > 1) {
        double prev_change = imbalance_history[i-1] - imbalance_history[i-2];
        if (change * prev_change < 0) {
          direction_changes++;
        }
      }
    }

    double oscillation_ratio = (double)direction_changes / (imbalance_history.size() - 2);
    
    if (oscillation_ratio > 0.4) {
      // decrease weight for next batch if oscillating, TODO: SHould this be safeguarded with max/ min value?
      optimization_history_factor_ = std::max(optimization_history_factor_ * 0.5, min_adaptation_factor);
    } else {
      // increase weight for faster convergence if too slow.
      optimization_history_factor_ = std::min(optimization_history_factor_ * 1.2, max_adaptation_factor);     // stable - accelerate slightly
    }

    // Revert to previous weights if load balancing did not improve and return
    auto min_it = std::min_element(imbalance_history.begin(), imbalance_history.end());
    int best_index = std::distance(imbalance_history.begin(), min_it);
    double best_imbalance = *min_it;
    rank_weights_ = weight_history[best_index];

    if (mpi::master){
      printf("Best imbalance during optimization was %.2f%% at iteration %d. New history factor: %.2f \n", 
        best_imbalance*100.0, best_index, optimization_history_factor_);
    }

    if (cnt_unconverged_optimizations_ == 5){
      cnt_unconverged_optimizations_ = 0;
      imbalance_tolerance_ = best_imbalance; // relax tolerance if not converging
      if (mpi::master){
        printf("Relaxing MPI load balancing tolerance to %.2f%% \n", imbalance_tolerance_*100.0);
      }
    }

    if (best_index == 0){
      return;
    }
  } else {
    cnt_unconverged_optimizations_ = 0;
    optimization_history_factor_ = 1.0; // reset history factor if converged
    if (mpi::master){
      printf("MPI load balancing converged after %d iterations. Max. imbalance: %.2f%% \n", it_outer, max_imbalance*100.0);
    }
  }

  simulation::time_load_balance_sr_transfer.start();
  redistribute_source_regions(domain);
  simulation::time_load_balance_sr_transfer.stop();

  //TODO: temporary
  // if (mpi::master){
  //   printf("Weights: ");
  //   for (int rank = 0; rank < mpi::n_procs; rank++) {
  //     printf("RANK %d: %.2f ", rank, rank_weights_[rank]);
  //   }
  //   printf("\n");
  // }
}

void DecompositionMap::update_load(FlatSourceDomain* domain, bool check_all_ranks){

  vector<uint64_t> n_hits(mpi::n_procs, 0);

  // Add up source region hits to respective rank depending of position of centroid of each source region
  #pragma omp parallel
  {
    // number of hits per thread
    vector<uint64_t> local_hits(mpi::n_procs, 0);

    #pragma omp for
      for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
        Position centroid = domain->source_regions_.centroid(sr);
        int owner = find_closest_rank(centroid, check_all_ranks);
        local_hits[owner] += domain->source_regions_.n_hits(sr);
      }

    // Combine results from different threads
    #pragma omp critical (combining_hits)
    {
      for (int i = 0; i < mpi::n_procs; i++) {
        n_hits[i] += local_hits[i];
      }
    }
  }

  // // Check if any new neighbors were discovered during this update
  // if (check_all_ranks){

  //   // Create a set from current neighbors for faster lookups
  //   std::unordered_set<int> existing_neighbors(
  //     mpi::decomp_map.my_neighbors.begin(), 
  //     mpi::decomp_map.my_neighbors.end()
  //   );

  //   for (int rank = 0; rank < mpi::n_procs; rank++) {
  //     // add to neighbor list if new neighbour found
  //     if (n_hits[rank] > 0 && rank != mpi::rank && existing_neighbors.count(rank) == 0){
  //       mpi::decomp_map.my_neighbors.push_back(rank);
  //     }
  //   }
  // }
  // Accumulate hits across all ranks
  MPI_Allreduce(MPI_IN_PLACE, n_hits.data(), mpi::n_procs, MPI_UINT64_T, MPI_SUM, mpi::intracomm);

  for (int rank = 0; rank < mpi::n_procs; rank++){
    rank_load_[rank] = (n_hits[rank]/static_cast<double>(n_hits_sum_));
  }
}

void DecompositionMap::redistribute_source_regions(FlatSourceDomain* domain) {

  std::unordered_map<int, vector<int>> sr_send_list;
  vector<int> num_sr_receiving(mpi::n_procs, 0);

  // local source region container that contains updated list
  SourceRegionContainer source_regions_new(negroups_, is_linear_);

  // Each rank identifies source regions that need to be transferred to new owner and updates subdomain map accordingly
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
    Position centroid = domain->source_regions_.centroid(sr);
    int owner = find_closest_rank(centroid, true);

    // If owner changed, write source region to list of outbound source regions
    if (owner != mpi::rank){
      sr_send_list[owner].push_back(sr);
    }
    // If owner did not change, add source region to new local container
    else {
      source_regions_new.push_back(domain->source_regions_.get_source_region_handle(sr));
    }
  }

  // Each rank communicates other ranks about ownership changes to update subdomain map
  for (int rank = 0; rank < mpi::n_procs; rank++) {

    // Send size
    int bcast_size = 0;
    if (rank == mpi::rank) {
      for (const auto& pair : sr_send_list) {
        bcast_size += pair.second.size();
      }
    }
    MPI_Bcast(&bcast_size, 1, MPI_INT, rank, mpi::intracomm);

    if (bcast_size > 0) {
      vector<int> rank_ids(bcast_size);
      vector<int64_t> base_ids(bcast_size);
      vector<int64_t> mesh_bins(bcast_size);
        
      // Current owner prepares communication and fills vectors to be sent
      if(rank == mpi::rank) {
        int i = 0;
        for (auto& pair : sr_send_list) {

          int receiver = pair.first;  // new owner
          vector<int>& sr_indices = pair.second;  // vector of source region indices
      
          // Iterate through all source regions for this key
          for (int sr_idx : sr_indices) {
            SourceRegionKey sr_key = domain->source_regions_.key(sr_idx);
            rank_ids[i] = receiver;
            base_ids[i] = sr_key.base_source_region_id;
            mesh_bins[i] = sr_key.mesh_bin;
            i++;
          }
        }
      }

      // Broadcast source region data
      MPI_Bcast(rank_ids.data(), bcast_size, MPI_INT, rank, mpi::intracomm);
      MPI_Bcast(base_ids.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);
      MPI_Bcast(mesh_bins.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);

      // Every rank updates subdomain map
      for (int j = 0; j < bcast_size; j++) {

        // Re-assemble source region key and extract new owner
        SourceRegionKey sr_key(base_ids[j], mesh_bins[j]);
        int rank_new = rank_ids[j];

        // Update subdomain_map with new owner
        subdomain_map_[sr_key] = rank_new;

        // If calling rank is new owner, increase count of source regions coming from sending rank (current broadcaster)
        if (mpi::rank == rank_new){
          num_sr_receiving[rank] += 1;
        }
      }
    }
  }

  // clear source_region_map_
  domain->source_region_map_.clear();

  // Send source region data to new owner
  for (auto& pair : sr_send_list) {
    int receiver = pair.first;              // destination rank
    vector<int>& sr_indices = pair.second;  // vector of source region indices

    // Iterate through all source regions for this key
    for (int sr_idx : sr_indices) {
      SourceRegionKey sr_key = domain->source_regions_.key(sr_idx);
      SourceRegionHandle srh = domain->source_regions_.get_source_region_handle(sr_idx);
      SourceRegion sr(srh);
      send_sr_data(receiver, sr);
    }
  }

  // // Update source regions in domain
  // domain->source_regions_ = source_regions_new;
  
  // int64_t start_sr_id = domain->source_regions_.n_source_regions();
  int64_t start_sr_id = source_regions_new.n_source_regions();

  // Receive source regions
  for (int sender = 0; sender < mpi::n_procs; sender++) {

    if (sender == mpi::rank) {
      if (num_sr_receiving[sender] > 0){
        fatal_error("Rank sends source regions to itself. Rank should not receive source regions from itself.");
      }
      continue; // skip self
    }

    int num_sr = num_sr_receiving[sender];
    for (int i = 0; i < num_sr; ++i) {
      SourceRegion sr_recv(negroups_, is_linear_);
      receive_sr_data(sender, sr_recv); //TODO: Use <MPI tags?
      source_regions_new.push_back(sr_recv);
      // domain->source_regions_.push_back(sr_recv);
    }
  }

  // // Reinitialise tallies
  // domain->convert_source_regions_to_tallies(start_sr_id);

  // Update source regions in domain
  domain->source_regions_ = source_regions_new;

  // Update source region map
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
    SourceRegionKey key = domain->source_regions_.key(sr);
    domain->source_region_map_[key] = sr;
  }

  // Reinitialise tallies
  domain->convert_source_regions_to_tallies(start_sr_id);
}

}// namespace openmc