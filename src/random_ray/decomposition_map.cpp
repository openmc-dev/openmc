#include "openmc/random_ray/decomposition_map.h"

#include "openmc/message_passing.h"
#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/random_lcg.h"
#include "openmc/mgxs_interface.h"
#include "openmc/timer.h"
#include "openmc/cell.h"


#include "openmc/simulation.h"

namespace openmc {

namespace mpi {
    DecompositionMap decomp_map;
}

// Constructor
DecompositionMap::DecompositionMap() {}

void DecompositionMap::initialize(){
  negroups_ = data::mg.num_energy_groups_;
  estimated_rank_load_fractions_.resize(mpi::n_procs, 0.0);
  measured_rank_load_fractions_.resize(mpi::n_procs, 0.0);
  estimated_rank_load_totals_.resize(mpi::n_procs, 0.0);
  target_load_ = 1.0/mpi::n_procs;
  rank_weights_.resize(mpi::n_procs, 1.0);

  spatial_box_ = dynamic_cast<SpatialBox*>(
    dynamic_cast<IndependentSource*>(RandomRay::ray_source_.get())->space());

  double x_length = spatial_box_->upper_right().x - spatial_box_->lower_left().x;
  double y_length = spatial_box_->upper_right().y - spatial_box_->lower_left().y;
  double z_length = spatial_box_->upper_right().z - spatial_box_->lower_left().z;

  max_domain_length_ = sqrt(x_length*x_length + y_length*y_length + z_length*z_length);

  is_linear_ = RandomRay::source_shape_ != RandomRaySourceShape::FLAT; 

  // Count the number of source regions, compute the cell offset
  // indices, and store the material type The reason for the offsets is that
  // some cell types may not have material fills, and therefore do not
  // produce FSRs. Thus, we cannot index into the global arrays directly
  // int base_source_regions = 0;
  for (const auto& c : model::cells) {
    if (c->type_ == Fill::MATERIAL) {
      n_base_sr_ += c->n_instances();
    }
  }

  num_base_source_region_RT_tot_.resize(n_base_sr_, 0);
  num_mesh_bin_RT_tot_.resize(n_base_sr_, 0);
  num_base_source_region_RT_batch_.resize(n_base_sr_, 0);
  num_mesh_bin_RT_batch_.resize(n_base_sr_, 0);
  ray_tracing_cost_.resize(n_base_sr_);
  volume_base_sr_.resize(n_base_sr_, 0.0);

  load_history_.resize(mpi::n_procs);
  for (auto& row : load_history_) {
      row.resize(load_history_size_, 0.0);
  }
}

void DecompositionMap::generate_rank_centers(){

  // Calculate grid points that are used for Voronoi cells
  int grid_points_total = grid_points_per_rank_ * mpi::n_procs;
  if (mpi::master)
    printf("Calculating %d grid points for Voronoi tessellation...\n", grid_points_total);
  calculate_grid_points(grid_points_total);

  // Initialize points with random positions
  initialize_voronoi_centers();

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
    // printf("The following Voronoi centres are being used:\n");
    // for (int rank = 0; rank < mpi::n_procs; rank++) {
    //   printf("  Rank %d: Point (%f, %f, %f)\n", rank, rank_centers_[rank].x, rank_centers_[rank].y, rank_centers_[rank].z);
    //   }
  }
  
}

void DecompositionMap::calculate_grid_points(int grid_points_total){

    // Calculate length along each dimension    
    vector <double> domain_length(3);
    domain_length[0] = spatial_box_->upper_right().x - spatial_box_->lower_left().x;
    domain_length[1] = spatial_box_->upper_right().y - spatial_box_->lower_left().y;
    domain_length[2] = spatial_box_->upper_right().z - spatial_box_->lower_left().z;

    double volume = domain_length[0] * domain_length[1] * domain_length[2];

    // For each dimension, determine grid points along that direction based on aspect ratio: 
    // domain_length / volume^(1/3) = grid_points_dimension / grid_points_total^(1/3).
    // Check if any dimension is so distorted that it would only receive minimum of 1 grid 
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

    // Adjust grid points in each dimension to match actual total number of grid points to the number of grid points requested
    double new_total = grid_points_per_dimension[0] * grid_points_per_dimension[1] * grid_points_per_dimension[2];
    double adjustment;
    if (excluded_dimension != -1) {
      // When one dimension is excluded, use square root for the other two dimensions
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
  
    // Calculate spacing between grid points in each dimension
    vector <double> delta_value(3, 0.0);

    for (int i = 0; i < 3; i++){
      if (grid_points_per_dimension[i] > 1){
        delta_value[i] = (domain_length[i] - 2*TINY_BIT) / (grid_points_per_dimension[i] - 1);
      }
    }

    // Initialize point at center of domain
    double x = spatial_box_->lower_left().x + domain_length[0] * 0.5;
    double y = spatial_box_->lower_left().y + domain_length[1] * 0.5;
    double z = spatial_box_->lower_left().z + domain_length[2] * 0.5;

    // Generate all grid points
    for (int i = 0; i < grid_points_per_dimension[0]; i++) {
        if (grid_points_per_dimension[0] > 1) {        
          x = spatial_box_->lower_left().x + TINY_BIT + i * delta_value[0];
        }
        for (int j = 0; j < grid_points_per_dimension[1]; j++) {
            if (grid_points_per_dimension[1] > 1) {        
              y = spatial_box_->lower_left().y + TINY_BIT + j * delta_value[1];
            }
            for (int k = 0; k < grid_points_per_dimension[2]; k++) {
                if (grid_points_per_dimension[2] > 1) {        
                  z = spatial_box_->lower_left().z + TINY_BIT + k * delta_value[2];
                }
                // Add grid point
                grid_points_.push_back({x, y, z});
            }
        }
    }

    // Check if mesh grid points are inside spatial domain, if not erase them
    //TODO: Maybe all grid points should just be sampled randomly inside the domain to avoid erasing points and decreasing total number of points?
    for (int i = grid_points_.size() - 1; i >= 0; i--){
      Position xi = grid_points_[i];

      bool is_inside_domain = RandomRay::ray_source_->satisfies_spatial_constraints(xi);

      if (!is_inside_domain){
        grid_points_.erase(grid_points_.begin() + i);
      }
    }

    if (mpi::master && grid_points_.size() < grid_points_total) {
      warning(fmt::format(
        "Spatial constraints reduced grid points for Voronoi tesselation from {} to {}.",
        grid_points_total, grid_points_.size()));
    }
}

// Places random points in the spatial domain. 
// Each point corresponds to the initial center of a rank.
void DecompositionMap::initialize_voronoi_centers(){
  rank_centers_.resize(mpi::n_procs);

  uint64_t seed = openmc_get_seed();
  int rank_cnt = 0;

  // Sample random positions to start with
  while(rank_cnt < mpi::n_procs){

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
    xi = (spatial_box_->lower_left() + shift) +
            xi * ((spatial_box_->upper_right() - shift) - (spatial_box_->lower_left() + shift));

    bool is_inside_domain = RandomRay::ray_source_->satisfies_spatial_constraints(xi) ;

    if (is_inside_domain){
      rank_centers_[rank_cnt] = xi;
      rank_cnt++;
    }

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
                // Power Voronoi diagram uses squared distances
                dist = dist*dist - rank_weights_[rank];

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

  // check if any points have been recorded in rank
  if (num_points == 0){
    fatal_error("Rank " + std::to_string(rank) + " has no Voronoi cell points. This indicates that the number of grid points for the Voronoi tesselation is too coarse. Requires source code modificaiton to fix.");
  }
  
  Position centroid = position_sum; 

  // Divide by number of points
  double n = static_cast<double>(num_points);
  centroid.x /= n;
  centroid.y /= n;
  centroid.z /= n;
  
  return centroid;
}

bool DecompositionMap::any_discovered_source_regions(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

  simulation::time_decomposition_handling.start();

  int flag = 0;
  if(discovered_source_regions.begin() != discovered_source_regions.end()) {
    flag = 1;
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, mpi::intracomm);

  return flag > 0;
  
  simulation::time_decomposition_handling.start();
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
        if (pair.second.scalars_.volume_ > 0.0) {
          bcast_size++;
        }
      }
    }

    MPI_Bcast(&bcast_size, 1, MPI_UINT64_T, rank, mpi::intracomm); 

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

          double ratio_resident = calculate_load_ratio(resident_rank);
          double ratio_challenger = calculate_load_ratio(challenger_rank);
          double resident_rank_load = ratio_resident * estimated_rank_load_totals_[resident_rank];
          double challenger_rank_load = ratio_challenger * estimated_rank_load_totals_[challenger_rank];

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

          // Broadcast load of exchanged source region such that each rank can update load balance
          double bcast_load = 0;
          if (mpi::rank == sender) {
            SourceRegion& contested_sr = discovered_source_regions[sr_key];
            double volume_sr = contested_sr.scalars_.volume_;
            bcast_load = C1_ * (contested_sr.scalars_.n_hits_/simulation::current_batch) * negroups_ + volume_sr * ray_tracing_cost_[sr_key.base_source_region_id];
          }
          MPI_Bcast(&bcast_load, 1, MPI_DOUBLE, sender, mpi::intracomm);

          double load_change_fraction = bcast_load / estimated_load_sum_; // load fraction update
          estimated_rank_load_fractions_[sender] -= load_change_fraction;
          estimated_rank_load_fractions_[receiver] += load_change_fraction;
          double load_change_total = bcast_load; // load total update
          estimated_rank_load_totals_[sender] -= load_change_total;
          estimated_rank_load_totals_[receiver] += load_change_total;

          // Communicate source region data and merge on receiver side
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
      "Check if num_vector_messages corresponds to number of transferred vectors."));
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

  // If not found in either map, check which rank owns source 
  // region beased on location
  int closest_rank = find_closest_rank(r, true);
  return closest_rank;
}

int DecompositionMap::find_closest_rank(Position r, bool test_all_ranks) {

  int closest_rank = C_NONE;
  double min_distance = INFTY;
  vector<int> test_ranks;

  if (test_all_ranks){
    test_ranks.resize(mpi::n_procs);
    std::iota(test_ranks.begin(), test_ranks.end(), 0); // fill with 0, 1, ..., n_procs-1
  } else {
    // convert unordered set to vector
    test_ranks=vector<int>(mpi::decomp_map.my_neighbors.begin(), mpi::decomp_map.my_neighbors.end());
    test_ranks.push_back(mpi::rank); // Always include own rank
  }
  
  // Find closest rank center
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

void DecompositionMap::calculate_rank_load(FlatSourceDomain* domain, double batch_transport_time){

  // Reset local volumes of base source regions, which might change when source regions change rank ownership
  std::fill(volume_base_sr_.begin(), volume_base_sr_.end(), 0.0);

  num_base_source_region_RT_tot_ = num_base_source_region_RT_batch_;
  num_mesh_bin_RT_tot_ = num_mesh_bin_RT_batch_;

  // Add volumes of newly discovered source regions
  vector<uint64_t> mesh_bins_per_base_sr_local(n_base_sr_, 0);
  for (const auto & [sr_key, sr] : domain->discovered_source_regions_) {
    volume_base_sr_[sr_key.base_source_region_id] += sr.scalars_.volume_;
  }

  // Add volumes of known source regions
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
    SourceRegionKey sr_key = domain->source_regions_.key(sr);
    uint64_t base_sr = sr_key.base_source_region_id;
    volume_base_sr_[base_sr] += domain->source_regions_.volume_t(sr);
    volume_base_sr_[base_sr] += domain->source_regions_.volume(sr);
  }

  // Calculate ray tracing cost per base source region
  #pragma omp parallel for
    for (uint64_t bsr = 0; bsr < n_base_sr_; bsr++) {
      if (volume_base_sr_[bsr] > 0.0) {
        ray_tracing_cost_[bsr] = (C2_ * static_cast<double>(num_base_source_region_RT_tot_[bsr]) +  C3_ * static_cast<double>(num_mesh_bin_RT_tot_[bsr]))/volume_base_sr_[bsr];
      } else {
        ray_tracing_cost_[bsr] = 0.0;
      }
    }
  
  // Accumulate load in known source regions
  double local_estimated_load = 0.0;
  #pragma omp parallel for reduction(+: local_estimated_load)
    for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
      SourceRegionKey sr_key = domain->source_regions_.key(sr);
      uint64_t base_sr = sr_key.base_source_region_id;
      double volume_sr = domain->source_regions_.volume_t(sr) + domain->source_regions_.volume(sr);
      double load_sr = C1_ * (domain->source_regions_.n_hits(sr)/simulation::current_batch) * negroups_ + volume_sr * ray_tracing_cost_[base_sr];
      local_estimated_load += load_sr;
    }

  // Add load of newly discovered source regions
  double load_sr;
  for (const auto & [sr_key, sr] : domain->discovered_source_regions_) {
    uint64_t base_sr = sr_key.base_source_region_id;
    double volume_sr = sr.scalars_.volume_;
    load_sr = C1_ * (sr.scalars_.n_hits_/simulation::current_batch) * negroups_ + volume_sr * ray_tracing_cost_[base_sr];
    local_estimated_load += load_sr;
  }

  // Communicate estimated load across ranks
  MPI_Allgather(&local_estimated_load, 1, MPI_DOUBLE, estimated_rank_load_totals_.data(), 1, MPI_DOUBLE, mpi::intracomm);
  estimated_load_sum_ = std::accumulate(estimated_rank_load_totals_.begin(), estimated_rank_load_totals_.end(), 0.0);

  // Communicate measured load across ranks
  MPI_Allgather(&batch_transport_time, 1, MPI_DOUBLE, measured_rank_load_fractions_.data(), 1, MPI_DOUBLE, mpi::intracomm);
  double measured_load_sum = std::accumulate(measured_rank_load_fractions_.begin(), measured_rank_load_fractions_.end(), 0.0);

  // Calculate fractions //TODO: maybe I do not need fractions?
  for (int rank = 0; rank < mpi::n_procs; rank++) {
    estimated_rank_load_fractions_[rank] = estimated_rank_load_totals_[rank] / estimated_load_sum_;
    measured_rank_load_fractions_[rank] = measured_rank_load_fractions_[rank] / measured_load_sum;
  }

   // Reset batch-wise counters
  fill(num_base_source_region_RT_batch_.begin(), num_base_source_region_RT_batch_.end(), 0);
  fill(num_mesh_bin_RT_batch_.begin(), num_mesh_bin_RT_batch_.end(), 0);

  // Average measured rank load over previous batches
  vector<double> averaged_measured_load_fractions(mpi::n_procs, 0.0);
  int n_batches_to_average = load_history_size_;
  //TODO: A lot here can be simplified when load balancing is only done in first 5 batches
  if (simulation::current_batch < load_history_size_){
    n_batches_to_average = simulation::current_batch;
  }
  // Save to history
  for (int rank = 0; rank < mpi::n_procs; rank++) {
    load_history_[rank][history_idx] = measured_rank_load_fractions_[rank];
    averaged_measured_load_fractions[rank] = std::accumulate(load_history_[rank].begin(), load_history_[rank].end(), 0.0) / static_cast<double>(n_batches_to_average);
  }

  // Circular index update
  history_idx = (history_idx + 1) % load_history_size_;

  // Calculate measured imbalance
  double max_load_measured = *std::max_element(averaged_measured_load_fractions.begin(), averaged_measured_load_fractions.end());
  // TODO: Can be removed if load balancing fixed for first 5 batches only
  max_load_imbalance_measured_ = (max_load_measured - target_load_) / target_load_;
}

void DecompositionMap::balance_load(FlatSourceDomain* domain){

  //TODO: The optimisation strategy is messy
  cnt_optimizations_total_ ++;

  int max_iterations = 200;
  int it_outer = 0;
  double adaptation_factor = 1;
  double min_adaptation_factor = 0.01;
  double max_adaptation_factor = 2;
  double avg_rank_distance = max_domain_length_/cbrt(mpi::n_procs); // rough estimate of average distance between ranks
  double weight_scale = avg_rank_distance * avg_rank_distance * optimization_history_factor_;
  double beta = 0.6; // momentum damping
  bool check_all_ranks = true;

  vector<double> weight_change(mpi::n_procs, 0.0);
  vector<double> combined_rank_load(mpi::n_procs, 0.0);
  vector<double> load_ratio(mpi::n_procs, 0.0);

  // Combine estimated load with measured load ratios
  for (int rank = 0; rank < mpi::n_procs; rank++) {
    load_ratio[rank] = calculate_load_ratio(rank);
    combined_rank_load[rank] = load_ratio[rank] * estimated_rank_load_totals_[rank];
  }

  double combined_load_sum = std::accumulate(combined_rank_load.begin(), combined_rank_load.end(), 0.0);

  for (int rank = 0; rank < mpi::n_procs; rank++) {
    combined_rank_load[rank] = (combined_rank_load[rank]/combined_load_sum);
  }
  double max_load = *std::max_element(combined_rank_load.begin(), combined_rank_load.end());
  double max_imbalance = (max_load - target_load_) / target_load_;
  double prev_imbalance = max_imbalance;

  // History tracking
  vector<double> imbalance_history;
  vector<vector<double>> weight_history;
  imbalance_history.push_back(max_imbalance);
  weight_history.push_back(rank_weights_);

  // Change weights to equalize load based on combined load estimates
  while (max_imbalance > imbalance_tolerance_ && it_outer < max_iterations){
    // if (mpi::master)  printf("MPI load balancing iteration %d, max. imbalance: %.2f%% \n", it_outer, max_imbalance*100.0); //TODO: Remove
    it_outer ++;

    for (int rank = 0; rank < mpi::n_procs; rank++) {
      double corr = ((combined_rank_load[rank] - target_load_) / target_load_) * weight_scale;
      weight_change[rank] = beta * weight_change[rank] + (1.0 - beta) * corr; // keep some inertia from previous changes to prevent oscillations
      double damping = std::clamp(max_imbalance / prev_imbalance, 0.1, 1.0); // dampening factor based on whether we are getting closer to convergence or not, prevents big jumps, if too big a change, more dampening is applied
      rank_weights_[rank] -= adaptation_factor * damping * weight_change[rank];
    }

    if (simulation::current_batch > 1){
      // Check distances to all ranks only in first batch, otherwise only recorded neighbors
      check_all_ranks = false;
    }

    // Calculate new load after weight update
    update_load(domain, check_all_ranks, combined_rank_load, load_ratio);
    double max_load = *std::max_element(combined_rank_load.begin(), combined_rank_load.end());
    max_imbalance = (max_load - target_load_) / target_load_;

    // Store imbalance history
    imbalance_history.push_back(max_imbalance);
    weight_history.push_back(rank_weights_);

    // Adaptive factor
    if (max_imbalance > prev_imbalance)
      adaptation_factor = std::max(adaptation_factor * 0.5, min_adaptation_factor);
    else
      adaptation_factor = std::min(adaptation_factor * 1.05, max_adaptation_factor);

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
      // decrease weight for next batch if oscillating, TODO: Should this be safeguarded with separate min/max values?
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
      printf("Best imbalance during optimization was %.2f%% at iteration %d. \n", 
        best_imbalance*100.0, best_index);
    }

    //TODO: remove if load balancing only applied to first 5 batches by default
    if (cnt_unconverged_optimizations_ == 5){
      // relax tolerance if not converging
      cnt_unconverged_optimizations_ = 0;
      imbalance_tolerance_ = best_imbalance; 
      if (mpi::master){
        printf("Relaxing MPI load balancing tolerance to %.2f%% \n", imbalance_tolerance_*100.0);
      }
    }

    if (best_index == 0){
      // if no improvement at all, just keep current decomposition and return
      return;
    }
  } else {
    cnt_unconverged_optimizations_ = 0;
    optimization_history_factor_ = 1.0; // reset history factor if converged
    if (mpi::master){
      printf("MPI load balancing converged after %d iterations. Max. imbalance: %.2f%% \n", it_outer, max_imbalance*100.0);
    }
  }

  redistribute_source_regions(domain);
  MPI_Barrier(mpi::intracomm);
}

void DecompositionMap::update_load(FlatSourceDomain* domain, bool check_all_ranks, vector<double>& combined_rank_load, vector<double>& load_ratio){

  vector<double> load(mpi::n_procs, 0);

  // Add up source region hits to respective rank depending of position of centroid of each source region
  #pragma omp parallel
  {
    // number of hits per thread
    vector<double> thread_load(mpi::n_procs, 0);

    #pragma omp for
      for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
        Position centroid = domain->source_regions_.centroid(sr);
        int owner = find_closest_rank(centroid, check_all_ranks);
        double volume_sr = domain->source_regions_.volume_t(sr);
        thread_load[owner] += load_ratio[owner] * (C1_ * (domain->source_regions_.n_hits(sr)/simulation::current_batch) * negroups_ + volume_sr * ray_tracing_cost_[domain->source_regions_.key(sr).base_source_region_id]);
      }

    // Combine results from different threads
    #pragma omp critical (combining_loads)
    {
      for (int i = 0; i < mpi::n_procs; i++) {
        load[i] += thread_load[i];
      }
    }
  }

  // Communicate new load estimates across ranks
  MPI_Allreduce(MPI_IN_PLACE, load.data(), mpi::n_procs, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
  double load_sum = std::accumulate(load.begin(), load.end(), 0.0);

  // Update new combined rank load fractions
  for (int rank = 0; rank < mpi::n_procs; rank++){
    combined_rank_load[rank] = load[rank]/load_sum;
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

  // Each rank informs other ranks about ownership changes to update subdomain map
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
    }
  }

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

double DecompositionMap::calculate_load_ratio(int rank){
  if (estimated_rank_load_fractions_[rank] > 0.0){
    return measured_rank_load_fractions_[rank]/estimated_rank_load_fractions_[rank];
  } else {
    return 1.0;
  }
}


}// namespace openmc