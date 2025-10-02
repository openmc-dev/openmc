#include "openmc/random_ray/decomposition_map.h"

#include "openmc/message_passing.h"
#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/random_lcg.h"
#include "openmc/mgxs_interface.h"
#include "openmc/timer.h" //TODO: temporary?

#include "openmc/simulation.h"

namespace openmc {

namespace mpi {
    DecompositionMap decomp_map;
}

// Constructor
DecompositionMap::DecompositionMap() {
  // negroups_ = data::mg.num_energy_groups_; //TODO: constructor gets called before num_energy_groups is known?
  // rank_load_.resize(mpi::n_procs, 0.0);
}

void DecompositionMap::generate_rank_centers(){

  spatial_box_ = dynamic_cast<SpatialBox*>(
    dynamic_cast<IndependentSource*>(RandomRay::ray_source_.get())->space());

  // Calculate grid points that are used for Voronoi cells
  int grid_points_per_dimension = ceil(grid_points_per_rank_ * cbrt(mpi::n_procs)); // number of grid points in each dimension
  calculate_grid_points(grid_points_per_dimension);

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
    // if (mpi::master){
    //   printf("ITERATION %d \n", it);
    // }

    // Reset error to determine maximum movement below
    err = 0.0;

    vector<Position> position_sum_per_rank(mpi::n_procs, Position(0,0,0));
    vector<int> num_points_per_rank(mpi::n_procs, 0);

    // Compute Voronoi cells by summing up all position values of mesh grid points 
    // that are closest to a Voronoi center
    calculate_voronoi(position_sum_per_rank, num_points_per_rank);

    for (int rank = 0; rank < mpi::n_procs; rank++) {

      if (mpi::master){
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

  if (mpi::master){
    for (int rank = 0; rank < mpi::n_procs; rank++) {
      printf("  Rank %d: Point (%f, %f, %f)\n", rank, rank_centers_[rank].x, rank_centers_[rank].y, rank_centers_[rank].z);
      }
  }

  if (mpi::master) {
    if (it == max_iterations) {
      warning("Lloyd's algorithm did not converge within the maximum number of iterations.");
    } else {
      printf("Lloyd's algorithm converged in %d iterations.\n", it);
    }
  }
  
}

void DecompositionMap::calculate_grid_points(int grid_points_per_dimension){

    // Calculate grid spacing
    double delta_x = (spatial_box_->upper_right().x - spatial_box_->lower_left().x) / (grid_points_per_dimension - 1);
    double delta_y = (spatial_box_->upper_right().y - spatial_box_->lower_left().y) / (grid_points_per_dimension - 1);
    double delta_z = (spatial_box_->upper_right().z - spatial_box_->lower_left().z) / (grid_points_per_dimension - 1);
    
    // Generate all grid points
    
    for (int i = 0; i < grid_points_per_dimension; i++) {
        double x = spatial_box_->lower_left().x + i * delta_x;
        for (int j = 0; j < grid_points_per_dimension; j++) {
            double y = spatial_box_->lower_left().y + j * delta_y;
            for (int k = 0; k < grid_points_per_dimension; k++) {
                double z = spatial_box_->lower_left().z + k * delta_z;
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
    double z = prn(&seed);

    Position xi {x, y, z};

    // make a small shift in position to avoid geometry floating point issues //TODO: necessary?
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
  if (position_sum.norm() == 0.0){
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

// Update subdomain list for each decomposition map.
void DecompositionMap::update(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

    negroups_ = data::mg.num_energy_groups_; //TODO: This should not be called here again an dagain

    // outside_source_regions_.clear();

    // bool test  = discovered_source_regions.contains(SourceRegionKey(4, 104171));
    // printf("Exchange sr data info");

    bool any_sr_to_exchange = any_discovered_source_regions(discovered_source_regions);

    // Check if any new regions discovered
    // if (any_discovered_source_regions(discovered_source_regions)){
    if (any_sr_to_exchange){
      // Exchange discovered cell data between ranks
      simulation::time_source_region_exchange.start();
      exchange_sr_info(discovered_source_regions);
      simulation::time_source_region_exchange.stop();
    }

      // test  = discovered_source_regions.contains(SourceRegionKey(4, 104171));
      // printf("2.9 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);

  // vector<int> rank_counts(mpi::n_procs, 0);

  // for (auto [sr_key, rank] : discovered_regions_map_){
  //   // subdomain_map_[sr_key] = rank;
  // }

}

bool DecompositionMap::any_discovered_source_regions(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

  // int local_new_srs = discovered_regions_map_.size();
  int flag = 0;

  if(discovered_source_regions.begin() != discovered_source_regions.end()) {
  // if (local_new_srs > 0) {
    flag = 1;
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, mpi::intracomm);
  
  return flag > 0;
}

void DecompositionMap::exchange_sr_info(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

  // Communicate maps
  for (int rank = 0; rank < mpi::n_procs; rank++) {

    // if (simulation::current_batch == 6) {
    //   bool test  = discovered_source_regions.contains(SourceRegionKey(4, 104171));
    //   printf("2.2 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);
    // }

    // Send size
    uint64_t bcast_size = 0;
    if (rank == mpi::rank) {
      // bcast_size = discovered_source_regions.size();
      // bcast_size = discovered_regions_map_.size();
      for (const auto& pair : discovered_source_regions) { 
        if (pair.second.scalars_.volume_ > 0.0) { //TODO: can this check be avoided?
          bcast_size++;
        }
      }
    }

    MPI_Bcast(&bcast_size, 1, MPI_UINT64_T, rank, mpi::intracomm); //TODO: MPI_UINT64_T moght not be available?

    // printf("Rank %d: bcast_size %ld \n", mpi::rank, bcast_size);
    
    // if (simulation::current_batch == 6) {
    //   bool test  = discovered_source_regions.contains(SourceRegionKey(4, 104171));
    //   printf("2.3 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);
    // }

    if (bcast_size > 0) {

      // vector<int> local_rank_ids(bcast_size);
      vector<int64_t> local_base_ids(bcast_size);
      vector<int64_t> local_mesh_bins(bcast_size);
        
      if(rank == mpi::rank) {
        // fill in vectors to be sent
        int i = 0;
        // for (auto [sr_key, rank] : discovered_regions_map_){
        for (const auto& pair : discovered_source_regions) {
        // for (auto [sr_key, rank] : discovered_regions_map_){
          // local_rank_ids[i] = rank;
          if (pair.second.scalars_.volume_ > 0.0) {
            SourceRegionKey sr_key = pair.first;
            local_base_ids[i] = sr_key.base_source_region_id;
            local_mesh_bins[i] = sr_key.mesh_bin;
            i++;
          }
        }
      }

    // if (simulation::current_batch == 6) {
    //   bool test  = discovered_source_regions.contains(SourceRegionKey(4, 104171));
    //   printf("2.4 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);
    // }

      // printf("Rank %d: Before bcast \n", mpi::rank);

      // Broadcast all data
      // MPI_Bcast(local_rank_ids.data(), bcast_size, MPI_INT, rank, mpi::intracomm);
      MPI_Bcast(local_base_ids.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);
      MPI_Bcast(local_mesh_bins.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);

      // printf("Rank %d: After bcast \n", mpi::rank);
      
      // Update subdomain map
      for (int j = 0; j < bcast_size; j++) {


        // printf("Rank %d: Received source region with base id %ld, mesh bin %ld from rank %d \n", mpi::rank, local_base_ids[j], local_mesh_bins[j], rank);
        // printf("Source region with base id %ld, mesh bin %ld received from rank %d \n", local_base_ids[j], local_mesh_bins[j], rank);

      //   if (local_base_ids[j] == 4 && local_mesh_bins[j] == 104171){
      //     printf("Received discovered source region with base id %ld, mesh bin %ld from rank %d \n", local_base_ids[j], local_mesh_bins[j], rank);
      // bool test  = discovered_source_regions.contains(SourceRegionKey(4, 104171));
      // printf("2.5 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);
      //   }
        SourceRegionKey sr_key(local_base_ids[j], local_mesh_bins[j]);

        // check if already in map or not, i.e. has someone else discovered that region already?
        if (subdomain_map_.find(sr_key) == subdomain_map_.end()){
      //     if (local_base_ids[j] == 4 && local_mesh_bins[j] == 104171){
            // printf("Writing source region with base id %ld, mesh bin %ld to map, assigned to rank %d \n", local_base_ids[j], local_mesh_bins[j], rank);
      // bool test  = discovered_source_regions.contains(SourceRegionKey(4, 104171));
      // printf("2.6 RANK %d: Source region (4, 104171) discovered? %d \n", mpi::rank, test);
      //     }
          // subdomain_map_[sr_key] = local_rank_ids[j];

          subdomain_map_[sr_key] = rank;
          
        } else{
          // printf("Source region with base id %ld, mesh bin %ld already in map, contested between rank %d and rank %d \n", local_base_ids[j], local_mesh_bins[j], rank, subdomain_map_[sr_key]);

          int resident_rank = subdomain_map_[sr_key];
          int challenger_rank = rank; // current broadcasting rank
          // int challenger_rank = local_rank_ids[j];
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

          // if (mpi::master){
            // printf("sender: %d, receiver: %d, resident: %d, challenger: %d \n", sender, receiver, challenger_rank, resident_rank);
            // printf("resident_rank_load: %f, challenger_rank_load: %f \n", resident_rank_load, challenger_rank_load);
          // }

          // Broadcast n_hits such that each rank can update load balance //TODO: n_hits must be batchwise
          // TODO: is load update here too costly?
          int bcast_n_hits = 0;
          if (mpi::rank == sender) {
            SourceRegion& contested_sr = discovered_source_regions[sr_key];
            bcast_n_hits = contested_sr.scalars_.n_hits_;
          }
          MPI_Bcast(&bcast_n_hits, 1, MPI_INT, sender, mpi::intracomm);

          rank_load_[sender] -= bcast_n_hits;
          rank_load_[receiver] += bcast_n_hits;

          // Merge source region data
          if (mpi::rank == sender || mpi::rank == receiver){
            transfer_sr_data(sender, receiver, sr_key, discovered_source_regions);
          }

          // if (mpi::master){
          //   printf(" Source regions transferred\n");
          // }
        }
      }
    }
  }

}

void DecompositionMap::transfer_sr_data(int sender, int receiver, SourceRegionKey sr_key, ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

  bool is_linear = RandomRay::source_shape_ != RandomRaySourceShape::FLAT;

  SourceRegion& sr = discovered_source_regions[sr_key]; //TODO: maybe this can be done more efficiently?

  int num_scalar_messages = 1;
  //! NOTE: update if new vector fields are added to SourceExchangeVectors struct
  int num_vector_messages = 4;
  if (is_linear){  
    num_vector_messages += 4;
  }
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    num_vector_messages += 1;
  }
  int num_requests = num_scalar_messages + num_vector_messages;

  vector<MPI_Request> requests(num_requests);
  // MPI_Request requests[num_requests];
  int req_idx = 0;

  if (mpi::rank == sender){ //TODO: maybe add check that only sends fields that are needed when source regions are freshly discovered?

    // SourceRegion sr = discovered_source_regions[sr_key]; //TODO: maybe this can be done more efficiently?

    // printf("Rank %d sending source region with base id %ld, mesh bin %ld to rank %d \n", sender, sr_key.base_source_region_id, sr_key.mesh_bin, receiver);

    // printf("RANK %d: Source Region exists %d \n", mpi::rank, discovered_source_regions.contains(sr_key));
    // printf("Size 1 flux new on rank %d: %lu \n", mpi::rank, sr.scalar_flux_old_.size());
    // printf("Size 2 flux new on rank %d: %lu \n", mpi::rank, sr.scalar_flux_new_.size());
    // printf("Size 3 flux new on rank %d: %lu \n", mpi::rank, sr.source_.size());
    // printf("Size 4 flux new on rank %d: %lu \n", mpi::rank, sr.external_source_.size());
    // printf("Size 5 flux new on rank %d: %lu \n", mpi::rank, sr.scalar_flux_final_.size());
    // printf("is linear: %d \n", is_linear);

    // Send scalar data to receiver
    MPI_Isend(&sr.scalars_, sizeof(ScalarSourceRegionFields), MPI_BYTE, receiver, 1, mpi::intracomm, &requests[req_idx]);
    req_idx ++;

    // Send vector data to receiver
    // Tags hardcoded to avoid confusion if new fields are not added sequentially
    MPI_Isend(sr.scalar_flux_old_.data(), negroups_, MPI_DOUBLE, receiver, 2, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;

    // printf("Size 1 flux new on rank %d: %lu \n", mpi::rank, sr.scalar_flux_old_.size());

    MPI_Isend(sr.scalar_flux_new_.data(), negroups_, MPI_DOUBLE, receiver, 3, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;

    // printf("Size 2 flux new on rank %d: %lu \n", mpi::rank, sr.scalar_flux_new_.size());

    MPI_Isend(sr.source_.data(), negroups_, MPI_FLOAT, receiver, 4, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;
    // printf("Size 3 flux new on rank %d: %lu \n", mpi::rank, sr.source_.size());

    if (settings::run_mode == RunMode::FIXED_SOURCE) {
        MPI_Isend(sr.external_source_.data(), negroups_, MPI_FLOAT, receiver, 5, mpi::intracomm, &requests[req_idx]); 
        req_idx ++;
        // printf("Size 4 flux new on rank %d: %lu \n", mpi::rank, sr.external_source_.size());
    }

    MPI_Isend(sr.scalar_flux_final_.data(), negroups_, MPI_DOUBLE, receiver, 6, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;
    // printf("Size 5 flux new on rank %d: %lu \n", mpi::rank, sr.scalar_flux_final_.size());

    if (is_linear){
      MPI_Isend(sr.source_gradients_.data(), 3*negroups_, MPI_DOUBLE, receiver, 7, mpi::intracomm, &requests[req_idx]); 
      req_idx ++;

      MPI_Isend(sr.flux_moments_old_.data(), 3*negroups_, MPI_DOUBLE, receiver, 8, mpi::intracomm, &requests[req_idx]); 
      req_idx ++;

      MPI_Isend(sr.flux_moments_new_.data(), 3*negroups_, MPI_DOUBLE, receiver, 9, mpi::intracomm, &requests[req_idx]); 
      req_idx ++;

      MPI_Isend(sr.flux_moments_t_.data(), 3*negroups_, MPI_DOUBLE, receiver, 10, mpi::intracomm, &requests[req_idx]); 
      req_idx ++;
    }
    // printf("Send complete");

    if (req_idx != num_requests){
      fatal_error(fmt::format("Number of MPI requests does not match number of messages sent."
        "Check if num_vector_messages is up to date in DecompositionMap::transfer_sr_data()."));
    }

    // Wait for all communication to complete
    // MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE);
    MPI_Waitall(num_requests, requests.data(), MPI_STATUSES_IGNORE);

    discovered_source_regions.erase(sr_key);
    // printf("Source region erased from rank %d \n", sender);

  } else if (mpi::rank == receiver){

    // printf("RANK %d: Source Region exists %d \n", mpi::rank, discovered_source_regions.contains(sr_key));

    // printf("Energy groups: %d \n", negroups_);

    SourceRegion sr_recv(negroups_, is_linear);

    // sr_recv.scalar_flux_old_.resize(negroups_);
    // sr_recv.scalar_flux_new_.resize(negroups_);
    // sr_recv.source_.resize(negroups_);
    // sr_recv.external_source_.resize(negroups_);
    // sr_recv.scalar_flux_final_.resize(negroups_);
    // if (is_linear){
    //   sr_recv.source_gradients_.resize(negroups_);
    //   sr_recv.flux_moments_old_.resize(negroups_);
    //   sr_recv.flux_moments_new_.resize(negroups_);
    //   sr_recv.flux_moments_t_.resize(negroups_);
    // }

    // printf("Size 1 flux new on rank %d: %lu \n", mpi::rank, sr_recv.scalar_flux_new_.size());

    // Receive scalar data from sender
    MPI_Recv(&sr_recv.scalars_, sizeof(ScalarSourceRegionFields), MPI_BYTE, sender, 1, mpi::intracomm, MPI_STATUS_IGNORE);

    // Receive vector data from sender
    MPI_Recv(sr_recv.scalar_flux_old_.data(), negroups_, MPI_DOUBLE, sender, 2, mpi::intracomm, MPI_STATUS_IGNORE); 
    MPI_Recv(sr_recv.scalar_flux_new_.data(), negroups_, MPI_DOUBLE, sender, 3, mpi::intracomm, MPI_STATUS_IGNORE); 
    
    // printf("Size 2 flux new on rank %d: %lu \n", mpi::rank, sr_recv.scalar_flux_new_.size());

    
    MPI_Recv(sr_recv.source_.data(), negroups_, MPI_FLOAT, sender, 4, mpi::intracomm, MPI_STATUS_IGNORE); 

    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      MPI_Recv(sr_recv.external_source_.data(), negroups_, MPI_FLOAT, sender, 5, mpi::intracomm, MPI_STATUS_IGNORE); 
    }

    MPI_Recv(sr_recv.scalar_flux_final_.data(), negroups_, MPI_DOUBLE, sender, 6, mpi::intracomm, MPI_STATUS_IGNORE); 

    if (is_linear){
      MPI_Recv(sr_recv.source_gradients_.data(), 3*negroups_, MPI_DOUBLE, sender, 7, mpi::intracomm, MPI_STATUS_IGNORE); 
      MPI_Recv(sr_recv.flux_moments_old_.data(), 3*negroups_, MPI_DOUBLE, sender, 8, mpi::intracomm, MPI_STATUS_IGNORE); 
      MPI_Recv(sr_recv.flux_moments_new_.data(), 3*negroups_, MPI_DOUBLE, sender, 9, mpi::intracomm, MPI_STATUS_IGNORE); 
      MPI_Recv(sr_recv.flux_moments_t_.data(), 3*negroups_, MPI_DOUBLE, sender, 10, mpi::intracomm, MPI_STATUS_IGNORE); 
    }
    // // Merge source region data %TODO: is it OK to make this method of source region class?
    // sr_recv.merge(sr, is_linear);

    // // clear old source region data from discovered regions map //TODO: deleting and replacing seems wasteful!
    // discovered_source_regions.erase(sr_key);

    // printf("Source region received on rank %d \n", receiver);
    // printf("RANK %d: Source Region exists %d \n", mpi::rank, discovered_source_regions.contains(sr_key));
  
    // Merge source region data %TODO: is it OK to make this method of source region class?
    sr.merge(sr_recv, is_linear);
    // TODO: will this feed back automatically into discovered_source_regions?

    // // add new merged source region to discovered regions map
    // discovered_source_regions.emplace(sr_key, sr_recv);

  }

  // // Wait for all communication to complete
  // MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE);

  // if (mpi::rank == sender){
  //   // clear old source region data from discovered regions map
  //   discovered_source_regions.erase(sr_key);
  //   printf("Source region erased from rank %d \n", sender);
  // }

}

int DecompositionMap::find_owner(SourceRegionKey sr_key, Position r, 
    ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions, int64_t ID){ //TODO: Remove ID
  
  // Check if SRK is in rank's subdomain
  auto it = subdomain_map_.find(sr_key);
  if (it != subdomain_map_.end()){
    // if (sr_key.base_source_region_id == 134 && sr_key.mesh_bin == 64056 || ID == 1083){
    //   printf("RANK %d: Identified owner of ray %ld at check 1 as rank %d\n", mpi::rank, ID, it->second);
    //   printf("RANK %d: Position (%f, %f, %f) \n", mpi::rank, r.x, r.y, r.z);
    // }
    return it->second;
  } 

  // Check if already recorded in newly discovered source regions
  discovered_source_regions.lock(sr_key);
  bool sr_key_discovered = discovered_source_regions.contains(sr_key);
  discovered_source_regions.unlock(sr_key);
  if (sr_key_discovered){
    // if(sr_key.base_source_region_id == 134 && sr_key.mesh_bin == 64056 || ID == 1083){
    //   printf("RANK %d: Identified owner of ray %ld at check 2 as rank %d\n", mpi::rank, ID, mpi::rank);
    //   printf("RANK %d: Position (%f, %f, %f) \n", mpi::rank, r.x, r.y, r.z);
    // }
    // printf("RANK %d: Identified owner at check 2 as rank %d\n", mpi::rank, mpi::rank);
    return mpi::rank;
  }

  // If not found in either map, check if my rank would own source 
  // region beased on location
  int closest_rank = find_closest_rank(r);
  // if (sr_key.base_source_region_id == 134 && sr_key.mesh_bin == 64056 || ID == 1083){
  //   printf("RANK %d: Identified owner of ray %ld at check 4 as rank %d\n", mpi::rank, ID, closest_rank);
  //   printf("RANK %d: Position (%f, %f, %f) \n", mpi::rank, r.x, r.y, r.z);
  // }
  // printf("RANK %d: Identified owner at check 4 as rank %d\n", mpi::rank, closest_rank);
  return closest_rank;
}

// bool DecompositionMap::is_closest_rank(Position r) {
int DecompositionMap::find_closest_rank(Position r) {
  // Determine which rank the position belongs to
  int closest_rank = C_NONE;
  double min_distance = INFTY;
  
  // Find closest rank center
  for (int rank = 0; rank < mpi::n_procs; rank++) {
      double dist = (r - rank_centers_[rank]).norm();
      if (dist < min_distance) {
          min_distance = dist;
          closest_rank = rank;
      }
  }

  if (mpi::master && closest_rank == C_NONE) {
    fatal_error("Could not find closest rank for new source region at position (" 
      + std::to_string(r.x) + ", " + std::to_string(r.y) + ", " + std::to_string(r.z) + ").");
  }

  // check whether it belongs to this rank's subdomain
  // bool is_in_my_subdomain = (closest_rank == mpi::rank);
  // return is_in_my_subdomain;
  return closest_rank;
}

// void DecompositionMap::calculate_rank_load(uint64_t total_geometric_intersections){
  
//   // Reset rank load
//   std::fill(rank_load_.begin(), rank_load_.end(), 0.0);

//   // Temporary print-outs
//   if (mpi::master) {
//     printf("Load distribution: ");
//   }

//   // Calculate rank load based on number of geometric intersections
//   uint64_t total_geometric_intersections_sum = 0;
//   MPI_Allreduce(&total_geometric_intersections, &total_geometric_intersections_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi::intracomm);
//   double load = (total_geometric_intersections/static_cast<double>(total_geometric_intersections_sum))*100.0;

//   for (int rank = 0; rank < mpi::n_procs; rank++) {
//     // Broadcast load
//     double bcast_load = 0.0;
//     if (rank == mpi::rank) {
//       bcast_load = load;
//     }
//     MPI_Bcast(&bcast_load, 1, MPI_DOUBLE, rank, mpi::intracomm);
//     rank_load_[rank] = bcast_load;
//     // Temporary print-outs
//     if (mpi::master) {
//       printf("RANK %d: %.2f%% ", rank, rank_load_[rank]);
//     }
//   }
//   // Temporary print-outs
//   if (mpi::master) {
//     printf("\n");
//   }
// }

void DecompositionMap::calculate_rank_load(FlatSourceDomain* domain){ //TODO: This uses accumulated value, not batch-wise values!
  
  rank_load_.resize(mpi::n_procs, 0.0); //TODO: that should be moved elsewhere

  // Reset rank load
  std::fill(rank_load_.begin(), rank_load_.end(), 0.0);

  // Add number of hits for each source region by going through all exisiting source regions.
  int64_t n_hits_rank = 0;
  
  // add source region hits
  #pragma omp parallel for reduction(+ : n_hits_rank)
    for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
      n_hits_rank += domain->source_regions_.n_hits(sr);
    }

  // add newly discovered source region hits
  for (const auto & [sr_key, sr] : domain->discovered_source_regions_) {
    n_hits_rank += sr.scalars_.n_hits_;
  }

  // Temporary print-outs
  if (mpi::master) {
    printf("Load distribution: ");
  }

  // Calculate rank load based on number of geometric intersections
  uint64_t n_hits_total = 0;
  MPI_Allreduce(&n_hits_rank, &n_hits_total, 1, MPI_INT64_T, MPI_SUM, mpi::intracomm);
  double load = (n_hits_rank/static_cast<double>(n_hits_total))*100.0;
  
  if(mpi::master){
    printf("Total hits: %ld \n", n_hits_total);
  }

  for (int rank = 0; rank < mpi::n_procs; rank++) {
    // Broadcast load
    double bcast_load = 0.0;
    if (rank == mpi::rank) {
      bcast_load = load;
    }
    MPI_Bcast(&bcast_load, 1, MPI_DOUBLE, rank, mpi::intracomm);
    rank_load_[rank] = bcast_load;
    // Temporary print-outs
    if (mpi::master) {
      printf("RANK %d: %.2f%% ", rank, rank_load_[rank]);
    }
  }
  // Temporary print-outs
  if (mpi::master) {
    printf("\n");
  }
}

}// namespace openmc