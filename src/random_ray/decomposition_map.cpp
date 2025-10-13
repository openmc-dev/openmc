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
DecompositionMap::DecompositionMap() {}

void DecompositionMap::initialize(){
  negroups_ = data::mg.num_energy_groups_;
  rank_load_.resize(mpi::n_procs, 0.0);
  target_load_ = 1.0/mpi::n_procs;
  // rank_weights_.resize(mpi::n_procs, 0.0);
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
  // int grid_points_per_dimension = ceil(grid_points_per_rank_ * cbrt(mpi::n_procs)); // number of grid points in each dimension
  int grid_points_total = grid_points_per_rank_ * mpi::n_procs;
  printf("Calculating %d grid points for Voronoi tessellation...\n", grid_points_total);
  // printf("Grid points per dimension: %d \n", grid_points_per_dimension);
  // calculate_grid_points(grid_points_per_dimension);
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

// void DecompositionMap::calculate_grid_points(int grid_points_per_dimension){
void DecompositionMap::calculate_grid_points(int grid_points_total){

    // Calculate grid spacing     
    vector <double> domain_length(3);
    domain_length[0] = spatial_box_->upper_right().x - spatial_box_->lower_left().x;
    domain_length[1] = spatial_box_->upper_right().y - spatial_box_->lower_left().y;
    domain_length[2] = spatial_box_->upper_right().z - spatial_box_->lower_left().z;
    // double x_length = spatial_box_->upper_right().x - spatial_box_->lower_left().x;
    // double y_length = spatial_box_->upper_right().y - spatial_box_->lower_left().y;
    // double z_length = spatial_box_->upper_right().z - spatial_box_->lower_left().z;

    double volume = domain_length[0] * domain_length[1] * domain_length[2];

    // printf("Spatial box dimensions: %f x %f x %f \n", x_length, y_length, z_length);
    printf("Spatial box volume: %f \n", volume);
    // printf("Grid points per dimension: %d \n", grid_points_per_dimension);

    // int grid_points_x = std::max(2, static_cast<int>(round((x_length / cbrt(volume)) * cbrt(grid_points_total))));
    // int grid_points_y = std::max(2, static_cast<int>(round((y_length / cbrt(volume)) * cbrt(grid_points_total))));
    // int grid_points_z = std::max(2, static_cast<int>(round((z_length / cbrt(volume)) * cbrt(grid_points_total))));

    int excluded_dimension = -1;
    vector<int> grid_points_per_dimension(3);
    for (int i = 0; i < 3; i++){
      double grid_points_estimate = ((domain_length[i] / cbrt(volume)) * cbrt(grid_points_total));
      // Check if any dimension distorted
      if (grid_points_estimate > 2){
        grid_points_per_dimension[i] = std::max(2, static_cast<int>(round(grid_points_estimate)));
      } else {
        excluded_dimension = i;
        grid_points_per_dimension[i] = 2;
      }
    }

    if (excluded_dimension != -1){
      double area = 1.0;
      
      for (int i = 0; i < 3; i++){
        if (i == excluded_dimension) continue;
        area *= domain_length[i];
      }

      // Recalculate grid points in non-distorted dimensions based on area
      for (int i = 0; i < 3; i++){
        if (i == excluded_dimension) continue;
        grid_points_per_dimension[i] = round((domain_length[i] / sqrt(area)) * sqrt(grid_points_total/2.0));
      }
    } 

    printf("Grid points in x, y, z: %d, %d, %d \n", grid_points_per_dimension[0], grid_points_per_dimension[1], grid_points_per_dimension[2]);

    // int total_grid_points = grid_points_per_dimension * grid_points_per_dimension * grid_points_per_dimension;

    // int grid_points_x = std::max(2, static_cast<int>(round((x_length / cbrt(volume)) * grid_points_per_dimension)));
    // int grid_points_y = std::max(2, static_cast<int>(round((y_length / cbrt(volume)) * grid_points_per_dimension)));
    // int grid_points_z = std::max(2, static_cast<int>(round((z_length / cbrt(volume)) * grid_points_per_dimension)));

    double new_total = grid_points_per_dimension[0] * grid_points_per_dimension[1] * grid_points_per_dimension[2];
    double adjustment;
    if (excluded_dimension != -1) {
      // When one dimension is fixed at 2, use square root for the other two dimensions
      adjustment = sqrt(grid_points_total / new_total);
    } else {
      // When all dimensions are free, use cube root
      adjustment = cbrt(grid_points_total / new_total);
    }

    // double adjustment = cbrt(grid_points_total / new_total);

    for (int i = 0; i < 3; i++){
      if (i == excluded_dimension) continue;
      grid_points_per_dimension[i] = round(grid_points_per_dimension[i] * adjustment);
    }
    
    // grid_points_x = std::max(2, static_cast<int>(round(grid_points_x * adjustment)));
    // grid_points_y = std::max(2, static_cast<int>(round(grid_points_y * adjustment)));
    // grid_points_z = std::max(2, static_cast<int>(round(grid_points_z * adjustment))); 

    // double exponent_x = (x_length / cbrt(volume));
    // double exponent_y = (y_length / cbrt(volume));
    // double exponent_z = (z_length / cbrt(volume));

    // int grid_points_x = std::pow(total_grid_points, exponent_x/3.0);
    // int grid_points_y = std::pow(total_grid_points, exponent_y/3.0);
    // int grid_points_z = std::pow(total_grid_points, exponent_z/3.0);



    printf("Grid points in x, y, z: %d, %d, %d \n", grid_points_per_dimension[0], grid_points_per_dimension[1], grid_points_per_dimension[2]);

    double delta_x = (spatial_box_->upper_right().x - spatial_box_->lower_left().x) / (grid_points_per_dimension[0] - 1);
    double delta_y = (spatial_box_->upper_right().y - spatial_box_->lower_left().y) / (grid_points_per_dimension[1] - 1);
    double delta_z = (spatial_box_->upper_right().z - spatial_box_->lower_left().z) / (grid_points_per_dimension[2] - 1);

    // double delta_x = (spatial_box_->upper_right().x - spatial_box_->lower_left().x) / (grid_points_per_dimension - 1);
    // double delta_y = (spatial_box_->upper_right().y - spatial_box_->lower_left().y) / (grid_points_per_dimension - 1);
    // double delta_z = (spatial_box_->upper_right().z - spatial_box_->lower_left().z) / (grid_points_per_dimension - 1);
    
    // Generate all grid points
    
    for (int i = 0; i < grid_points_per_dimension[0]; i++) {
        double x = spatial_box_->lower_left().x + i * delta_x;
        for (int j = 0; j < grid_points_per_dimension[1]; j++) {
            double y = spatial_box_->lower_left().y + j * delta_y;
            for (int k = 0; k < grid_points_per_dimension[2]; k++) {
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

  // check if any points have been recorded in rank
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

// Update subdomain list for each decomposition map.
void DecompositionMap::update(ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions){

    // negroups_ = data::mg.num_energy_groups_; //TODO: This should not be called here again an dagain

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

  // bool is_linear = RandomRay::source_shape_ != RandomRaySourceShape::FLAT; 

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

          double load_change = bcast_n_hits / static_cast<double>(n_hits_sum_);
          rank_load_[sender] -= load_change;
          rank_load_[receiver] += load_change;

          // Merge source region data
          // if (mpi::rank == sender || mpi::rank == receiver){
          //   SourceRegion& sr = discovered_source_regions[sr_key];
          //   SourceRegion sr_recv(negroups_, is_linear);
          //   transfer_sr_data(sender, receiver, sr, sr_recv, is_linear);

            if (mpi::rank == sender){
              SourceRegion& sr = discovered_source_regions[sr_key];
              // SourceRegion sr_recv(negroups_, is_linear);
              // transfer_sr_data(sender, receiver, sr, sr_recv, is_linear);
              send_sr_data(receiver, sr);

              // clear old source region data from discovered regions map
              discovered_source_regions.erase(sr_key);
              // printf("Source region erased from rank %d \n", sender);
            }
            if (mpi::rank == receiver){
              SourceRegion& sr = discovered_source_regions[sr_key];
              SourceRegion sr_recv(negroups_, is_linear_);
              // transfer_sr_data(sender, receiver, sr, sr_recv, is_linear);
              receive_sr_data(sender, sr_recv);

              sr.merge(sr_recv, is_linear_);
            }
            
          // }

          // if (mpi::master){
          //   printf(" Source regions transferred\n");
          // }
        }
      }
    }
  }

}

// void DecompositionMap::transfer_sr_data(int sender, int receiver, SourceRegionKey sr_key, ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
//     discovered_source_regions){

void DecompositionMap::transfer_sr_data(int sender, int receiver, SourceRegion& sr_send, SourceRegion& sr_recv){

  // bool is_linear = RandomRay::source_shape_ != RandomRaySourceShape::FLAT;

  // SourceRegion& sr_send = discovered_source_regions[sr_key];

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

  if (mpi::rank == sender){ //TODO: maybe add check that only sends fields that are needed when source regions are freshly discovered?

    // Send scalar data to receiver
    MPI_Isend(&sr_send.scalars_, sizeof(ScalarSourceRegionFields), MPI_BYTE, receiver, 1, mpi::intracomm, &requests[req_idx]);
    req_idx ++;

    // Send vector data to receiver
    // Tags hardcoded to avoid confusion if new fields are not added sequentially
    MPI_Isend(sr_send.scalar_flux_old_.data(), negroups_, MPI_DOUBLE, receiver, 2, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;

    // printf("Size 1 flux new on rank %d: %lu \n", mpi::rank, sr_send.scalar_flux_old_.size());

    MPI_Isend(sr_send.scalar_flux_new_.data(), negroups_, MPI_DOUBLE, receiver, 3, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;

    // printf("Size 2 flux new on rank %d: %lu \n", mpi::rank, sr_send.scalar_flux_new_.size());

    MPI_Isend(sr_send.source_.data(), negroups_, MPI_FLOAT, receiver, 4, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;
    // printf("Size 3 flux new on rank %d: %lu \n", mpi::rank, sr_send.source_.size());

    if (settings::run_mode == RunMode::FIXED_SOURCE) {
        MPI_Isend(sr_send.external_source_.data(), negroups_, MPI_FLOAT, receiver, 5, mpi::intracomm, &requests[req_idx]); 
        req_idx ++;
        // printf("Size 4 flux new on rank %d: %lu \n", mpi::rank, sr_send.external_source_.size());
    }

    MPI_Isend(sr_send.scalar_flux_final_.data(), negroups_, MPI_DOUBLE, receiver, 6, mpi::intracomm, &requests[req_idx]); 
    req_idx ++;
    // printf("Size 5 flux new on rank %d: %lu \n", mpi::rank, sr_send.scalar_flux_final_.size());

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
    // printf("Send complete");

    if (req_idx != num_requests){
      fatal_error(fmt::format("Number of MPI requests does not match number of messages sent."
        "Check if num_vector_messages is up to date in DecompositionMap::transfer_sr_data()."));
    }

    // Wait for all communication to complete
    // MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE);
    MPI_Waitall(num_requests, requests.data(), MPI_STATUSES_IGNORE);

    // discovered_source_regions.erase(sr_key);
    // printf("Source region erased from rank %d \n", sender);

  } else if (mpi::rank == receiver){

    // printf("RANK %d: Source Region exists %d \n", mpi::rank, discovered_source_regions.contains(sr_key));

    // printf("Energy groups: %d \n", negroups_);

    // SourceRegion sr_recv(negroups_, is_linear);

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
  
    // Merge source region data %TODO: is it OK to make this method of source region class?
    // sr_send.merge(sr_recv, is_linear);
    // TODO: will this feed back automatically into discovered_source_regions?
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

  // printf("Size 1 flux new on rank %d: %lu \n", mpi::rank, sr_send.scalar_flux_old_.size());

  MPI_Isend(sr_send.scalar_flux_new_.data(), negroups_, MPI_DOUBLE, receiver, 3, mpi::intracomm, &requests[req_idx]); 
  req_idx ++;

  // printf("Size 2 flux new on rank %d: %lu \n", mpi::rank, sr_send.scalar_flux_new_.size());

  MPI_Isend(sr_send.source_.data(), negroups_, MPI_FLOAT, receiver, 4, mpi::intracomm, &requests[req_idx]); 
  req_idx ++;
  // printf("Size 3 flux new on rank %d: %lu \n", mpi::rank, sr_send.source_.size());

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
      MPI_Isend(sr_send.external_source_.data(), negroups_, MPI_FLOAT, receiver, 5, mpi::intracomm, &requests[req_idx]); 
      req_idx ++;
      // printf("Size 4 flux new on rank %d: %lu \n", mpi::rank, sr_send.external_source_.size());
  }

  MPI_Isend(sr_send.scalar_flux_final_.data(), negroups_, MPI_DOUBLE, receiver, 6, mpi::intracomm, &requests[req_idx]); 
  req_idx ++;
  // printf("Size 5 flux new on rank %d: %lu \n", mpi::rank, sr_send.scalar_flux_final_.size());

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
  // printf("Send complete");

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
      // Power Voronoi diagram uses sqaured distances
      // dist *= dist;
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
  
  // rank_load_.resize(mpi::n_procs, 0.0); //TODO: that should be moved elsewhere

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
  double load = (n_hits_rank/static_cast<double>(n_hits_total));
  n_hits_sum_ = n_hits_total;
  
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
      printf("RANK %d: %.2f%% ", rank, rank_load_[rank]*100.0);
    }
  }
  // Temporary print-outs
  if (mpi::master) {
    printf("\n");
  }

  double max_load = *std::max_element(rank_load_.begin(), rank_load_.end());
  double avg_load = std::accumulate(rank_load_.begin(), rank_load_.end(), 0.0)/mpi::n_procs;
  max_imbalance_ = (max_load - avg_load) / avg_load;
}

void DecompositionMap::balance_load(FlatSourceDomain* domain){

  // if (mpi::master) printf("Starting MPI load balancing... \n");

  int max_iterations = 1000;
  // double max_load = *std::max_element(rank_load_.begin(), rank_load_.end());
  // double avg_load = std::accumulate(rank_load_.begin(), rank_load_.end(), 0.0)/mpi::n_procs;
  // double max_imbalance = (max_load - avg_load) / avg_load;
  double max_imbalance = max_imbalance_;

  // double tol_outer = 0.01;
  int it_outer = 0;
  double adaptation_factor = 1; //0.5 //0.1;
  double min_factor = 0.01;
  double max_factor = 2;
  double prev_imbalance = max_imbalance;
  double avg_rank_distance = max_domain_length_/cbrt(mpi::n_procs); // rough estimate of average distance between ranks //TODO: cubic root might not be adaquate for problems that are not box like
  double weight_scale = avg_rank_distance * avg_rank_distance; //TODO: maybe adjust dynamically dependent on wether previous batch has converged
  // double weight_scale = std::accumulate(rank_weights_.begin(), rank_weights_.end(), 0.0) / mpi::n_procs; // scale weights based on average load
  double beta = 0.6; // momentum damping

  vector<double> weight_change(mpi::n_procs, 0.0);

  vector<double> old_weights = rank_weights_;
  double max_imbalance_old = max_imbalance_;

  std::unordered_map<int, vector<int>> sr_send; // contains regions to be send

  // if (mpi::master) printf("Max. load: %.2f%%, Avg. load: %.2f%%, Max. imbalance: %.2f%% \n", max_load*100.0, avg_load*100.0, max_imbalance*100.0);

  while (max_imbalance > imbalance_tolerance_ && it_outer < max_iterations){

    if (mpi::master)  printf("MPI load balancing iteration %d, max. imbalance: %.2f%% \n", it_outer, max_imbalance*100.0);

    // Optimise rank load rank by rank
    // for (int rank = 0; rank < mpi::n_procs; rank++) {
    //   double correction = adaptation_factor * ((rank_load_[rank] - target_load_) / target_load_); // * weight_scale;
    //   rank_weights_[rank] -= correction;
    // }

    for (int rank = 0; rank < mpi::n_procs; rank++) {
      // weight_scale = std::accumulate(rank_weights_.begin(), rank_weights_.end(), 0.0) / mpi::n_procs;
      double corr = ((rank_load_[rank] - target_load_) / target_load_) * weight_scale;
      weight_change[rank] = beta * weight_change[rank] + (1.0 - beta) * corr; // keep some inertia from previous changes to prevent oscillations
      // double damping = std::min(1.0, max_imbalance / imbalance_tolerance_); // dampening factor based on how far we are from convergence
      double damping = std::clamp(max_imbalance / prev_imbalance, 0.1, 1.0); // dampening factor based on whether we are getting closer to convergence or not, prevents big jumps, if too big a change, more dampening is applied
      rank_weights_[rank] -= adaptation_factor * damping * weight_change[rank];
    }

    update_load(domain);
    it_outer ++;

    double max_load = *std::max_element(rank_load_.begin(), rank_load_.end());
    double avg_load = std::accumulate(rank_load_.begin(), rank_load_.end(), 0.0) / mpi::n_procs;
    max_imbalance = (max_load - avg_load) / avg_load;

    // Adaptive factor
    if (max_imbalance > prev_imbalance)
        adaptation_factor = std::max(adaptation_factor * 0.5, min_factor); // don’t go below min
    else
        adaptation_factor = std::min(adaptation_factor * 1.05, max_factor);     // stable - accelerate slightly

    prev_imbalance = max_imbalance;

  }

  if (it_outer == max_iterations){
    if(mpi::master) {
      warning("MPI load balancing has not converged after "
              + std::to_string(max_iterations) + " iterations."); 
    }
    if (max_imbalance_old < max_imbalance){
      rank_weights_ = old_weights;
      return;
    }
  } else {
    if (mpi::master){
      printf("MPI load balancing converged after %d iterations. Max. imbalance: %.2f%% \n", it_outer, max_imbalance*100.0);
    }
  }

  // if (mpi::master)  printf("Finished rebalancing \n");

  // fatal_error("STOP");

  // update decomposition map
  redistribute_source_regions(domain);


  if (mpi::master){
    printf("Weights: ");
    for (int rank = 0; rank < mpi::n_procs; rank++) {
      printf("RANK %d: %.2f ", rank, rank_weights_[rank]);
    }
    printf("\n");
  }

}

void DecompositionMap::update_load(FlatSourceDomain* domain){

  vector<uint64_t> n_hits(mpi::n_procs, 0);

  // add source region hits
  #pragma omp parallel
  {
    // number of hits per thread
    vector<uint64_t> local_hits(mpi::n_procs, 0);

    // if (mpi::master)  printf("Calculating rank loads... \n");

    #pragma omp for
      for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {

        Position centroid = domain->source_regions_.centroid(sr);
        // int owner = C_NONE;
        // double min_distance = INFTY;
        
        // // Find closest rank center
        // for (int rank = 0; rank < mpi::n_procs; rank++) {
        //     double dist = (centroid - rank_centers_[rank]).norm();
        //     // Power Voronoi diagram uses squared distances
        //     dist = dist*dist - rank_weights_[rank];
        //     if (dist < min_distance) {
        //         min_distance = dist;
        //         owner = rank;
        //         // if (min_distance < 0.0){
        //         //   break; // no need to continue searching
        //         // }
        //     }
        // }

        int owner = find_closest_rank(centroid);

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

  // if (mpi::master)  printf("Total hits before allreduce: %d \n", std::accumulate(n_hits.begin(), n_hits.end(), 0));

  // Accumulate rank loads across all ranks
  MPI_Allreduce(MPI_IN_PLACE, n_hits.data(), mpi::n_procs, MPI_UINT64_T, MPI_SUM, mpi::intracomm);

  for (int rank = 0; rank < mpi::n_procs; rank++){
    rank_load_[rank] = (n_hits[rank]/static_cast<double>(n_hits_sum_));
  }

  // if (mpi::master)  printf("Total hits after allreduce: %d \n", std::accumulate(n_hits.begin(), n_hits.end(), 0));
}

void DecompositionMap::redistribute_source_regions(FlatSourceDomain* domain) {

  // printf("1 Rank %d: Redistributing source regions... \n", mpi::rank);

  std::unordered_map<int, vector<int>> sr_send_list;
  vector<int> num_sr_receiving(mpi::n_procs, 0);

  // local source region container that contains updated list
  SourceRegionContainer source_regions_new(negroups_, is_linear_);

  // Each rank identifies source regions that need to be transferred to new owner and updates subdomain map accordingly
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {

    Position centroid = domain->source_regions_.centroid(sr);
    // int owner = C_NONE;
    // double min_distance = INFTY;
    
    // // Find closest rank center //TODO: Should this be put in update_load() to avoid repeating the same loop? But then I might right  
    // for (int rank = 0; rank < mpi::n_procs; rank++) {
    //     double dist = (centroid - rank_centers_[rank]).norm();
    //     // Power Voronoi diagram uses squared distances
    //     dist = dist*dist - rank_weights_[rank];
    //     if (dist < min_distance) {
    //       min_distance = dist;
    //       owner = rank;
    //     }
    // }

    int owner = find_closest_rank(centroid);

    // If owner changed, write source region to list of outbound source regions
    if (owner != mpi::rank){
      sr_send_list[owner].push_back(sr);
    } else {
      // If owner did not change, add source region to new local container
      source_regions_new.push_back(domain->source_regions_.get_source_region_handle(sr));
    }
  }

  // printf("2 Rank %d: Identified %lu ranks to send source regions to... \n", mpi::rank, sr_send_list.size());

  // Each rank communicates other ranks about ownership changes to update subdomain map
  for (int rank = 0; rank < mpi::n_procs; rank++) {

    // Send size
    int bcast_size = 0;
    if (rank == mpi::rank) {
      for (const auto& pair : sr_send_list) {
        bcast_size += pair.second.size();
      }
      // printf("2.1 Rank %d: Broadcasting %d source regions of %ld total source regions to update subdomain map... \n", mpi::rank, bcast_size, domain->n_source_regions());
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

        //TODO: Do sending here?
      }
    }
  }

  // printf("3 Rank %d: Updated subdomain map... \n", mpi::rank);

  // clear source_region_map_
  domain->source_region_map_.clear();

  // printf("3.1 Rank %d: Updated subdomain map... \n", mpi::rank);

  // Send source regions
  for (auto& pair : sr_send_list) {

    // printf("Rank %d: Test-2 \n", mpi::rank);
    int receiver = pair.first;              // destination rank
    vector<int>& sr_indices = pair.second;  // vector of source region indices
    // printf("Rank %d: Test-1 \n", mpi::rank);

    // Sort indices in descending order to make sure the correct source regions are erased //TODO: Is this safe?
    std::sort(sr_indices.begin(), sr_indices.end(), std::greater<int>());
    // printf("Rank %d: Test0 \n", mpi::rank);

    // Iterate through all source regions for this key
    for (int sr_idx : sr_indices) {
      SourceRegionKey sr_key = domain->source_regions_.key(sr_idx);
      // printf("Rank %d: Test1 \n", mpi::rank);
      SourceRegionHandle srh = domain->source_regions_.get_source_region_handle(sr_idx);
      // printf("Rank %d: Test2 \n", mpi::rank);
      SourceRegion sr(srh);
      // printf("Rank %d: Test3 \n", mpi::rank);
      send_sr_data(receiver, sr);
      // printf("Rank %d: Test4 \n", mpi::rank);
      // source_regions_new.erase(sr_idx); // remove source region from local container //TODO: Erase seems ugly?
      // printf("Rank %d: Test5 \n", mpi::rank);
    }
  }

  // printf("4 Rank %d: Sent source regions... \n", mpi::rank);

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

  // printf("5 Rank %d: Received source regions... \n", mpi::rank);

  int num_sr_old = domain->n_source_regions();

  // Update source regions in domain
  domain->source_regions_ = source_regions_new;
  
  int num_sr_new = domain->n_source_regions();

  // Update source region map
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
    SourceRegionKey key = domain->source_regions_.key(sr);
    domain->source_region_map_[key] = sr;
  }

  // printf("6 Rank %d: Updated source region map... \n", mpi::rank);
  // printf("Rank %d: Source regions before: %d, after: %d \n", mpi::rank, num_sr_old, num_sr_new);
}

}// namespace openmc