#include "openmc/random_ray/decomposition_map.h"

#include "openmc/message_passing.h"
#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/random_lcg.h"

namespace openmc {

namespace mpi {
    DecompositionMap decomp_map;
}

// Constructor
DecompositionMap::DecompositionMap() {
  // generate_rank_centers();
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
    if (mpi::master){
      printf("ITERATION %d \n", it);
    }

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
      printf("Lloyd's algorithm converged in %d iterations.", it);
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

// Method to initialize subdomain list
void DecompositionMap::create(FlatSourceDomain* domain){
  vector <int> start(mpi::n_procs);
  vector <int> end(mpi::n_procs);

  int sr_per_rank = ((domain->n_source_regions() + mpi::n_procs - 1) / mpi::n_procs);

  int sr_cnt = 0;
  int rank = 0;
  for (int i = 0; i < domain->n_source_regions(); i++) {
    if (sr_cnt >= sr_per_rank){
      sr_cnt = 0;
      rank++;
    }
    subdomain_map_[i] = rank;
    sr_cnt++;
  }  


  // int start = mpi::rank * sr_per_rank; // + 1; //TODO: does cell ID numbering starts at 1
  // int end = (mpi::rank + 1) * sr_per_rank; // + 1;

  // printf("Rank %d: Handling source regions %d to %d\n", mpi::rank, start, end-1);

  // for (int i = start; i < end; i++) {
  //   // subdomain_map_[SourceRegionKey(i, 0)] = mpi::rank;
  //   subdomain_map_[i] = mpi::rank;
  // }
}

// Update subdomain list for each decomposition map.
void DecompositionMap::update(){

  // Volume data needs to be moved around as well

}

// bool DecompositionMap::is_SRK_in_domain(SourceRegionKey sr_key){
bool DecompositionMap::is_SRK_in_domain(int sr_key){
  if (subdomain_map_.find(sr_key) != subdomain_map_.end()){
    if (subdomain_map_[sr_key] == mpi::rank){
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

// int DecompositionMap::n_source_regions(){
//   int count = 0;
  
//   for (const auto& [key, val] : subdomain_map_) {
//     if (val == mpi::rank) {
//       count++;
//     }
//   }

//   return count;
// }

}// namespace openmc
