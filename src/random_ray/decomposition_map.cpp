#include "openmc/random_ray/decomposition_map.h"

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/random_lcg.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"
#include "openmc/vector.h"
#include <algorithm>
#include <cstring>
#include <numeric>

#ifdef OPENMC_MPI

namespace openmc {

namespace mpi {
DecompositionMap decomp_map;
}

// Constructor
DecompositionMap::DecompositionMap() {}

void DecompositionMap::initialize()
{
  negroups_ = data::mg.num_energy_groups_;
  estimated_rank_load_fractions_.resize(mpi::n_procs, 0.0);
  measured_rank_load_fractions_.resize(mpi::n_procs, 0.0);
  estimated_rank_load_totals_.resize(mpi::n_procs, 0.0);
  target_load_ = 1.0 / mpi::n_procs;
  rank_weights_.resize(mpi::n_procs, 1.0);

  auto* independent_source =
    dynamic_cast<IndependentSource*>(RandomRay::ray_source_.get());
  spatial_box_ = independent_source
                   ? dynamic_cast<SpatialBox*>(independent_source->space())
                   : nullptr;
  if (!spatial_box_) {
    fatal_error("Domain decomposition of the random ray solver requires the "
                "ray source to be an independent source with a box spatial "
                "distribution.");
  }

  double x_length =
    spatial_box_->upper_right().x - spatial_box_->lower_left().x;
  double y_length =
    spatial_box_->upper_right().y - spatial_box_->lower_left().y;
  double z_length =
    spatial_box_->upper_right().z - spatial_box_->lower_left().z;

  max_domain_length_ =
    sqrt(x_length * x_length + y_length * y_length + z_length * z_length);

  is_linear_ = RandomRay::source_shape_ != RandomRaySourceShape::FLAT;
  is_adaptive_ =
    is_adaptive_family(FlatSourceDomain::resolved_volume_estimator_);

  // Size in bytes of the record that carries one source region between ranks.
  // It depends only on the settings above and on the run mode, so it is
  // measured once per solve on a freshly constructed region.
  SourceRegion sr(negroups_, is_linear_);
  sr_record_bytes_ = process_sr_fields(sr, nullptr, SRProcessMode::SIZE);

  // Count the number of source regions in the model
  n_base_sr_ = 0;
  for (const auto& c : model::cells) {
    if (c->type_ == Fill::MATERIAL) {
      n_base_sr_ += c->n_instances();
    }
  }

  num_base_source_region_RT_.resize(n_base_sr_, 0);
  num_mesh_bin_RT_.resize(n_base_sr_, 0);
  ray_tracing_cost_.resize(n_base_sr_);
  volume_base_sr_.resize(n_base_sr_, 0.0);
}

void DecompositionMap::generate_rank_centers()
{

  // Calculate grid points that are used for Voronoi cells
  int grid_points_total = grid_points_per_rank_ * mpi::n_procs;
  write_message(6, "Calculating {} grid points for Voronoi tessellation...",
    grid_points_total);
  calculate_grid_points(grid_points_total);

  // Initialize points with random positions
  initialize_voronoi_centers();

  double err = 1.0;
  double precision = 1e-3; // corresponding to 0.001 cm position change
  int it = 0;
  int max_iterations = 100;

  // Lloyd's algorithm to move Voronoi centers to centroids of their cells
  // https://en.wikipedia.org/wiki/Lloyd%27s_algorithm
  while (err > precision && it < max_iterations) {
    // Reset error to determine maximum movement below
    err = 0.0;

    vector<Position> position_sum_per_rank(mpi::n_procs, Position(0, 0, 0));
    vector<int> num_points_per_rank(mpi::n_procs, 0);

    // Compute Voronoi cells by summing up all position values of mesh grid
    // points that are closest to a Voronoi center
    calculate_voronoi(position_sum_per_rank, num_points_per_rank);

    for (int rank = 0; rank < mpi::n_procs; rank++) {

      // Calculate centroid of the cell
      Position centroid = calculate_centroids(
        position_sum_per_rank[rank], num_points_per_rank[rank], rank);

      // Calculate movement for convergence check
      double movement = (centroid - rank_centers_[rank]).norm();

      // Move rank center to centroid
      rank_centers_[rank] = centroid;

      // Record maximum movement
      if (movement > err) {
        err = movement;
      }
    }

    it++;
  }

  // The centers must be bitwise identical on every rank: find_closest_rank()
  // is evaluated independently on different ranks for the same physical
  // point and they have to agree, or ownership diverges and a ray can
  // ping-pong forever. The Lloyd accumulation above uses omp atomic on
  // doubles, whose summation order is not reproducible, so broadcast rank
  // 0's result rather than relying on every rank landing on the same bits.
  MPI_Bcast(rank_centers_.data(), rank_centers_.size() * sizeof(Position),
    MPI_BYTE, 0, mpi::intracomm);

  if (mpi::master) {
    if (it == max_iterations) {
      warning("Lloyd's algorithm did not converge within the maximum number of "
              "iterations.");
    } else {
      write_message(6, "Lloyd's algorithm converged in {} iterations.", it);
    }
  }
}

// Calculate grid points needed for calculating Voronoi cells
void DecompositionMap::calculate_grid_points(int grid_points_total)
{

  // Calculate length along each dimension
  vector<double> domain_length(3);
  domain_length[0] =
    spatial_box_->upper_right().x - spatial_box_->lower_left().x;
  domain_length[1] =
    spatial_box_->upper_right().y - spatial_box_->lower_left().y;
  domain_length[2] =
    spatial_box_->upper_right().z - spatial_box_->lower_left().z;

  double volume = domain_length[0] * domain_length[1] * domain_length[2];

  // For each dimension, determine grid points along that direction based on
  // aspect ratio: domain_length / volume^(1/3) = grid_points_dimension /
  // grid_points_total^(1/3). Check if any dimension is so distorted that it
  // would only receive minimum of 1 grid point and flag that direction to
  // correct total number of grid points.
  int excluded_dimension = -1;
  vector<int> grid_points_per_dimension(3);
  for (int i = 0; i < 3; i++) {
    double grid_points_estimate =
      ((domain_length[i] / cbrt(volume)) * cbrt(grid_points_total));

    if (grid_points_estimate > 1) {
      grid_points_per_dimension[i] = round(grid_points_estimate);
    } else {
      excluded_dimension = i;
      grid_points_per_dimension[i] = 1;
    }
  }

  // If problem is 2D, exclude z direction
  if (RandomRay::geom_dim_ == RandomRayGeomDim::TWO_DIM) {
    excluded_dimension = 2;
    grid_points_per_dimension[2] = 1;
  }

  // If one dimension is excluded, recalculate grid points in other two
  // dimensions based on area.
  if (excluded_dimension != -1) {
    double area = 1.0;

    for (int i = 0; i < 3; i++) {
      if (i == excluded_dimension)
        continue;
      area *= domain_length[i];
    }

    for (int i = 0; i < 3; i++) {
      if (i == excluded_dimension)
        continue;
      grid_points_per_dimension[i] =
        round((domain_length[i] / sqrt(area)) * sqrt(grid_points_total));
    }
  }

  // Adjust grid points in each dimension to match actual total number of grid
  // points to the number of grid points requested
  double new_total = grid_points_per_dimension[0] *
                     grid_points_per_dimension[1] *
                     grid_points_per_dimension[2];
  double adjustment;
  if (excluded_dimension != -1) {
    // When one dimension is excluded, use square root for the other two
    // dimensions
    adjustment = sqrt(grid_points_total / new_total);
  } else {
    // When all dimensions are used, use cubic root
    adjustment = cbrt(grid_points_total / new_total);
  }

  // Multiply adjustment factor
  for (int i = 0; i < 3; i++) {
    if (i == excluded_dimension)
      continue;
    grid_points_per_dimension[i] =
      round(grid_points_per_dimension[i] * adjustment);
  }

  // Calculate spacing between grid points in each dimension
  vector<double> delta_value(3, 0.0);

  for (int i = 0; i < 3; i++) {
    if (grid_points_per_dimension[i] > 1) {
      delta_value[i] =
        (domain_length[i] - 2 * TINY_BIT) / (grid_points_per_dimension[i] - 1);
    }
  }

  // Initialize point at center of domain (in case of only 1 grid point in a
  // dimension)
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

  // Check if mesh grid points are inside spatial domain, which can be different
  // from a box. If not, erase them.
  for (int i = grid_points_.size() - 1; i >= 0; i--) {
    Position xi = grid_points_[i];

    bool is_inside_domain =
      RandomRay::ray_source_->satisfies_spatial_constraints(xi);

    if (!is_inside_domain) {
      grid_points_.erase(grid_points_.begin() + i);
    }
  }

  if (mpi::master && grid_points_.size() < grid_points_total) {
    warning(fmt::format("Spatial constraints reduced grid points for Voronoi "
                        "tesselation from {} to {}.",
      grid_points_total, grid_points_.size()));
  }
}

// Places random points in the spatial domain.
// Each point corresponds to the initial center of a rank.
void DecompositionMap::initialize_voronoi_centers()
{
  rank_centers_.resize(mpi::n_procs);

  uint64_t seed = openmc_get_seed();
  int rank_cnt = 0;

  // Sample random positions to start with
  while (rank_cnt < mpi::n_procs) {

    double x = prn(&seed);
    double y = prn(&seed);
    double z = 0.0;
    if (RandomRay::geom_dim_ == RandomRayGeomDim::THREE_DIM) {
      z = prn(&seed);
    } else {
      z = 0.5; // Mid-plane in z direction
    }

    Position xi {x, y, z};

    // Make a small shift in position to avoid geometry floating point issues at
    // boundaries
    Position shift {FP_COINCIDENT, FP_COINCIDENT, FP_COINCIDENT};
    xi = (spatial_box_->lower_left() + shift) +
         xi * ((spatial_box_->upper_right() - shift) -
                (spatial_box_->lower_left() + shift));

    bool is_inside_domain =
      RandomRay::ray_source_->satisfies_spatial_constraints(xi);

    if (is_inside_domain) {
      rank_centers_[rank_cnt] = xi;
      rank_cnt++;
    }
  }
}

// Determine the distance of each mesh grid point to all rank centers.
// Sum up the positions of all mesh grid points that are closest to a given rank
// center for computation of centroid later. Record number of grid points per
// rank center.
void DecompositionMap::calculate_voronoi(
  vector<Position>& position_sum_per_rank, vector<int>& num_points_per_rank)
{

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
      dist = dist * dist - rank_weights_[rank];

      if (dist < min_distance) {
        min_distance = dist;
        closest_rank = rank;
      }
    }

    if (closest_rank == C_NONE) {
      fatal_error("Could not find closest rank for Voronoi cell point " +
                  std::to_string(p) + ".");
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

Position DecompositionMap::calculate_centroids(
  const Position position_sum, const int num_points, int rank)
{

  // check if any points have been recorded in rank
  if (num_points == 0) {
    fatal_error("Rank " + std::to_string(rank) +
                " has no Voronoi cell points. This indicates that the number "
                "of grid points for the Voronoi tesselation is too coarse. "
                "Requires source code modificaiton to fix.");
  }

  Position centroid = position_sum;

  // Divide by number of points
  double n = static_cast<double>(num_points);
  centroid.x /= n;
  centroid.y /= n;
  centroid.z /= n;

  return centroid;
}

bool DecompositionMap::any_discovered_source_regions(
  ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions)
{

  simulation::time_decomposition_handling.start();

  int flag = 0;
  if (discovered_source_regions.begin() != discovered_source_regions.end()) {
    flag = 1;
  }

  MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, mpi::intracomm);

  simulation::time_decomposition_handling.stop();

  return flag > 0;
}

void DecompositionMap::exchange_sr_info(
  ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions)
{

  // Communicate maps
  for (int rank = 0; rank < mpi::n_procs; rank++) {

    // Send size
    uint64_t bcast_size = 0;
    if (rank == mpi::rank) {
      for (const auto& pair : discovered_source_regions) {
        // Only broadcast source regions that have non-zero volume, i.e. regions
        // discovered during active phase of ray
        if (pair.second.scalars_.volume_ > 0.0) {
          bcast_size++;
        }
      }
    }

    MPI_Bcast(&bcast_size, 1, MPI_UINT64_T, rank, mpi::intracomm);

    if (bcast_size > 0) {

      vector<int64_t> local_base_ids(bcast_size);
      vector<int64_t> local_mesh_bins(bcast_size);

      if (rank == mpi::rank) {
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
      MPI_Bcast(
        local_base_ids.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);
      MPI_Bcast(
        local_mesh_bins.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);

      // Update subdomain map
      for (int j = 0; j < bcast_size; j++) {

        SourceRegionKey sr_key(local_base_ids[j], local_mesh_bins[j]);

        // check if already in map or not, i.e. has someone else discovered that
        // region already?
        if (subdomain_map_.find(sr_key) == subdomain_map_.end()) {
          subdomain_map_[sr_key] = rank;
        } else {
          int resident_rank = subdomain_map_[sr_key]; // current owner
          int challenger_rank = rank; // current broadcasting rank

          double ratio_resident = calculate_load_ratio(resident_rank);
          double ratio_challenger = calculate_load_ratio(challenger_rank);
          double resident_rank_load =
            ratio_resident * estimated_rank_load_totals_[resident_rank];
          double challenger_rank_load =
            ratio_challenger * estimated_rank_load_totals_[challenger_rank];

          // If load of challenger rank is lower, assign source region to that
          // rank, otherwise resident keeps it
          int sender;
          int receiver;
          if (challenger_rank_load < resident_rank_load) {
            subdomain_map_[sr_key] = challenger_rank;
            sender = resident_rank;
            receiver = challenger_rank;
          } else {
            // resident keeps it
            sender = challenger_rank;
            receiver = resident_rank;
          }

          // Broadcast load of exchanged source region such that each rank can
          // update load balance
          double bcast_load = 0;
          if (mpi::rank == sender) {
            SourceRegion& contested_sr = discovered_source_regions[sr_key];
            double volume_sr = contested_sr.scalars_.volume_;
            bcast_load =
              C1_ *
                (static_cast<double>(contested_sr.scalars_.n_hits_) /
                  simulation::current_batch) *
                negroups_ +
              volume_sr * ray_tracing_cost_[sr_key.base_source_region_id];
          }
          MPI_Bcast(&bcast_load, 1, MPI_DOUBLE, sender, mpi::intracomm);

          double load_change_fraction =
            bcast_load / estimated_load_sum_; // load fraction update
          estimated_rank_load_fractions_[sender] -= load_change_fraction;
          estimated_rank_load_fractions_[receiver] += load_change_fraction;
          double load_change_total = bcast_load; // load total update
          estimated_rank_load_totals_[sender] -= load_change_total;
          estimated_rank_load_totals_[receiver] += load_change_total;

          // Communicate source region data and merge on receiver side
          if (mpi::rank == sender) {
            SourceRegion& sr = discovered_source_regions[sr_key];
            send_sr_data(receiver, sr);

            // clear old source region data from discovered regions map
            discovered_source_regions.erase(sr_key);
          }
          if (mpi::rank == receiver) {
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

//==============================================================================
// Source region transfer between ranks
//
// A source region is moved as one fixed-size record of raw bytes (the scalar
// fields, then the group-wise arrays in a fixed order). The record layout and
// size depend only on run-wide settings (energy groups, linear source, source
// gradient limiter, adaptive volume estimator, fixed source mode), so sender
// and receiver agree on it without any metadata. When adding a per-region
// field, add it to ScalarSourceRegionFields (carried automatically) or to
// process_sr_fields() (and allocate it in the SourceRegion constructors).
//==============================================================================

// Start receiving `bytes` bytes from rank `sender` into `buf`, as messages of
// at most MAX_MPI_MESSAGE_BYTES. The sender splits its buffer the same way
// (start_chunked_send), and MPI delivers the pieces in the order they were
// sent.
void DecompositionMap::start_chunked_receive(
  char* buf, size_t bytes, int sender, vector<MPI_Request>& requests) const
{
  for (size_t off = 0; off < bytes; off += MAX_MPI_MESSAGE_BYTES) {
    int n = static_cast<int>(std::min(MAX_MPI_MESSAGE_BYTES, bytes - off));
    requests.emplace_back();
    MPI_Irecv(
      buf + off, n, MPI_BYTE, sender, 1, mpi::intracomm, &requests.back());
  }
}

// Start sending `bytes` bytes of `buf` to rank `receiver`, as messages of at
// most MAX_MPI_MESSAGE_BYTES
void DecompositionMap::start_chunked_send(
  char* buf, size_t bytes, int receiver, vector<MPI_Request>& requests) const
{
  for (size_t off = 0; off < bytes; off += MAX_MPI_MESSAGE_BYTES) {
    int n = static_cast<int>(std::min(MAX_MPI_MESSAGE_BYTES, bytes - off));
    requests.emplace_back();
    MPI_Isend(
      buf + off, n, MPI_BYTE, receiver, 1, mpi::intracomm, &requests.back());
  }
}

// Copy the `n` bytes of one field into the record at `offset` (PACK) or out of
// it (UNPACK), or only count them (SIZE), and advance the offset
void DecompositionMap::process_bytes(
  void* field, size_t n, char* buf, size_t& offset, SRProcessMode mode) const
{
  if (mode == SRProcessMode::PACK) {
    std::memcpy(buf + offset, field, n);
  } else if (mode == SRProcessMode::UNPACK) {
    std::memcpy(field, buf + offset, n);
  }
  offset += n;
}

// The same for one group-wise array, after checking that it has the expected
// number of entries (a field that is transferred but not allocated by the
// SourceRegion constructors would otherwise be read or written out of bounds)
template<typename T>
void DecompositionMap::process_array(
  vector<T>& v, size_t n, char* buf, size_t& offset, SRProcessMode mode) const
{
  if (v.size() != n) {
    fatal_error(fmt::format(
      "A transferred source region array has {} entries instead of {}. Check "
      "that the SourceRegion constructors allocate every field listed in "
      "DecompositionMap::process_sr_fields().",
      v.size(), n));
  }
  process_bytes(v.data(), n * sizeof(T), buf, offset, mode);
}

size_t DecompositionMap::process_sr_fields(
  SourceRegion& sr, char* buf, SRProcessMode mode) const
{
  // The one list of what moves with a source region between ranks, in record
  // order. Scalars ride along in ScalarSourceRegionFields.
  size_t offset = 0;
  size_t ng = negroups_;
  process_bytes(
    &sr.scalars_, sizeof(ScalarSourceRegionFields), buf, offset, mode);
  process_array(sr.scalar_flux_old_, ng, buf, offset, mode);
  process_array(sr.scalar_flux_new_, ng, buf, offset, mode);
  process_array(sr.source_, ng, buf, offset, mode);
  if (settings::run_mode == RunMode::FIXED_SOURCE)
    process_array(sr.external_source_, ng, buf, offset, mode);
  process_array(sr.scalar_flux_final_, ng, buf, offset, mode);
  if (is_linear_) {
    process_array(sr.source_gradients_, ng, buf, offset, mode);
    process_array(sr.flux_moments_old_, ng, buf, offset, mode);
    process_array(sr.flux_moments_new_, ng, buf, offset, mode);
    process_array(sr.flux_moments_t_, ng, buf, offset, mode);
  }
  if (is_linear_ && FlatSourceDomain::source_gradient_limiter_)
    process_bytes(&sr.extent_, sizeof(BoundingBox), buf, offset, mode);
  if (is_adaptive_)
    process_array(sr.scalar_flux_t_, ng, buf, offset, mode);

  // The size of the record in bytes
  return offset;
}

void DecompositionMap::send_sr_data(int receiver, SourceRegion& sr_send)
{
  vector<char> buf(sr_record_bytes_);
  process_sr_fields(sr_send, buf.data(), SRProcessMode::PACK);
  MPI_Send(buf.data(), static_cast<int>(buf.size()), MPI_BYTE, receiver, 1,
    mpi::intracomm);
}

void DecompositionMap::receive_sr_data(int sender, SourceRegion& sr_recv)
{
  vector<char> buf(sr_record_bytes_);
  MPI_Recv(buf.data(), static_cast<int>(buf.size()), MPI_BYTE, sender, 1,
    mpi::intracomm, MPI_STATUS_IGNORE);
  process_sr_fields(sr_recv, buf.data(), SRProcessMode::UNPACK);
}

int DecompositionMap::find_owner(SourceRegionKey sr_key, Position r,
  ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>&
    discovered_source_regions)
{

  // Check if source region key is in subdomain map
  auto it = subdomain_map_.find(sr_key);
  if (it != subdomain_map_.end()) {
    return it->second;
  }

  // Check if already recorded in newly discovered source regions
  discovered_source_regions.lock(sr_key);
  bool sr_key_discovered = discovered_source_regions.contains(sr_key);
  discovered_source_regions.unlock(sr_key);
  if (sr_key_discovered) {
    return mpi::rank;
  }

  // If not found in either map, check which rank owns source
  // region beased on location
  int closest_rank = find_closest_rank(r, true);
  return closest_rank;
}

int DecompositionMap::find_closest_rank(Position r, bool test_all_ranks)
{

  int closest_rank = C_NONE;
  double min_distance = INFTY;

  // Distance function corresponding to weighted power Voronoi diagram. This
  // runs inside the transport loop (via find_owner) and inside the load
  // balancing iteration, so it deliberately allocates nothing.
  auto test_rank = [&](int rank) {
    double dist = (r - rank_centers_[rank]).norm();
    dist = dist * dist - rank_weights_[rank];
    if (dist < min_distance) {
      min_distance = dist;
      closest_rank = rank;
    }
  };

  if (test_all_ranks) {
    for (int rank = 0; rank < mpi::n_procs; rank++) {
      test_rank(rank);
    }
  } else {
    // Only the recorded neighbors, plus this rank itself
    for (int rank : mpi::decomp_map.my_neighbors_) {
      test_rank(rank);
    }
    test_rank(mpi::rank);
  }

  if (closest_rank == C_NONE) {
    fatal_error(
      "Could not find closest rank for new source region at position (" +
      std::to_string(r.x) + ", " + std::to_string(r.y) + ", " +
      std::to_string(r.z) + ").");
  }

  return closest_rank;
}

void DecompositionMap::calculate_rank_load(
  FlatSourceDomain* domain, double batch_transport_time)
{

  // Reset local volumes of base source regions, which might change when source
  // regions change rank ownership
  std::fill(volume_base_sr_.begin(), volume_base_sr_.end(), 0.0);

  // Add volumes of newly discovered source regions
  for (const auto& [sr_key, sr] : domain->discovered_source_regions_) {
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
      ray_tracing_cost_[bsr] =
        (C2_ * static_cast<double>(num_base_source_region_RT_[bsr]) +
          C3_ * static_cast<double>(num_mesh_bin_RT_[bsr])) /
        volume_base_sr_[bsr];
    } else {
      ray_tracing_cost_[bsr] = 0.0;
    }
  }

  // Accumulate load of known source regions
  double local_estimated_load = 0.0;
#pragma omp parallel for reduction(+ : local_estimated_load)
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
    SourceRegionKey sr_key = domain->source_regions_.key(sr);
    uint64_t base_sr = sr_key.base_source_region_id;

    // Calculate volume of unique source region to weight volume-dependent ray
    // tracing cost
    double volume_sr =
      domain->source_regions_.volume_t(sr) + domain->source_regions_.volume(sr);

    // Calculate load of source region
    double load_sr =
      C1_ *
        (static_cast<double>(domain->source_regions_.n_hits(sr)) /
          simulation::current_batch) *
        negroups_ +
      volume_sr * ray_tracing_cost_[base_sr];

    // Accumulate to local estimated load
    local_estimated_load += load_sr;
  }

  // Accumulate load of newly discovered source regions
  for (const auto& [sr_key, sr] : domain->discovered_source_regions_) {
    uint64_t base_sr = sr_key.base_source_region_id;
    double volume_sr = sr.scalars_.volume_;
    double load_sr =
      C1_ *
        (static_cast<double>(sr.scalars_.n_hits_) / simulation::current_batch) *
        negroups_ +
      volume_sr * ray_tracing_cost_[base_sr];
    local_estimated_load += load_sr;
  }

  // Communicate estimated load across ranks
  MPI_Allgather(&local_estimated_load, 1, MPI_DOUBLE,
    estimated_rank_load_totals_.data(), 1, MPI_DOUBLE, mpi::intracomm);
  estimated_load_sum_ = std::accumulate(estimated_rank_load_totals_.begin(),
    estimated_rank_load_totals_.end(), 0.0);

  // Communicate measured load across ranks
  MPI_Allgather(&batch_transport_time, 1, MPI_DOUBLE,
    measured_rank_load_fractions_.data(), 1, MPI_DOUBLE, mpi::intracomm);
  double measured_load_sum =
    std::accumulate(measured_rank_load_fractions_.begin(),
      measured_rank_load_fractions_.end(), 0.0);

  // Calculate fractions
  for (int rank = 0; rank < mpi::n_procs; rank++) {
    estimated_rank_load_fractions_[rank] =
      estimated_rank_load_totals_[rank] / estimated_load_sum_;
    measured_rank_load_fractions_[rank] =
      measured_rank_load_fractions_[rank] / measured_load_sum;
  }

  // Reset ray trace counters
  fill(num_base_source_region_RT_.begin(), num_base_source_region_RT_.end(), 0);
  fill(num_mesh_bin_RT_.begin(), num_mesh_bin_RT_.end(), 0);
}

void DecompositionMap::balance_load(FlatSourceDomain* domain)
{

  // Optimization parameters
  int max_iterations = 200;
  int it_outer = 0;
  double adaptation_factor = 1;
  double min_adaptation_factor = 0.01;
  double max_adaptation_factor = 2;
  double avg_rank_distance =
    max_domain_length_ /
    cbrt(mpi::n_procs); // rough estimate of average distance between ranks
  double weight_scale =
    avg_rank_distance * avg_rank_distance * optimization_history_factor_;
  double beta = 0.6; // momentum damping
  bool check_all_ranks = true;

  vector<double> weight_change(mpi::n_procs, 0.0);
  vector<double> combined_rank_load(mpi::n_procs, 0.0);
  vector<double> load_ratio(mpi::n_procs, 0.0);

  // Combine estimated load with measured load ratios
  for (int rank = 0; rank < mpi::n_procs; rank++) {
    load_ratio[rank] = calculate_load_ratio(rank);
    combined_rank_load[rank] =
      load_ratio[rank] * estimated_rank_load_totals_[rank];
  }

  double combined_load_sum =
    std::accumulate(combined_rank_load.begin(), combined_rank_load.end(), 0.0);

  for (int rank = 0; rank < mpi::n_procs; rank++) {
    combined_rank_load[rank] = (combined_rank_load[rank] / combined_load_sum);
  }
  double max_load =
    *std::max_element(combined_rank_load.begin(), combined_rank_load.end());
  double max_imbalance = (max_load - target_load_) / target_load_;
  double prev_imbalance = max_imbalance;

  // History tracking
  vector<double> imbalance_history;
  vector<vector<double>> weight_history;
  imbalance_history.push_back(max_imbalance);
  weight_history.push_back(rank_weights_);

  // Change weights to equalize load based on combined load estimates
  while (max_imbalance > imbalance_tolerance_ && it_outer < max_iterations) {

    it_outer++;

    for (int rank = 0; rank < mpi::n_procs; rank++) {
      double corr = ((combined_rank_load[rank] - target_load_) / target_load_) *
                    weight_scale;
      weight_change[rank] =
        beta * weight_change[rank] +
        (1.0 - beta) * corr; // keep some inertia from previous changes to
                             // prevent oscillations
      double damping = std::clamp(max_imbalance / prev_imbalance, 0.1,
        1.0); // dampening factor based on whether we are getting closer to
              // convergence or not, prevents big jumps, if too big a change,
              // more dampening is applied
      rank_weights_[rank] -= adaptation_factor * damping * weight_change[rank];
    }

    if (simulation::current_batch > 1) {
      // Check distances to all ranks only in first batch, otherwise only
      // recorded neighbors
      check_all_ranks = false;
    }

    // Calculate new load after weight update
    update_load(domain, check_all_ranks, combined_rank_load, load_ratio);
    max_load =
      *std::max_element(combined_rank_load.begin(), combined_rank_load.end());
    max_imbalance = (max_load - target_load_) / target_load_;

    // Store imbalance history
    imbalance_history.push_back(max_imbalance);
    weight_history.push_back(rank_weights_);

    // Adaptive factor
    if (max_imbalance > prev_imbalance)
      adaptation_factor =
        std::max(adaptation_factor * 0.5, min_adaptation_factor);
    else
      adaptation_factor =
        std::min(adaptation_factor * 1.05, max_adaptation_factor);

    prev_imbalance = max_imbalance;
  }

  // Check convergence and adjust history optimization factor depending on
  // failure mode to enable better convergence in the following batch
  if (it_outer == max_iterations) {
    if (mpi::master) {
      warning("MPI load balancing has not converged after " +
              std::to_string(max_iterations) + " iterations.");
    }

    // Check if oscillating or simply slow convergence
    int direction_changes = 0;

    for (int i = 1; i < imbalance_history.size(); i++) {
      // Calculate change
      double change = imbalance_history[i] - imbalance_history[i - 1];

      // Count direction changes
      if (i > 1) {
        double prev_change =
          imbalance_history[i - 1] - imbalance_history[i - 2];
        if (change * prev_change < 0) {
          direction_changes++;
        }
      }
    }

    // Check for oscillations
    double oscillation_ratio =
      (double)direction_changes / (imbalance_history.size() - 2);

    if (oscillation_ratio > 0.4) {
      // decrease weight for next batch if oscillating
      optimization_history_factor_ =
        std::max(optimization_history_factor_ * 0.5, min_adaptation_factor);
    } else {
      // increase weight for faster convergence if too slow
      optimization_history_factor_ =
        std::min(optimization_history_factor_ * 1.2, max_adaptation_factor);
    }

    // Check which iteration had the best imbalance and revert to those weights
    auto min_it =
      std::min_element(imbalance_history.begin(), imbalance_history.end());
    int best_index = std::distance(imbalance_history.begin(), min_it);
    double best_imbalance = *min_it;
    rank_weights_ = weight_history[best_index];

    write_message(6,
      "Best imbalance during optimization was {:.2f}% at iteration {}.",
      best_imbalance * 100.0, best_index);

    if (best_index == 0) {
      // if no improvement at all, just keep current decomposition and return
      // without redistributing
      return;
    }
  } else {
    optimization_history_factor_ = 1.0; // reset history factor if converged
    write_message(6,
      "MPI load balancing converged after {} iterations. Max. imbalance: "
      "{:.2f}%",
      it_outer, max_imbalance * 100.0);
  }

  // Redistribute source regions according to new weights determined in
  // optimization
  redistribute_source_regions(domain);
}

void DecompositionMap::update_load(FlatSourceDomain* domain,
  bool check_all_ranks, vector<double>& combined_rank_load,
  vector<double>& load_ratio)
{

  vector<double> load(mpi::n_procs, 0);

// Add up source region hits to respective rank depending of position of
// centroid of each source region
#pragma omp parallel
  {
    // number of hits per thread
    vector<double> thread_load(mpi::n_procs, 0);

#pragma omp for
    for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
      Position centroid = domain->source_regions_.centroid(sr);
      int owner = find_closest_rank(centroid, check_all_ranks);
      double volume_sr = domain->source_regions_.volume_t(sr);
      thread_load[owner] +=
        load_ratio[owner] *
        (C1_ *
            (static_cast<double>(domain->source_regions_.n_hits(sr)) /
              simulation::current_batch) *
            negroups_ +
          volume_sr * ray_tracing_cost_[domain->source_regions_.key(sr)
                                          .base_source_region_id]);
    }

// Combine results from different threads
#pragma omp critical(combining_loads)
    {
      for (int i = 0; i < mpi::n_procs; i++) {
        load[i] += thread_load[i];
      }
    }
  }

  // Communicate new load estimates across ranks
  MPI_Allreduce(MPI_IN_PLACE, load.data(), mpi::n_procs, MPI_DOUBLE, MPI_SUM,
    mpi::intracomm);
  double load_sum = std::accumulate(load.begin(), load.end(), 0.0);

  // Update new combined rank load fractions
  for (int rank = 0; rank < mpi::n_procs; rank++) {
    combined_rank_load[rank] = load[rank] / load_sum;
  }
}

void DecompositionMap::redistribute_source_regions(FlatSourceDomain* domain)
{

  // Map of source regions to be sent to other ranks
  std::unordered_map<int, vector<int>> sr_send_list;
  // Number of source regions to be received from other ranks
  vector<int> num_sr_receiving(mpi::n_procs, 0);

  // Local source region container that contains updated list
  SourceRegionContainer source_regions_new =
    domain->source_regions_.empty_like();

  // Each rank identifies source regions that need to be transferred to new
  // owner and updates subdomain map accordingly
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
    Position centroid = domain->source_regions_.centroid(sr);
    int owner = find_closest_rank(centroid, true);

    // If owner changed, write source region to list of outbound source regions
    if (owner != mpi::rank) {
      sr_send_list[owner].push_back(sr);
    }
    // If owner did not change, add source region to new local container
    else {
      source_regions_new.push_back(
        domain->source_regions_.get_source_region_handle(sr));
    }
  }

  // Each rank informs other ranks about ownership changes to update subdomain
  // map
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
      if (rank == mpi::rank) {
        int i = 0;
        for (auto& pair : sr_send_list) {

          int receiver = pair.first; // new owner
          vector<int>& sr_indices =
            pair.second; // vector of source region indices

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
      MPI_Bcast(
        mesh_bins.data(), bcast_size, MPI_INT64_T, rank, mpi::intracomm);

      // Every rank updates subdomain map
      for (int j = 0; j < bcast_size; j++) {

        // Re-assemble source region key and extract new owner
        SourceRegionKey sr_key(base_ids[j], mesh_bins[j]);
        int rank_new = rank_ids[j];

        // Update subdomain_map with new owner
        subdomain_map_[sr_key] = rank_new;

        // If calling rank is new owner, increase count of source regions coming
        // from sending rank (current broadcaster)
        if (mpi::rank == rank_new) {
          num_sr_receiving[rank] += 1;
        }
      }
    }
  }

  // Clear source_region_map_
  domain->source_region_map_.clear();

  // Move the source region data to the new owners: one packed buffer per
  // pair of ranks (fixed-size records, see process_sr_fields()), sent in
  // messages of at most MAX_MPI_MESSAGE_BYTES, instead of one message per
  // field and region. Every rank already knows how many regions it receives
  // from each sender (num_sr_receiving, from the broadcasts above), so all
  // receives are posted first, then all sends, then one Waitall: no ordering
  // between ranks is needed. The messages are few and large, which keeps MPI
  // resource use small (posting tens of thousands of small messages made Open
  // MPI's memory grow without bound).
  if (num_sr_receiving[mpi::rank] > 0) {
    fatal_error("Rank sends source regions to itself. Rank should not "
                "receive source regions from itself.");
  }
  const size_t record = sr_record_bytes_;
  vector<MPI_Request> requests;

  // Post receives
  vector<vector<char>> recv_buffers(mpi::n_procs);
  for (int sender = 0; sender < mpi::n_procs; sender++) {
    if (num_sr_receiving[sender] == 0)
      continue;
    recv_buffers[sender].resize(num_sr_receiving[sender] * record);
    start_chunked_receive(recv_buffers[sender].data(),
      recv_buffers[sender].size(), sender, requests);
  }

  // Pack and post sends
  vector<vector<char>> send_buffers;
  send_buffers.reserve(sr_send_list.size());
  for (auto& pair : sr_send_list) {
    int receiver = pair.first;             // destination rank
    vector<int>& sr_indices = pair.second; // vector of source region indices
    send_buffers.emplace_back(sr_indices.size() * record);
    char* p = send_buffers.back().data();
    for (int sr_idx : sr_indices) {
      SourceRegion sr(domain->source_regions_.get_source_region_handle(sr_idx));
      process_sr_fields(sr, p, SRProcessMode::PACK);
      p += record;
    }
    start_chunked_send(send_buffers.back().data(), send_buffers.back().size(),
      receiver, requests);
  }

  if (!requests.empty()) {
    MPI_Waitall(
      static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
  }
  send_buffers.clear();

  // Record starting source region ID for tally reinitialization later
  int64_t start_sr_id = source_regions_new.n_source_regions();

  // Unpack, in order of sending rank (as the region order before)
  for (int sender = 0; sender < mpi::n_procs; sender++) {
    char* p = recv_buffers[sender].data();
    for (int i = 0; i < num_sr_receiving[sender]; ++i) {
      SourceRegion sr_recv(negroups_, is_linear_);
      process_sr_fields(sr_recv, p, SRProcessMode::UNPACK);
      p += record;
      source_regions_new.push_back(sr_recv);
    }
    vector<char>().swap(recv_buffers[sender]);
  }

  // Update source regions in domain to new container
  domain->source_regions_ = source_regions_new;

  // Update source region map
  for (int64_t sr = 0; sr < domain->n_source_regions(); sr++) {
    SourceRegionKey key = domain->source_regions_.key(sr);
    domain->source_region_map_[key] = sr;
  }

  // Reinitialise tallies
  domain->convert_source_regions_to_tallies(start_sr_id);
}

double DecompositionMap::calculate_load_ratio(int rank)
{
  if (estimated_rank_load_fractions_[rank] > 0.0) {
    return measured_rank_load_fractions_[rank] /
           estimated_rank_load_fractions_[rank];
  } else {
    return 1.0;
  }
}

} // namespace openmc
#endif // OPENMC_MPI
