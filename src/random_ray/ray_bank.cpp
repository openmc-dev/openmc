#include "openmc/random_ray/ray_bank.h"

#include <cstring>

#include "openmc/random_ray/random_ray.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/message_passing.h"
#include "openmc/random_ray/decomposition_map.h"
#include "openmc/mgxs_interface.h"
#include "openmc/timer.h"
#include "openmc/geometry.h"
#include <numeric>

namespace openmc {

// Constructor
RayBank::RayBank() {
  negroups_ = data::mg.num_energy_groups_;
  num_messages_receiving_.resize(mpi::n_procs, 0);
  num_messages_sending_.resize(mpi::n_procs, 0);
}

// Buffer ray that has left my subdomain
void RayBank::buffer_ray_data_to_send(RandomRay& ray, FlatSourceDomain* domain){

    // Get rank to send ray to
    int rank = ray.owner_rank_;

    if (rank == mpi::rank){
      warning(fmt::format(
        "Ray {} at position ({:.5e}, {:.5e}, {:.5e})"
        "is being sent to the same rank {}. This may indicate an error in the decomposition map."
        , ray.exchange_data_.ray_id, ray.exchange_data_.position.x, ray.exchange_data_.position.y, ray.exchange_data_.position.z, rank));
    }

    // Get or create the send buffer for this rank
    auto& buffers = ray_send_buffer_[rank];
    
    // Get the maximum coordinate levels
    const int n_coord_max = model::n_coord_levels;
    
    // Reserve space if this is the first ray for this rank
    // This is just a heuristic to reduce reallocations
    if (buffers.count == 0) {
      buffers.ray_data.reserve(reserved_buffer_size_);
      buffers.angular_flux.reserve(reserved_buffer_size_ * negroups_);
      buffers.coord.reserve(reserved_buffer_size_ * n_coord_max);
      buffers.cell_last.reserve(reserved_buffer_size_ * n_coord_max);


      buffers.ray_data.reserve(32);
      buffers.angular_flux.reserve(32 * negroups_);
      buffers.coord.reserve(32 * n_coord_max);
      buffers.cell_last.reserve(32 * n_coord_max);
    }
    
    // Pack RayExchangeData directly into send buffer
    RayBufferContainer& rbc = ray.exchange_data_;
    RayExchangeData rd;
    
    rd.position = rbc.position;
    rd.direction = rbc.direction;
    rd.distance_travelled = rbc.distance_travelled;
    rd.surface = rbc.surface;
    rd.is_active = rbc.is_active;
    rd.ray_id = rbc.ray_id;
    rd.n_event = rbc.n_event;
    rd.n_coord = rbc.n_coord;
    rd.cell_instance = rbc.cell_instance;
    rd.n_coord_last = rbc.n_coord_last;
    rd.material = rbc.material;
    rd.material_last = rbc.material_last;
    rd.sqrtkT = rbc.sqrtkT;
    rd.sqrtkT_last = rbc.sqrtkT_last;
    
#ifdef OPENMC_DAGMC_ENABLED
    rd.last_dir = rbc.last_dir;
    rd.n_handles = rbc.n_handles;
    std::memcpy(rd.handles, rbc.handles, MAX_N_HANDLES * sizeof(moab::EntityHandle));
#endif
    
    buffers.ray_data.push_back(rd);
    
    // Pack angular flux array directly
    buffers.angular_flux.insert(buffers.angular_flux.end(), 
                                rbc.angular_flux.begin(), 
                                rbc.angular_flux.end());
    
    // Pack coord and cell_last arrays directly
    buffers.coord.insert(buffers.coord.end(), 
                        rbc.coord.begin(), 
                        rbc.coord.end());
    buffers.cell_last.insert(buffers.cell_last.end(), 
                            rbc.cell_last.begin(), 
                            rbc.cell_last.end());
    
    buffers.count++;
}

// Update ray bank
void RayBank::update(FlatSourceDomain* domain){

    // Empty ray list because rays have either died or are in buffer to be sent to other ranks
    my_ray_list_.resize(0);

    // Communicate number of rays to be sent/received between ranks
    communicate_message_metadata();

    // Send and receive ray data between MPI ranks
    communicate_rays();

    // Add received rays to ray list of that rank
    update_my_ray_list(domain);
}

int RayBank::ray_bank_size(){
  int ray_bank_size = my_ray_list_.size();
  return ray_bank_size;
}

// Tells each rank how many rays to receive from who and how many rays to send to who
void RayBank::communicate_message_metadata() {  

  // Ensure all values are zero in vector for receiving counts
  fill(num_messages_receiving_.begin(), num_messages_receiving_.end(), 0);
  fill(num_messages_sending_.begin(), num_messages_sending_.end(), 0);

  // Fill the sending counts
  for (auto& [rank, buffers] : ray_send_buffer_) {
    num_messages_sending_[rank] = buffers.count;
  }

  // Exchange message counts with all ranks
  MPI_Alltoall(num_messages_sending_.data(), 1, MPI_INT,
               num_messages_receiving_.data(), 1, MPI_INT,
               mpi::intracomm);

  total_receiving_rays_ = accumulate(num_messages_receiving_.begin(), 
                                     num_messages_receiving_.end(), 0);

}

void RayBank::communicate_rays(){

    // Get the maximum coordinate levels
    const int n_coord_max = model::n_coord_levels;

    // Allocate receiving buffers
    received_ray_data_.resize(total_receiving_rays_);
    received_angular_flux_data_.resize(total_receiving_rays_ * negroups_);
    received_coord_.resize(total_receiving_rays_ * n_coord_max);
    received_cell_last_.resize(total_receiving_rays_ * n_coord_max);

    // Calculate total number of MPI requests needed
    // 4 messages per sending rank + 4 messages per receiving rank
    int num_send_ranks = ray_send_buffer_.size();
    int num_recv_ranks = 0;
    for (int i = 0; i < mpi::n_procs; i++) {
      if (num_messages_receiving_[i] > 0) num_recv_ranks++;
    }
    
    int total_requests = num_send_ranks * 4 + num_recv_ranks * 4;
    vector<MPI_Request> requests(total_requests);
    int req_idx = 0;

    // Post all non-blocking receives first to allow for potential overlap 
    // of communication and packing of send buffers
    int recv_offset = 0;
    for (int sending_rank = 0; sending_rank < mpi::n_procs; sending_rank++) {
      int num_rays_receiving = num_messages_receiving_[sending_rank];
      if (num_rays_receiving == 0) continue;
      
      MPI_Irecv(&received_ray_data_[recv_offset], 
                num_rays_receiving * sizeof(RayExchangeData), 
                MPI_BYTE, sending_rank, 1, mpi::intracomm, &requests[req_idx++]);
      
      MPI_Irecv(&received_angular_flux_data_[recv_offset * negroups_], 
                num_rays_receiving * negroups_, 
                MPI_FLOAT, sending_rank, 2, mpi::intracomm, &requests[req_idx++]);
      
      MPI_Irecv(&received_coord_[recv_offset * n_coord_max], 
                num_rays_receiving * n_coord_max * sizeof(LocalCoord), 
                MPI_BYTE, sending_rank, 3, mpi::intracomm, &requests[req_idx++]);
      
      MPI_Irecv(&received_cell_last_[recv_offset * n_coord_max], 
                num_rays_receiving * n_coord_max, 
                MPI_INT, sending_rank, 4, mpi::intracomm, &requests[req_idx++]);

      recv_offset += num_rays_receiving;
    }

    // Post all non-blocking sends!
    for (auto& [receiving_rank, buffers] : ray_send_buffer_) {
      int num_rays = buffers.count;
      
      // Check neighbor list and add if not already known
      // Filter out rays that are sampled elsewhere
      bool has_transported_rays = false;
      for (int i = 0; i < num_rays; i++) {
        if (buffers.ray_data[i].distance_travelled > 0.0 || buffers.ray_data[i].is_active) {
          has_transported_rays = true;
          break;
        }
      }
      if (has_transported_rays) {
        mpi::decomp_map.my_neighbors_.insert(receiving_rank);
      }
      
      // Send all 4 data arrays - already packed, no intermediate copying!
      MPI_Isend(buffers.ray_data.data(), 
                num_rays * sizeof(RayExchangeData), 
                MPI_BYTE, receiving_rank, 1, mpi::intracomm, &requests[req_idx++]);
      
      MPI_Isend(buffers.angular_flux.data(), 
                num_rays * negroups_, 
                MPI_FLOAT, receiving_rank, 2, mpi::intracomm, &requests[req_idx++]);
      
      MPI_Isend(buffers.coord.data(), 
                num_rays * n_coord_max * sizeof(LocalCoord), 
                MPI_BYTE, receiving_rank, 3, mpi::intracomm, &requests[req_idx++]);
      
      MPI_Isend(buffers.cell_last.data(), 
                num_rays * n_coord_max, 
                MPI_INT, receiving_rank, 4, mpi::intracomm, &requests[req_idx++]);
    }

    // Wait for all communication to complete
    MPI_Waitall(req_idx, requests.data(), MPI_STATUSES_IGNORE);
    
    // Clear send buffers
    ray_send_buffer_.clear();
}

void RayBank::update_my_ray_list(FlatSourceDomain* domain){

  my_ray_list_.resize(received_ray_data_.size());

  const int n_coord_max = model::n_coord_levels;

  // Add re-initialized random ray objects to my_ray_list
  #pragma omp parallel for
    for (int i = 0; i < received_ray_data_.size(); i++) {

      // Re-initialize rays with received data, including full geometry state
      my_ray_list_[i].restart_ray(domain, received_ray_data_[i], 
                                   &received_angular_flux_data_[i * negroups_],
                                   &received_coord_[i * n_coord_max],
                                   &received_cell_last_[i * n_coord_max]);
    }

  // Clear received data vectors
  received_ray_data_.resize(0);
  received_angular_flux_data_.resize(0);
  received_coord_.resize(0);
  received_cell_last_.resize(0);
}

bool RayBank::is_any_ray_alive(){

  int local_rays_alive = ray_bank_size();
  int flag = 0;
  if (local_rays_alive > 0) {
    flag = 1;
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_INT, MPI_MAX, mpi::intracomm);
  
  return flag > 0;
}

}// namespace openmc
