#include "openmc/random_ray/ray_bank.h"

#include "openmc/random_ray/random_ray.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/message_passing.h"
#include "openmc/random_ray/decomposition_map.h"
#include "openmc/mgxs_interface.h"
#include "openmc/timer.h"

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

    // Add ray data to buffer
    ray_send_buffer_[rank].push_back(ray.exchange_data_);
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

// Tells each rank how many rays to receive from whom
void RayBank::communicate_message_metadata() {  

  // Ensure all values are zero in vector for receiving counts
  fill(num_messages_receiving_.begin(), num_messages_receiving_.end(), 0);
  fill(num_messages_sending_.begin(), num_messages_sending_.end(), 0);

  // Fill the sending counts
  for (auto& [rank, rays] : ray_send_buffer_) {
    num_messages_sending_[rank] = rays.size();
  }

  // Exchange message counts with all ranks
  MPI_Alltoall(num_messages_sending_.data(), 1, MPI_INT,
               num_messages_receiving_.data(), 1, MPI_INT,
               mpi::intracomm);

  total_sending_rays_ = accumulate(num_messages_sending_.begin(), 
                                     num_messages_sending_.end(), 0);

  total_receiving_rays_ = accumulate(num_messages_receiving_.begin(), 
                                     num_messages_receiving_.end(), 0);

}

void RayBank::communicate_rays(){

    // Each ray requires 2 sends (data + angular flux)
    int num_requests = ray_send_buffer_.size() * 2;
    vector<MPI_Request> requests(num_requests);
    int req_idx = 0;

    // Define one-dimensional arrays to be sent and received and allocate size
    // Sending
    vector<RayExchangeData> ray_data;
    vector<float> angular_flux_data;
    ray_data.resize(total_sending_rays_);
    angular_flux_data.resize(total_sending_rays_ * negroups_);  

    // Receiving
    received_ray_data_.resize(total_receiving_rays_);
    received_angular_flux_data_.resize(total_receiving_rays_ * negroups_);

    // Indices to keep track of where to pack/unpack data in 1D arrays
    int vector_send_idx = 0;
    int vector_receive_idx = 0;

    // Send ray data to neighbouring ranks
    for (auto [receiving_rank, rays] : ray_send_buffer_) {

      int num_rays_sending = num_messages_sending_[receiving_rank];

      for (int i = 0; i < num_rays_sending; i++) { 
        //TODO: get rid of this! Entire rayExchangeData should be buffered
        //TODO: OMP?
        // Pack slimmed down data container for MPI send
        RayExchangeData exchange_data;
        exchange_data.position = rays[i].position;
        exchange_data.direction = rays[i].direction;
        exchange_data.distance_travelled = rays[i].distance_travelled;
        exchange_data.is_active = rays[i].is_active;
        exchange_data.ray_id = rays[i].ray_id;
        exchange_data.surface = rays[i].surface;
        ray_data[vector_send_idx + i] = exchange_data;
        // Angular flux array
        for (int g = 0; g < negroups_; g++){
          angular_flux_data[(vector_send_idx + i) * negroups_ + g] = rays[i].angular_flux[g];
        }
        // Check neighbor list and add if not already known (insert does this check automatically). 
        // Only check active rays to filter out rays that are sampled in wrong subdomain.
        // TODO: Maybe this is not efficient here. If load balancing constrained to first 
        // 5 batches, maybe this should be moved elsewhere
        if (rays[i].distance_travelled > 0.0 || rays[i].is_active) {
            mpi::decomp_map.my_neighbors.insert(receiving_rank);
        }
      }

      MPI_Isend(&ray_data[vector_send_idx], num_rays_sending * sizeof(RayExchangeData), MPI_BYTE, receiving_rank, 1, mpi::intracomm, &requests[req_idx]);
      MPI_Isend(&angular_flux_data[vector_send_idx * negroups_], num_rays_sending * negroups_, MPI_FLOAT, receiving_rank, 2, mpi::intracomm, &requests[req_idx+1]); 

      vector_send_idx += num_rays_sending;
      req_idx += 2;
      }

    //TODO: Post Irecv before Isend?
    //TODO: OMP?
    // Receive ray data from neighbouring ranks 
    for (int sending_rank = 0; sending_rank < mpi::n_procs; sending_rank++) {
      int num_rays_receiving = num_messages_receiving_[sending_rank];
      if (num_rays_receiving == 0) continue;
      MPI_Recv(&received_ray_data_[vector_receive_idx], num_rays_receiving * sizeof(RayExchangeData), MPI_BYTE, sending_rank, 1, mpi::intracomm, MPI_STATUS_IGNORE);
      MPI_Recv(&received_angular_flux_data_[vector_receive_idx * negroups_], num_rays_receiving * negroups_, MPI_FLOAT, sending_rank, 2, mpi::intracomm, MPI_STATUS_IGNORE);

      vector_receive_idx += num_rays_receiving;
    }

    // Wait for all communication to complete
    MPI_Waitall(num_requests, requests.data(), MPI_STATUSES_IGNORE);
    
    // Empty buffered_ray_data
    ray_send_buffer_.clear();
}

void RayBank::update_my_ray_list(FlatSourceDomain* domain){

  my_ray_list_.resize(received_ray_data_.size());

  // Add rays to ray list by restarting them from received data
  #pragma omp parallel for
    for (int i = 0; i < received_ray_data_.size(); i++) {
      my_ray_list_[i].restart_ray(domain, received_ray_data_[i], &received_angular_flux_data_[i * negroups_]);
    }

  // Clear received data vectors
  received_ray_data_.resize(0);
  received_angular_flux_data_.resize(0);
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
