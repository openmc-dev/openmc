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
}

// Initialize list of ray bank that each MPI rank will handle
void RayBank::add_ray_to_bank(RandomRay& ray){
    my_ray_list_.push_back(ray);
}

// Buffer ray that has left my subdomain
void RayBank::buffer_ray_data_to_send(RandomRay& ray, FlatSourceDomain* domain){

    // Get rank to send ray to
    int rank = ray.owner_rank_;

    // printf("RANK %d: Put ray %ld into map, new owner: %d, %d \n", mpi::rank, ray.exchange_data_.ray_id, ray.owner_rank_, rank);

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

    // for (auto& [rank, rays] : ray_send_buffer_) {
    //   for (int r = 0; r < rays.size(); r++){
    //     printf("RANK %d: 2: Buffered ray ID %ld to be sent to rank: %d \n", mpi::rank, rays[r].ray_id, rank);
    //   }
    // }

    // Empty ray list because rays have either died or are in buffer to be sent to other ranks
    reset_my_ray_list();

    // for (auto& [rank, rays] : ray_send_buffer_) {
    //   for (int r = 0; r < rays.size(); r++){
    //     printf("RANK %d: 3: Buffered ray ID %ld to be sent to rank: %d \n", mpi::rank, rays[r].ray_id, rank);
    //   }
    // }

    // Communicate how many rays will be sent to each rank

    

    simulation::time_mpi_imbalance.start();
    MPI_Barrier(mpi::intracomm);
    simulation::time_mpi_imbalance.stop();

    simulation::time_comms_metadata.start();
    communicate_message_metadata();
    simulation::time_comms_metadata.stop();

    // Send and receive rays between MPI ranks
    simulation::time_ray_comms.start();
    communicate_rays();
    simulation::time_ray_comms.stop();


    simulation::time_unpack_data.start();
    // Add received rays to ray list of that rank
    update_my_ray_list(domain);
    simulation::time_unpack_data.stop();

}

// Clears my_ray_list, but keeps memory allocation in place
void RayBank::reset_my_ray_list(){
  my_ray_list_.resize(0);
}

int RayBank::ray_bank_size(){
  int ray_bank_size = my_ray_list_.size();
  return ray_bank_size;
}

// Tells each rank how many rays to receive from whom
void RayBank::communicate_message_metadata() {  
  vector<int> num_messages_sending(mpi::n_procs, 0);

  // Ensure all values are zero
  fill(num_messages_receiving_.begin(), num_messages_receiving_.end(), 0);

  total_sending_rays_ = 0;

  // Fill the sending counts //TODO: OMP?
  for (auto& [rank, rays] : ray_send_buffer_) {
    num_messages_sending[rank] = rays.size();
    total_sending_rays_ += num_messages_sending[rank];
    // for (int r = 0; r < rays.size(); r++){
    //   printf("RANK %d: Metadata, ray ID %ld to be sent to rank: %d \n", mpi::rank, rays[r].ray_id, rank);
    // }
  }

  // Exchange message counts with all ranks
  MPI_Alltoall(num_messages_sending.data(), 1, MPI_INT,
               num_messages_receiving_.data(), 1, MPI_INT,
               mpi::intracomm);

  total_receiving_rays_ = accumulate(num_messages_receiving_.begin(), 
                                     num_messages_receiving_.end(), 0);

  // // Store results
  // for (int rank = 0; rank < mpi::n_procs; rank++) {
  //   // Skip over non-sending ranks
  //   if (num_messages_receiving[rank] == 0) continue;
  //   // Save number of incoming rays from each rank sending to me in incoming_ray_data map
  //   incoming_ray_data_[rank] = num_messages_receiving[rank]; //TODO: Could also just be vector
  //   total_receiving_rays_ += num_messages_receiving[rank];
  // }

}

void RayBank::communicate_rays(){

    // Each ray requires 2 sends (data + angular flux)
    int num_requests = ray_send_buffer_.size() * 2;
    vector<MPI_Request> requests(num_requests); // heap
    // MPI_Request requests[num_requests]; // stack
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

      int num_rays_sending = rays.size();

      for (int i = 0; i < num_rays_sending; i++) {
        // Pack slimmed down data container for MPI send
        RayExchangeData exchange_data;
        exchange_data.position = rays[i].position;
        exchange_data.direction = rays[i].direction;
        exchange_data.distance_travelled = rays[i].distance_travelled;
        exchange_data.is_active = rays[i].is_active;
        exchange_data.ray_id = rays[i].ray_id;
        ray_data[vector_send_idx + i] = exchange_data;
        // Angular flux array
        for (int g = 0; g < negroups_; g++){
          angular_flux_data[(vector_send_idx + i) * negroups_ + g] = rays[i].angular_flux[g];
        }
      }
      // printf("RANK %d: Sending ray %ld to rank: %d \n", mpi::rank, ray_data[vector_send_idx].ray_id, receiving_rank);

      MPI_Isend(&ray_data[vector_send_idx], num_rays_sending * sizeof(RayExchangeData), MPI_BYTE, receiving_rank, 1, mpi::intracomm, &requests[req_idx]);
      MPI_Isend(&angular_flux_data[vector_send_idx * negroups_], num_rays_sending * negroups_, MPI_FLOAT, receiving_rank, 2, mpi::intracomm, &requests[req_idx+1]); 

      vector_send_idx += num_rays_sending;
      req_idx += 2;
      }

    //TODO: Post Irecv before Isend?
    // Receive ray data from neighbouring ranks //TODO: OMP?
    for (int sending_rank = 0; sending_rank < mpi::n_procs; sending_rank++) {
      int num_rays_receiving = num_messages_receiving_[sending_rank];
      if (num_rays_receiving == 0) continue;
      MPI_Recv(&received_ray_data_[vector_receive_idx], num_rays_receiving * sizeof(RayExchangeData), MPI_BYTE, sending_rank, 1, mpi::intracomm, MPI_STATUS_IGNORE);
      MPI_Recv(&received_angular_flux_data_[vector_receive_idx * negroups_], num_rays_receiving * negroups_, MPI_FLOAT, sending_rank, 2, mpi::intracomm, MPI_STATUS_IGNORE);
      // printf("RANK %d: Receiving ray %ld from rank: %d \n", mpi::rank, received_ray_data_[vector_receive_idx].ray_id, sending_rank);

      vector_receive_idx += num_rays_receiving;
    }

    // Wait for all communication to complete
    // MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE);
    MPI_Waitall(num_requests, requests.data(), MPI_STATUSES_IGNORE);
    
    // Empty buffered_ray_data
    ray_send_buffer_.clear();
}

void RayBank::update_my_ray_list(FlatSourceDomain* domain){

  //TODO: OMP?

  // Temporary vector containing angular flux data for re-initialization of random rays
  vector <float> angular_flux_vec;
  angular_flux_vec.resize(negroups_);

  // Add re-initialized random ray objects to my_ray_list
  for (int i = 0; i < received_ray_data_.size(); i++) {

    for (int g = 0; g < negroups_; g++) {
      angular_flux_vec[g] = received_angular_flux_data_[i * negroups_ + g];
    }

    // Re-initialize rays with received data
    RandomRay ray_received(domain, received_ray_data_[i], angular_flux_vec);
    my_ray_list_.push_back(ray_received);
  }

  // clear received data vectors //TODO: maybe resize(0)?
  received_ray_data_.clear();
  received_angular_flux_data_.clear();

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
