#ifndef OPENMC_RAY_BANK_H
#define OPENMC_RAY_BANK_H

#include "openmc/vector.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/random_ray/flat_source_domain.h"
#ifdef OPENMC_MPI
#include <mpi.h>
#endif

namespace openmc {

// // Forward declaration
// class FlatSourceDomain;
// class RandomRay;
// struct RayBufferContainer;

class RayBank {
public:
  //----------------------------------------------------------------------------
  // Constructors
  RayBank();

  //----------------------------------------------------------------------------
  // Methods
  void buffer_ray_data_to_send(RandomRay& ray, FlatSourceDomain* domain);
  void update(FlatSourceDomain* domain);
  int ray_bank_size();
  void communicate_rays();
  void communicate_message_metadata();
  void update_my_ray_list(FlatSourceDomain* domain);
  bool is_any_ray_alive(); 

  //----------------------------------------------------------------------------
  // Public data members
  vector<RandomRay> my_ray_list_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  int total_receiving_rays_;
  int negroups_;

  // Per-rank send buffers for ray data, including geometry state and angular flux
  struct RankSendBuffers {
    vector<RayExchangeData> ray_data;
    vector<float> angular_flux;
    vector<LocalCoord> coord;
    vector<int> cell_last;
    int count = 0;  // Number of rays buffered for this rank
  };
  std::unordered_map<int, RankSendBuffers> ray_send_buffer_;
  int reserved_buffer_size_ = 32; // Initial reserved size for send buffers, can be tuned based on expected ray counts

  // Vector that contains the number of rays to be received from each rank
  vector<int> num_messages_receiving_;
  vector<int> num_messages_sending_;

  // vectors that received ray data
  vector<RayExchangeData> received_ray_data_;
  vector<float> received_angular_flux_data_;
  vector<LocalCoord> received_coord_;
  vector<int> received_cell_last_;

}; // class DecompositionMap

} // namespace openmc

#endif // OPENMC_RAY_BANK_H
