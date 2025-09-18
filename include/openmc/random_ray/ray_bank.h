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
  void add_ray_to_bank(RandomRay& ray);
  void buffer_ray_data_to_send(RandomRay& ray, FlatSourceDomain* domain);
  void update(FlatSourceDomain* domain);
  int ray_bank_size();
  void reset_my_ray_list();
  void communicate_rays();
  void communicate_message_metadata();
  void update_my_ray_list(FlatSourceDomain* domain);
  bool is_any_ray_alive(); 

  //----------------------------------------------------------------------------
  // Static data members


  //----------------------------------------------------------------------------
  // Public data members
  vector<RandomRay> my_ray_list_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  int total_sending_rays_;
  int total_receiving_rays_;
  int negroups_;
  // FlatSourceDomain* domain_ {nullptr}; // pointer to domain that has flat source //TODO: maybe include domain in constructor
  //                                      // data needed for ray transport

  // Map that contains the rank to which rays are buffered to be sent
  std::unordered_map<int, vector<RayBufferContainer>> ray_send_buffer_;

  // Map that contains the number of rays to be received from each rank
  vector<int> num_messages_receiving_;

  // vectors that received ray data
  vector<RayExchangeData> received_ray_data_;
  vector<float> received_angular_flux_data_;

}; // class DecompositionMap

} // namespace openmc

#endif // OPENMC_RAY_BANK_H
