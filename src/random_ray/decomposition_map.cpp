#include "openmc/random_ray/decomposition_map.h"

#include "openmc/message_passing.h"
#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"

namespace openmc {

namespace mpi {
    DecompositionMap decomp_map;
}

// Constructor
DecompositionMap::DecompositionMap() {}

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

int DecompositionMap::n_source_regions(){
  int count = 0;
  
  for (const auto& [key, val] : subdomain_map_) {
    if (val == mpi::rank) {
      count++;
    }
  }

  return count;
}

}// namespace openmc
