#ifndef OPENMC_DECOMPOSITION_MAP_H
#define OPENMC_DECOMPOSITION_MAP_H

#include "openmc/vector.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {

class DecompositionMap;

namespace mpi {
    extern DecompositionMap decomp_map;
} // namespace mpi

class DecompositionMap {
public:
  //----------------------------------------------------------------------------
  // Constructors
  DecompositionMap();

  //----------------------------------------------------------------------------
  // Methods
  // void initialize(FlatSourceDomain* domain);
  void update();
  void create(FlatSourceDomain* domain);
  bool is_SRK_in_domain(int sr_key);
  // bool is_SRK_in_domain(SourceRegionKey sr_key);
  int n_source_regions();

  //----------------------------------------------------------------------------
  // Static data members


  //----------------------------------------------------------------------------
  // Public data members

  // Map that relates a SourceRegionKey to the index of the MPI rank that
  // contains that source region in its subdomain.
  // std::unordered_map<SourceRegionKey, int, SourceRegionKey::HashFunctor>
  //   subdomain_map_;

  std::unordered_map<int, int> subdomain_map_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  // vector<int> my_subdomain_list_;

}; // class DecompositionMap

} // namespace openmc

#endif // OPENMC_DECOMPOSITION_MAP_H
