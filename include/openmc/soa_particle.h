#ifndef OPENMC_SOA_PARTICLE_H
#define OPENMC_SOA_PARTICLE_H

// Get things like NuBank, ParticleType, etc.
#include "openmc/error.h"
#include "openmc/particle_data.h"
#include "openmc/vector.h"

#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace openmc {

//============================================================================
//! Defines a handle on structure-of-array particle data layout
//============================================================================

void allocate_soa_data();

namespace soa {

// Max number of secondary particles stored at once by a particle
// Note: may need to change this!
constexpr int max_secondary_particles = 5;

// Max possible number of neutrons emitted by one fission event
constexpr int max_fissions = 5;

// Cached numbers of things required for accessing SOA data
// This lets us not include the headers included by soa_particle.cpp
// in practically every other OpenMC, reducing compile time.
extern int n_nuclides;
extern int n_elements;
extern int n_coord_levels;
extern int n_tally_derivs;
extern int n_tally_filters;

namespace gpu {
extern __constant__ int n_nuclides;
extern __constant__ int n_elements;
extern __constant__ int n_coord_levels;
extern __constant__ int n_tally_derivs;
extern __constant__ int n_tally_filters;
} // namespace gpu

extern vector<NuclideMicroXS> neutron_xs;
extern vector<ElementMicroXS> photon_xs;
extern vector<MacroXS> macro_xs;
extern vector<int64_t> id;
extern vector<ParticleType> type;
extern vector<int> n_coord;
extern vector<int> cell_instance;
extern vector<LocalCoord> coord;
extern vector<int> n_coord_last;
extern vector<int> cell_last;
extern vector<double> E;
extern vector<double> E_last;
extern vector<int> g;
extern vector<int> g_last;
extern vector<double> wgt;
extern vector<double> mu;
extern vector<char> alive; // std::vector<bool> can't return references
extern vector<Position> r_last_current;
extern vector<Position> r_last;
extern vector<Direction> u_last;
extern vector<double> wgt_last;
extern vector<double> wgt_absorb;
extern vector<char> fission;
extern vector<TallyEvent> event;
extern vector<int> event_nuclide;
extern vector<int> event_mt;
extern vector<int> delayed_group;
extern vector<int> n_bank;
extern vector<int> n_bank_second;
extern vector<double> wgt_bank;
extern vector<int> n_delayed_bank; // MAX_DELAYED_GROUPS pitch
extern vector<int> surface;
extern vector<int> cell_born;
extern vector<int> material;
extern vector<int> material_last;
extern vector<BoundaryInfo> boundary;
extern vector<double> sqrtkT;
extern vector<double> sqrtkT_last;
extern vector<int> n_collision;
extern vector<char> write_track;
extern vector<uint64_t> seeds; // N_STREAMS pitch
extern vector<int> stream;

extern vector<int> secondary_bank_current_indx;
extern vector<SourceSite> secondary_bank; // max_secondary_particles pitch

extern vector<int64_t> current_work;
extern vector<double> flux_derivs;         // tally_derivs.size() pitch
extern vector<FilterMatch> filter_matches; // tally_filters.size() pitch

extern vector<int> nu_bank_current_indx;
extern vector<NuBank> nu_bank;

extern vector<double> keff_tally_absorption;
extern vector<double> keff_tally_collision;
extern vector<double> keff_tally_tracklength;
extern vector<double> keff_tally_leakage;
extern vector<char> trace;
extern vector<double> collision_distance;
extern vector<int> n_event;
extern vector<int64_t> n_progeny;

namespace gpu {
extern __constant__ NuclideMicroXS* neutron_xs;
extern __constant__ ElementMicroXS* photon_xs;
extern __constant__ MacroXS* macro_xs;
extern __constant__ int64_t* id;
extern __constant__ ParticleType* type;
extern __constant__ int* n_coord;
extern __constant__ int* cell_instance;
extern __constant__ LocalCoord* coord;
extern __constant__ int* n_coord_last;
extern __constant__ int* cell_last;
extern __constant__ double* E;
extern __constant__ double* E_last;
extern __constant__ int* g;
extern __constant__ int* g_last;
extern __constant__ double* wgt;
extern __constant__ double* mu;
extern __constant__ char* alive; // vector<bool* can't return references
extern __constant__ Position* r_last_current;
extern __constant__ Position* r_last;
extern __constant__ Direction* u_last;
extern __constant__ double* wgt_last;
extern __constant__ double* wgt_absorb;
extern __constant__ char* fission;
extern __constant__ TallyEvent* event;
extern __constant__ int* event_nuclide;
extern __constant__ int* event_mt;
extern __constant__ int* delayed_group;
extern __constant__ int* n_bank;
extern __constant__ int* n_bank_second;
extern __constant__ double* wgt_bank;
extern __constant__ int* n_delayed_bank; // MAX_DELAYED_GROUPS pitch
extern __constant__ int* surface;
extern __constant__ int* cell_born;
extern __constant__ int* material;
extern __constant__ int* material_last;
extern __constant__ BoundaryInfo* boundary;
extern __constant__ double* sqrtkT;
extern __constant__ double* sqrtkT_last;
extern __constant__ int* n_collision;
extern __constant__ char* write_track;
extern __constant__ uint64_t* seeds; // N_STREAMS pitch
extern __constant__ int* stream;
extern __constant__ int* secondary_bank_current_indx;
extern __constant__ SourceSite* secondary_bank; // max_secondary_particles pitch
extern __constant__ int64_t* current_work;
extern __constant__ double* flux_derivs;         // tally_derivs.size() pitch
extern __constant__ FilterMatch* filter_matches; // tally_filters.size() pitch
extern __constant__ int* nu_bank_current_indx;
extern __constant__ NuBank* nu_bank;
extern __constant__ double* keff_tally_absorption;
extern __constant__ double* keff_tally_collision;
extern __constant__ double* keff_tally_tracklength;
extern __constant__ double* keff_tally_leakage;
extern __constant__ char* trace;
extern __constant__ double* collision_distance;
extern __constant__ int* n_event;
extern __constant__ int64_t* n_progeny;
} // namespace gpu
} // namespace soa

inline int thread_idx()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

class ParticleHandle {

public:
  ParticleHandle() : p(thread_idx()) {}
  HD ParticleHandle(const int& indx) : p(indx) {}

private:
  const int p; // index of particle in arrays

public:
  //==========================================================================
  // Methods and accessors

  HD NuclideMicroXS& neutron_xs(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::neutron_xs[p * soa::gpu::n_nuclides + i];
#else
    return soa::neutron_xs[p * soa::n_nuclides + i];
#endif
  }
  const NuclideMicroXS& neutron_xs(const int& i) const
  {
    // TODO experiment with __ldg here and return by value
    // (and elsewhere in this class)
#ifdef __CUDA_ARCH__
    return soa::gpu::neutron_xs[p * soa::gpu::n_nuclides + i];
#else
    return soa::neutron_xs[p * soa::n_nuclides + i];
#endif
  }
  HD ElementMicroXS& photon_xs(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::photon_xs[p * soa::gpu::n_elements + i];
#else
    return soa::photon_xs[p * soa::n_elements + i];
#endif
  }
  HD MacroXS& macro_xs()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::macro_xs[p];
#else
    return soa::macro_xs[p];
#endif
  }
  HD const MacroXS& macro_xs() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::macro_xs[p];
#else
    return soa::macro_xs[p];
#endif
  }

  HD int64_t& id()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::id[p];
#else
    return soa::id[p];
#endif
  }
  HD const int64_t& id() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::id[p];
#else
    return soa::id[p];
#endif
  }
  HD ParticleType& type()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::type[p];
#else
    return soa::type[p];
#endif
  }
  HD const ParticleType& type() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::type[p];
#else
    return soa::type[p];
#endif
  }

  HD int& n_coord()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_coord[p];
#else
    return soa::n_coord[p];
#endif
  }
  HD const int& n_coord() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_coord[p];
#else
    return soa::n_coord[p];
#endif
  }
  HD int& cell_instance()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::cell_instance[p];
#else
    return soa::cell_instance[p];
#endif
  }
  HD const int& cell_instance() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::cell_instance[p];
#else
    return soa::cell_instance[p];
#endif
  }
  HD LocalCoord& coord(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::coord[p * soa::gpu::n_coord_levels + i];
#else
    return soa::coord[p * soa::n_coord_levels + i];
#endif
  }
  HD const LocalCoord& coord(const int& i) const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::coord[p * soa::gpu::n_coord_levels + i];
#else
    return soa::coord[p * soa::n_coord_levels + i];
#endif
  }

  HD int& n_coord_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_coord_last[p];
#else
    return soa::n_coord_last[p];
#endif
  }
  HD const int& n_coord_last() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_coord_last[p];
#else
    return soa::n_coord_last[p];
#endif
  }
  HD int& cell_last(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::cell_last[p * soa::gpu::n_coord_levels + i];
#else
    return soa::cell_last[p * soa::n_coord_levels + i];
#endif
  }
  HD const int& cell_last(const int& i) const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::cell_last[p * soa::gpu::n_coord_levels + i];
#else
    return soa::cell_last[p * soa::n_coord_levels + i];
#endif
  }

  HD double& E()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::E[p];
#else
    return soa::E[p];
#endif
  }
  HD const double& E() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::E[p];
#else
    return soa::E[p];
#endif
  }
  HD double& E_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::E_last[p];
#else
    return soa::E_last[p];
#endif
  }
  HD const double& E_last() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::E_last[p];
#else
    return soa::E_last[p];
#endif
  }
  HD int& g()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::g[p];
#else
    return soa::g[p];
#endif
  }
  HD const int& g() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::g[p];
#else
    return soa::g[p];
#endif
  }
  HD int& g_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::g_last[p];
#else
    return soa::g_last[p];
#endif
  }
  HD const int& g_last() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::g_last[p];
#else
    return soa::g_last[p];
#endif
  }

  HD double& wgt()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::wgt[p];
#else
    return soa::wgt[p];
#endif
  }
  HD double& mu()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::mu[p];
#else
    return soa::mu[p];
#endif
  }
  HD const double& mu() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::mu[p];
#else
    return soa::mu[p];
#endif
  }
  HD char& alive()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::alive[p];
#else
    return soa::alive[p];
#endif
  }

  HD Position& r_last_current()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::r_last_current[p];
#else
    return soa::r_last_current[p];
#endif
  }
  HD const Position& r_last_current() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::r_last_current[p];
#else
    return soa::r_last_current[p];
#endif
  }
  HD Position& r_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::r_last[p];
#else
    return soa::r_last[p];
#endif
  }
  HD const Position& r_last() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::r_last[p];
#else
    return soa::r_last[p];
#endif
  }
  HD Position& u_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::u_last[p];
#else
    return soa::u_last[p];
#endif
  }
  HD const Position& u_last() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::u_last[p];
#else
    return soa::u_last[p];
#endif
  }
  HD double& wgt_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::wgt_last[p];
#else
    return soa::wgt_last[p];
#endif
  }
  HD const double& wgt_last() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::wgt_last[p];
#else
    return soa::wgt_last[p];
#endif
  }
  HD double& wgt_absorb()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::wgt_absorb[p];
#else
    return soa::wgt_absorb[p];
#endif
  }
  HD const double& wgt_absorb() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::wgt_absorb[p];
#else
    return soa::wgt_absorb[p];
#endif
  }

  HD char& fission()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::fission[p];
#else
    return soa::fission[p];
#endif
  }
  HD TallyEvent& event()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::event[p];
#else
    return soa::event[p];
#endif
  }
  HD const TallyEvent& event() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::event[p];
#else
    return soa::event[p];
#endif
  }
  HD int& event_nuclide()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::event_nuclide[p];
#else
    return soa::event_nuclide[p];
#endif
  }
  HD const int& event_nuclide() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::event_nuclide[p];
#else
    return soa::event_nuclide[p];
#endif
  }
  HD int& event_mt()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::event_mt[p];
#else
    return soa::event_mt[p];
#endif
  }
  HD int& delayed_group()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::delayed_group[p];
#else
    return soa::delayed_group[p];
#endif
  }

  HD int& n_bank()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_bank[p];
#else
    return soa::n_bank[p];
#endif
  }
  HD int& n_bank_second()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_bank_second[p];
#else
    return soa::n_bank_second[p];
#endif
  }
  HD double& wgt_bank()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::wgt_bank[p];
#else
    return soa::wgt_bank[p];
#endif
  }
  HD int* n_delayed_bank()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_delayed_bank + p * MAX_DELAYED_GROUPS;
#else
    return soa::n_delayed_bank.data() + p * MAX_DELAYED_GROUPS;
#endif
  }
  HD int& n_delayed_bank(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_delayed_bank[p * MAX_DELAYED_GROUPS + i];
#else
    return soa::n_delayed_bank[p * MAX_DELAYED_GROUPS + i];
#endif
  }

  HD int& surface()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::surface[p];
#else
    return soa::surface[p];
#endif
  }
  HD const int& surface() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::surface[p];
#else
    return soa::surface[p];
#endif
  }
  HD int& cell_born()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::cell_born[p];
#else
    return soa::cell_born[p];
#endif
  }
  HD const int& cell_born() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::cell_born[p];
#else
    return soa::cell_born[p];
#endif
  }
  HD int& material()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::material[p];
#else
    return soa::material[p];
#endif
  }
  HD const int& material() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::material[p];
#else
    return soa::material[p];
#endif
  }
  HD int& material_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::material_last[p];
#else
    return soa::material_last[p];
#endif
  }

  HD BoundaryInfo& boundary()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::boundary[p];
#else
    return soa::boundary[p];
#endif
  }

  HD double& sqrtkT()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::sqrtkT[p];
#else
    return soa::sqrtkT[p];
#endif
  }
  HD const double& sqrtkT() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::sqrtkT[p];
#else
    return soa::sqrtkT[p];
#endif
  }
  HD double& sqrtkT_last()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::sqrtkT_last[p];
#else
    return soa::sqrtkT_last[p];
#endif
  }

  HD int& n_collision()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_collision[p];
#else
    return soa::n_collision[p];
#endif
  }
  HD const int& n_collision() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_collision[p];
#else
    return soa::n_collision[p];
#endif
  }

  HD char& write_track()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::write_track[p];
#else
    return soa::write_track[p];
#endif
  }
  HD uint64_t& seeds(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::seeds[p * N_STREAMS + i];
#else
    return soa::seeds[p * N_STREAMS + i];
#endif
  }
  HD uint64_t* seeds()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::seeds + p * N_STREAMS;
#else
    return soa::seeds.data() + p * N_STREAMS;
#endif
  }
  HD int& stream()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::stream[p];
#else
    return soa::stream[p];
#endif
  }

  HD SourceSite& secondary_bank(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::secondary_bank[p * soa::max_secondary_particles + i];
#else
    return soa::secondary_bank[p * soa::max_secondary_particles + i];
#endif
  }
  HD int64_t& current_work()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::current_work[p];
#else
    return soa::current_work[p];
#endif
  }
  HD const int64_t& current_work() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::current_work[p];
#else
    return soa::current_work[p];
#endif
  }
  HD double& flux_derivs(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::flux_derivs[soa::gpu::n_tally_derivs * p + i];
#else
    return soa::flux_derivs[soa::n_tally_derivs * p + i];
#endif
  }
  HD const double& flux_derivs(const int& i) const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::flux_derivs[soa::gpu::n_tally_derivs * p + i];
#else
    return soa::flux_derivs[soa::n_tally_derivs * p + i];
#endif
  }
  HD FilterMatch& filter_matches(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::filter_matches[soa::gpu::n_tally_filters * p + i];
#else
    return soa::filter_matches[soa::n_tally_filters * p + i];
#endif
  }
  HD FilterMatch* filter_matches()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::filter_matches + soa::gpu::n_tally_filters * p;
#else
    return &soa::filter_matches[soa::n_tally_filters * p];
#endif
  }
  HD void reset_filter_matches()
  {
#ifdef __CUDA_ARCH__
    auto start = soa::gpu::n_tally_filters * p;
    for (int i = 0; i < soa::gpu::n_tally_filters; ++i) {
      soa::gpu::filter_matches[start + i].bins_present_ = false;
    }
#else
    auto start = soa::n_tally_filters * p;
    for (int i = 0; i < soa::n_tally_filters; ++i) {
      soa::filter_matches[start + i].bins_present_ = false;
    }
#endif
  }

  vector<vector<Position>> tracks()
  {
    fatal_error("tracks cannot be written in SOA mode");
    return vector<vector<Position>>();
  }

  HD NuBank& nu_bank(const int& i)
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::nu_bank[soa::max_fissions * p + i];
#else
    return soa::nu_bank[soa::max_fissions * p + i];
#endif
  }
  HD NuBank& nu_bank_emplace_back()
  {
#ifdef __CUDA_ARCH__
    return nu_bank(soa::gpu::nu_bank_current_indx[p]++);
#else
    return nu_bank(soa::nu_bank_current_indx[p]++);
    assert(soa::nu_bank_current_indx[p] <= soa::max_fissions);
#endif
  }
  HD NuBank& nu_bank_back()
  {
#ifdef __CUDA_ARCH__
    return nu_bank(soa::gpu::nu_bank_current_indx[p] - 1);
#else
    return nu_bank(soa::nu_bank_current_indx[p] - 1);
#endif
  }
  HD void nu_bank_clear()
  {
#ifdef __CUDA_ARCH__
    soa::gpu::nu_bank_current_indx[p] = 0;
#else
    soa::nu_bank_current_indx[p] = 0;
#endif
  }

  HD double& keff_tally_absorption()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::keff_tally_absorption[p];
#else
    return soa::keff_tally_absorption[p];
#endif
  }
  HD double& keff_tally_collision()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::keff_tally_collision[p];
#else
    return soa::keff_tally_collision[p];
#endif
  }
  HD double& keff_tally_tracklength()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::keff_tally_tracklength[p];
#else
    return soa::keff_tally_tracklength[p];
#endif
  }
  HD double& keff_tally_leakage()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::keff_tally_leakage[p];
#else
    return soa::keff_tally_leakage[p];
#endif
  }

  HD char& trace()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::trace[p];
#else
    return soa::trace[p];
#endif
  }
  HD double& collision_distance()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::collision_distance[p];
#else
    return soa::collision_distance[p];
#endif
  }
  HD int& n_event()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_event[p];
#else
    return soa::n_event[p];
#endif
  }

  HD int64_t& n_progeny()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::n_progeny[p];
#else
    return soa::n_progeny[p];
#endif
  }

  // Accessors for position in global coordinates
  HD Position& r()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::coord[p * soa::gpu::n_coord_levels].r;
#else
    return soa::coord[p * soa::n_coord_levels].r;
#endif
  }
  HD const Position& r() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::coord[p * soa::gpu::n_coord_levels].r;
#else
    return soa::coord[p * soa::n_coord_levels].r;
#endif
  }

  // Accessors for position in local coordinates
  HD Position& r_local()
  {
#ifdef __CUDA_ARCH__
    const auto& n_coord = soa::gpu::n_coord[p];
    return soa::gpu::coord[p * soa::gpu::n_coord_levels + n_coord - 1].r;
#else
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].r;
#endif
  }
  const Position& r_local() const
  {
#ifdef __CUDA_ARCH__
    const auto& n_coord = soa::gpu::n_coord[p];
    return soa::gpu::coord[p * soa::gpu::n_coord_levels + n_coord - 1].r;
#else
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].r;
#endif
  }

  // Accessors for direction in global coordinates
  HD Direction& u()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::coord[p * soa::gpu::n_coord_levels].u;
#else
    return soa::coord[p * soa::n_coord_levels].u;
#endif
  }
  HD const Direction& u() const
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::coord[p * soa::gpu::n_coord_levels].u;
#else
    return soa::coord[p * soa::n_coord_levels].u;
#endif
  }

  // Accessors for direction in local coordinates
  HD Direction& u_local()
  {
#ifdef __CUDA_ARCH__
    const auto& n_coord = soa::gpu::n_coord[p];
    return soa::gpu::coord[p * soa::gpu::n_coord_levels + n_coord - 1].u;
#else
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].u;
#endif
  }
  HD const Direction& u_local() const
  {
#ifdef __CUDA_ARCH__
    const auto& n_coord = soa::gpu::n_coord[p];
    return soa::gpu::coord[p * soa::gpu::n_coord_levels + n_coord - 1].u;
#else
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].u;
#endif
  }

  //! Gets the pointer to the particle's current PRN seed
  HD uint64_t* current_seed()
  {
#ifdef __CUDA_ARCH__
    const auto& stream = soa::gpu::stream[p];
    return soa::gpu::seeds + p * N_STREAMS + stream;
#else
    const auto& stream = soa::stream[p];
    return soa::seeds.data() + p * N_STREAMS + stream;
#endif
  }
  HD const uint64_t* current_seed() const
  {
#ifdef __CUDA_ARCH__
    const auto& stream = soa::gpu::stream[p];
    return soa::gpu::seeds + p * N_STREAMS + stream;
#else
    const auto& stream = soa::stream[p];
    return soa::seeds.data() + p * N_STREAMS + stream;
#endif
  }

  //! Force recalculation of neutron xs by setting last energy to zero
  HD void invalidate_neutron_xs()
  {
#ifdef __CUDA_ARCH__
    for (int i_nuc = 0; i_nuc < soa::gpu::n_nuclides; ++i_nuc) {
      soa::gpu::neutron_xs[p * soa::gpu::n_nuclides + i_nuc].last_E = 0.0;
    }
#else
    for (int i_nuc = 0; i_nuc < soa::n_nuclides; ++i_nuc) {
      soa::neutron_xs[p * soa::n_nuclides + i_nuc].last_E = 0.0;
    }
#endif
  }

  //! resets all coordinate levels for the particle
  HD void clear()
  {
#ifdef __CUDA_ARCH__
    for (int i_level = 0; i_level < soa::gpu::n_coord_levels; ++i_level) {
      soa::gpu::coord[p * soa::gpu::n_coord_levels + i_level].reset();
    }
    soa::gpu::n_coord[p] = 1;
#else
    for (int i_level = 0; i_level < soa::n_coord_levels; ++i_level) {
      soa::coord[p * soa::n_coord_levels + i_level].reset();
    }
#endif
  }

  HD void zero_delayed_bank()
  {
#ifdef __CUDA_ARCH__
    for (int i_delayed = 0; i_delayed < MAX_DELAYED_GROUPS; ++i_delayed) {
      soa::gpu::n_delayed_bank[p * MAX_DELAYED_GROUPS + i_delayed] = 0.0;
    }
#else
    for (int i_delayed = 0; i_delayed < MAX_DELAYED_GROUPS; ++i_delayed) {
      soa::n_delayed_bank[p * MAX_DELAYED_GROUPS + i_delayed] = 0.0;
    }
#endif
  }

  HD void zero_flux_derivs()
  {
#ifdef __CUDA_ARCH__
    for (int i_deriv = 0; i_deriv < soa::gpu::n_tally_derivs; ++i_deriv) {
      soa::gpu::flux_derivs[p * soa::gpu::n_tally_derivs + i_deriv] = 0.0;
    }
#else
    for (int i_deriv = 0; i_deriv < soa::n_tally_derivs; ++i_deriv) {
      soa::flux_derivs[p * soa::n_tally_derivs + i_deriv] = 0.0;
    }
#endif
  }

  // Methods for manipulating the secondary bank

  HD int secondary_bank_size()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::secondary_bank_current_indx[p];
#else
    return soa::secondary_bank_current_indx[p];
#endif
  }

  HD void secondary_bank_pop_back()
  {
#ifdef __CUDA_ARCH__
    int& secondary_indx = soa::gpu::secondary_bank_current_indx[p];
    secondary_indx--;
    assert(soa::gpu::secondary_bank_current_indx[p] >= 0);
#else
    int& secondary_indx = soa::secondary_bank_current_indx[p];
    secondary_indx--;
    assert(soa::secondary_bank_current_indx[p] >= 0);
#endif
  }

  HD void secondary_bank_push_back(SourceSite& site)
  {
#ifdef __CUDA_ARCH__
    // Check we're not writing into the next particle's space
    assert(
      soa::gpu::secondary_bank_current_indx[p] < soa::max_secondary_particles);

    secondary_bank(soa::gpu::secondary_bank_current_indx[p]) = site;
    soa::gpu::secondary_bank_current_indx[p]++;
#else
    assert(soa::secondary_bank_current_indx[p] < soa::max_secondary_particles);

    secondary_bank(soa::secondary_bank_current_indx[p]) = site;
    soa::secondary_bank_current_indx[p]++;
#endif
  }

  HD SourceSite& secondary_bank_back()
  {
#ifdef __CUDA_ARCH__
    assert(soa::gpu::secondary_bank_current_indx[p] > 0);
    return secondary_bank(soa::gpu::secondary_bank_current_indx[p] - 1);
#else
    assert(soa::secondary_bank_current_indx[p] > 0);
    return secondary_bank(soa::secondary_bank_current_indx[p] - 1);
#endif
  }

  HD SourceSite* secondary_bank_end()
  {
#ifdef __CUDA_ARCH__
    return &secondary_bank(soa::gpu::secondary_bank_current_indx[p]);
#else
    return &secondary_bank(soa::secondary_bank_current_indx[p]);
#endif
  }

  HD bool secondary_bank_empty()
  {
#ifdef __CUDA_ARCH__
    return soa::gpu::secondary_bank_current_indx[p] == 0;
#else
    return soa::secondary_bank_current_indx[p] == 0;
#endif
  }

  HD SourceSite& secondary_bank_emplace_back()
  {
#ifdef __CUDA_ARCH__
    return secondary_bank(soa::gpu::secondary_bank_current_indx[p]++);
#else
    return secondary_bank(soa::secondary_bank_current_indx[p]++);
#endif
  }

  // Applies defaults as defined in particle_data.h
  HD void initialize_values()
  {
    type() = ParticleType::neutron;
    n_coord() = 1;
    n_coord_last() = 1;
    g() = 0;
    wgt() = 1.0;
    alive() = true;
    wgt_last() = 1.0;
    wgt_absorb() = 1.0;
    fission() = false;
    delayed_group() = 0;
    n_bank() = 0;
    n_bank_second() = 0;
    wgt_bank() = 0.0;
    surface() = 0;
    cell_born() = -1;
    material() = -1;
    material_last() = -1;
    sqrtkT() = -1.0;
    sqrtkT_last() = 0.0;
    n_collision() = 0;
    write_track() = false;
    keff_tally_absorption() = 0.0;
    keff_tally_collision() = 0.0;
    keff_tally_tracklength() = 0.0;
    keff_tally_leakage() = 0.0;
    trace() = false;
    n_event() = 0;
    n_progeny() = 0;

#ifdef __CUDA_ARCH__
    soa::gpu::secondary_bank_current_indx[p] = 0;
    soa::gpu::nu_bank_current_indx[p] = 0;
#else
    soa::secondary_bank_current_indx[p] = 0;
    soa::nu_bank_current_indx[p] = 0;
#endif
  }
};

} // namespace openmc

#endif // OPENMC_SOA_PARTICLE_H
