#ifndef OPENMC_SOA_PARTICLE_H
#define OPENMC_SOA_PARTICLE_H

// Get things like NuBank, ParticleType, etc.
#include "openmc/error.h"
#include "openmc/particle_data.h"

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

extern std::vector<NuclideMicroXS> neutron_xs;
extern std::vector<ElementMicroXS> photon_xs;
extern std::vector<MacroXS> macro_xs;
extern std::vector<int64_t> id;
extern std::vector<ParticleType> type;
extern std::vector<int> n_coord;
extern std::vector<int> cell_instance;
extern std::vector<LocalCoord> coord;
extern std::vector<int> n_coord_last;
extern std::vector<int> cell_last;
extern std::vector<double> E;
extern std::vector<double> E_last;
extern std::vector<int> g;
extern std::vector<int> g_last;
extern std::vector<double> wgt;
extern std::vector<double> mu;
extern std::vector<char> alive; // vector<bool> can't return references
extern std::vector<Position> r_last_current;
extern std::vector<Position> r_last;
extern std::vector<Direction> u_last;
extern std::vector<double> wgt_last;
extern std::vector<double> wgt_absorb;
extern std::vector<char> fission;
extern std::vector<TallyEvent> event;
extern std::vector<int> event_nuclide;
extern std::vector<int> event_mt;
extern std::vector<int> delayed_group;
extern std::vector<int> n_bank;
extern std::vector<int> n_bank_second;
extern std::vector<double> wgt_bank;
extern std::vector<int> n_delayed_bank; // MAX_DELAYED_GROUPS pitch
extern std::vector<int> surface;
extern std::vector<int> cell_born;
extern std::vector<int> material;
extern std::vector<int> material_last;
extern std::vector<BoundaryInfo> boundary;
extern std::vector<double> sqrtkT;
extern std::vector<double> sqrtkT_last;
extern std::vector<int> n_collision;
extern std::vector<char> write_track;
extern std::vector<uint64_t> seeds; // N_STREAMS pitch
extern std::vector<int> stream;

extern std::vector<int> secondary_bank_current_indx;
extern std::vector<ParticleBank>
  secondary_bank; // max_secondary_particles pitch

extern std::vector<int64_t> current_work;
extern std::vector<double> flux_derivs;         // tally_derivs.size() pitch
extern std::vector<FilterMatch> filter_matches; // tally_filters.size() pitch

extern std::vector<int> nu_bank_current_indx;
extern std::vector<NuBank> nu_bank;

extern std::vector<double> keff_tally_absorption;
extern std::vector<double> keff_tally_collision;
extern std::vector<double> keff_tally_tracklength;
extern std::vector<double> keff_tally_leakage;
extern std::vector<char> trace;
extern std::vector<double> collision_distance;
extern std::vector<int> n_event;
#ifdef DAGMC
extern std::vector<moab::DagMC::RayHistory> history;
extern std::vector<Direction> last_dir;
#endif
extern std::vector<int64_t> n_progeny;
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
  ParticleHandle(const int& indx) : p(indx) {}

private:
  const int p; // index of particle in arrays

public:
  //==========================================================================
  // Methods and accessors

  NuclideMicroXS& neutron_xs(const int& i)
  {
    return soa::neutron_xs[p * soa::n_nuclides + i];
  }
  const NuclideMicroXS& neutron_xs(const int& i) const
  {
    return soa::neutron_xs[p * soa::n_nuclides + i];
  }
  ElementMicroXS& photon_xs(const int& i)
  {
    return soa::photon_xs[p * soa::n_elements + i];
  }
  MacroXS& macro_xs() { return soa::macro_xs[p]; }
  const MacroXS& macro_xs() const { return soa::macro_xs[p]; }

  int64_t& id() { return soa::id[p]; }
  const int64_t& id() const { return soa::id[p]; }
  ParticleType& type() { return soa::type[p]; }
  const ParticleType& type() const { return soa::type[p]; }

  int& n_coord() { return soa::n_coord[p]; }
  const int& n_coord() const { return soa::n_coord[p]; }
  int& cell_instance() { return soa::cell_instance[p]; }
  const int& cell_instance() const { return soa::cell_instance[p]; }
  LocalCoord& coord(const int& i)
  {
    return soa::coord[p * soa::n_coord_levels + i];
  }
  const LocalCoord& coord(const int& i) const
  {
    return soa::coord[p * soa::n_coord_levels + i];
  }

  int& n_coord_last() { return soa::n_coord_last[p]; }
  const int& n_coord_last() const { return soa::n_coord_last[p]; }
  int& cell_last(const int& i)
  {
    return soa::cell_last[p * soa::n_coord_levels + i];
  }
  const int& cell_last(const int& i) const
  {
    return soa::cell_last[p * soa::n_coord_levels + i];
  }

  double& E() { return soa::E[p]; }
  const double& E() const { return soa::E[p]; }
  double& E_last() { return soa::E_last[p]; }
  const double& E_last() const { return soa::E_last[p]; }
  int& g() { return soa::g[p]; }
  const int& g() const { return soa::g[p]; }
  int& g_last() { return soa::g_last[p]; }
  const int& g_last() const { return soa::g_last[p]; }

  double& wgt() { return soa::wgt[p]; }
  double& mu() { return soa::mu[p]; }
  const double& mu() const { return soa::mu[p]; }
  char& alive() { return soa::alive[p]; }

  Position& r_last_current() { return soa::r_last_current[p]; }
  const Position& r_last_current() const { return soa::r_last_current[p]; }
  Position& r_last() { return soa::r_last[p]; }
  const Position& r_last() const { return soa::r_last[p]; }
  Position& u_last() { return soa::u_last[p]; }
  const Position& u_last() const { return soa::u_last[p]; }
  double& wgt_last() { return soa::wgt_last[p]; }
  const double& wgt_last() const { return soa::wgt_last[p]; }
  double& wgt_absorb() { return soa::wgt_absorb[p]; }
  const double& wgt_absorb() const { return soa::wgt_absorb[p]; }

  char& fission() { return soa::fission[p]; }
  TallyEvent& event() { return soa::event[p]; }
  const TallyEvent& event() const { return soa::event[p]; }
  int& event_nuclide() { return soa::event_nuclide[p]; }
  const int& event_nuclide() const { return soa::event_nuclide[p]; }
  int& event_mt() { return soa::event_mt[p]; }
  int& delayed_group() { return soa::delayed_group[p]; }

  int& n_bank() { return soa::n_bank[p]; }
  int& n_bank_second() { return soa::n_bank_second[p]; }
  double& wgt_bank() { return soa::wgt_bank[p]; }
  int* n_delayed_bank()
  {
    return soa::n_delayed_bank.data() + p * MAX_DELAYED_GROUPS;
  }
  int& n_delayed_bank(const int& i)
  {
    return soa::n_delayed_bank[p * MAX_DELAYED_GROUPS + i];
  }

  int& surface() { return soa::surface[p]; }
  const int& surface() const { return soa::surface[p]; }
  int& cell_born() { return soa::cell_born[p]; }
  const int& cell_born() const { return soa::cell_born[p]; }
  int& material() { return soa::material[p]; }
  const int& material() const { return soa::material[p]; }
  int& material_last() { return soa::material_last[p]; }

  BoundaryInfo& boundary() { return soa::boundary[p]; }

  double& sqrtkT() { return soa::sqrtkT[p]; }
  const double& sqrtkT() const { return soa::sqrtkT[p]; }
  double& sqrtkT_last() { return soa::sqrtkT_last[p]; }

  int& n_collision() { return soa::n_collision[p]; }
  const int& n_collision() const { return soa::n_collision[p]; }

  char& write_track() { return soa::write_track[p]; }
  uint64_t& seeds(const int& i) { return soa::seeds[p * N_STREAMS + i]; }
  uint64_t* seeds() { return soa::seeds.data() + p * N_STREAMS; }
  int& stream() { return soa::stream[p]; }

  ParticleBank& secondary_bank(const int& i)
  {
    return soa::secondary_bank[p * soa::max_secondary_particles + i];
  }
  int64_t& current_work() { return soa::current_work[p]; }
  const int64_t& current_work() const { return soa::current_work[p]; }
  double& flux_derivs(const int& i)
  {
    return soa::flux_derivs[soa::n_tally_derivs * p + i];
  }
  const double& flux_derivs(const int& i) const
  {
    return soa::flux_derivs[soa::n_tally_derivs * p + i];
  }
  FilterMatch& filter_matches(const int& i)
  {
    return soa::filter_matches[soa::n_tally_filters * p + i];
  }
  FilterMatch* filter_matches()
  {
    return &soa::filter_matches[soa::n_tally_filters * p];
  }
  void reset_filter_matches()
  {
    auto start = soa::n_tally_filters * p;
    for (int i = 0; i < soa::n_tally_filters; ++i) {
      soa::filter_matches[start + i].bins_present_ = false;
    }
  }

  std::vector<std::vector<Position>> tracks()
  {
    fatal_error("tracks cannot be written in SOA mode");
    return std::vector<std::vector<Position>>();
  }

  NuBank& nu_bank(const int& i)
  {
    return soa::nu_bank[soa::max_fissions * p + i];
  }
  NuBank& nu_bank_emplace_back()
  {
    return nu_bank(soa::nu_bank_current_indx[p]++);
    assert(soa::nu_bank_current_indx[p] <= soa::max_fissions);
  }
  NuBank& nu_bank_back()
  {
    assert(soa::nu_bank_current_indx[p] > 0);
    return nu_bank(soa::nu_bank_current_indx[p] - 1);
  }
  void nu_bank_clear() { soa::nu_bank_current_indx[p] = 0; }

  double& keff_tally_absorption() { return soa::keff_tally_absorption[p]; }
  double& keff_tally_collision() { return soa::keff_tally_collision[p]; }
  double& keff_tally_tracklength() { return soa::keff_tally_tracklength[p]; }
  double& keff_tally_leakage() { return soa::keff_tally_leakage[p]; }

  char& trace() { return soa::trace[p]; }
  double& collision_distance() { return soa::collision_distance[p]; }
  int& n_event() { return soa::n_event[p]; }

#ifdef DAGMC
  moab::DagMC::RayHistory& rayhistory() { return soa::history[p]; }
  Direction& last_dir() { return soa::last_dir[p]; }
#endif

  int64_t& n_progeny() { return soa::n_progeny[p]; }

  // Accessors for position in global coordinates
  Position& r() { return soa::coord[p * soa::n_coord_levels].r; }
  const Position& r() const { return soa::coord[p * soa::n_coord_levels].r; }

  // Accessors for position in local coordinates
  Position& r_local()
  {
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].r;
  }
  const Position& r_local() const
  {
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].r;
  }

  // Accessors for direction in global coordinates
  Direction& u() { return soa::coord[p * soa::n_coord_levels].u; }
  const Direction& u() const { return soa::coord[p * soa::n_coord_levels].u; }

  // Accessors for direction in local coordinates
  Direction& u_local()
  {
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].u;
  }
  const Direction& u_local() const
  {
    const auto& n_coord = soa::n_coord[p];
    return soa::coord[p * soa::n_coord_levels + n_coord - 1].u;
  }

  //! Gets the pointer to the particle's current PRN seed
  uint64_t* current_seed()
  {
    const auto& stream = soa::stream[p];
    return soa::seeds.data() + p * N_STREAMS + stream;
  }
  const uint64_t* current_seed() const
  {
    const auto& stream = soa::stream[p];
    return soa::seeds.data() + p * N_STREAMS + stream;
  }

  //! Force recalculation of neutron xs by setting last energy to zero
  void invalidate_neutron_xs()
  {
    for (int i_nuc = 0; i_nuc < soa::n_nuclides; ++i_nuc) {
      soa::neutron_xs[p * soa::n_nuclides + i_nuc].last_E = 0.0;
    }
  }

  //! resets all coordinate levels for the particle
  void clear()
  {
    for (int i_level = 0; i_level < soa::n_coord_levels; ++i_level) {
      soa::coord[p * soa::n_coord_levels + i_level].reset();
    }
    soa::n_coord[p] = 1;
  }

  void zero_delayed_bank()
  {
    for (int i_delayed = 0; i_delayed < MAX_DELAYED_GROUPS; ++i_delayed) {
      soa::n_delayed_bank[p * MAX_DELAYED_GROUPS + i_delayed] = 0.0;
    }
  }

  void zero_flux_derivs()
  {
    for (int i_deriv = 0; i_deriv < soa::n_tally_derivs; ++i_deriv) {
      soa::flux_derivs[p * soa::n_tally_derivs + i_deriv] = 0.0;
    }
  }

  // Methods for manipulating the secondary bank

  int secondary_bank_size() { return soa::secondary_bank_current_indx[p]; }

  void secondary_bank_pop_back()
  {
    int& secondary_indx = soa::secondary_bank_current_indx[p];
    secondary_indx--;
    assert(soa::secondary_bank_current_indx[p] >= 0);
  }

  void secondary_bank_push_back(ParticleBank& site)
  {
    // Check we're not writing into the next particle's space
    assert(soa::secondary_bank_current_indx[p] < soa::max_secondary_particles);

    secondary_bank(soa::secondary_bank_current_indx[p]) = site;
    soa::secondary_bank_current_indx[p]++;
  }

  ParticleBank& secondary_bank_back()
  {
    assert(soa::secondary_bank_current_indx[p] > 0);
    return secondary_bank(soa::secondary_bank_current_indx[p] - 1);
  }

  ParticleBank* secondary_bank_end()
  {
    return &secondary_bank(soa::secondary_bank_current_indx[p]);
  }

  bool secondary_bank_empty()
  {
    return soa::secondary_bank_current_indx[p] == 0;
  }

  ParticleBank& secondary_bank_emplace_back()
  {
    return secondary_bank(soa::secondary_bank_current_indx[p]++);
  }

  // Applies defaults as defined in particle_data.h
  void initialize_values()
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

    soa::secondary_bank_current_indx[p] = 0;
    soa::nu_bank_current_indx[p] = 0;
  }
};

} // namespace openmc

#endif // OPENMC_SOA_PARTICLE_H
