#include "openmc/soa_particle.h"
#include "openmc/settings.h"
namespace openmc {

void allocate_soa_data()
{
  using namespace soa;

  // Pick the number of particles in flight
  int particles_in_flight;
  if (settings::event_based) {
    particles_in_flight = settings::max_particles_in_flight;
  } else {
#ifdef _OPENMP
    particles_in_flight = omp_get_max_threads();
#else
    particles_in_flight = 1;
#endif
  }

  // Allocate sufficient room in all arrays
  neutron_xs.resize(particles_in_flight);
  photon_xs.resize(particles_in_flight);
  macro_xs.resize(particles_in_flight);
  id.resize(particles_in_flight);
  type.resize(particles_in_flight);
  n_coord.resize(particles_in_flight);
  cell_instance.resize(particles_in_flight);
  coord.resize(particles_in_flight);
  n_coord_last.resize(particles_in_flight);
  cell_last.resize(particles_in_flight);
  E.resize(particles_in_flight);
  E_last.resize(particles_in_flight);
  g.resize(particles_in_flight);
  g_last.resize(particles_in_flight);
  wgt.resize(particles_in_flight);
  mu.resize(particles_in_flight);
  alive.resize(particles_in_flight);
  r_last_current.resize(particles_in_flight);
  r_last.resize(particles_in_flight);
  u_last.resize(particles_in_flight);
  wgt_last.resize(particles_in_flight);
  wgt_absorb.resize(particles_in_flight);
  fission.resize(particles_in_flight);
  event.resize(particles_in_flight);
  event_nuclide.resize(particles_in_flight);
  event_mt.resize(particles_in_flight);
  delayed_group.resize(particles_in_flight);
  n_bank.resize(particles_in_flight);
  n_bank_second.resize(particles_in_flight);
  wgt_bank.resize(particles_in_flight);
  n_delayed_bank.resize(particles_in_flight * MAX_DELAYED_GROUPS);
  surface.resize(particles_in_flight);
  cell_born.resize(particles_in_flight);
  material.resize(particles_in_flight);
  material_last.resize(particles_in_flight);
  boundary.resize(particles_in_flight);
  sqrtkT.resize(particles_in_flight);
  sqrtkT_last.resize(particles_in_flight);
  n_collision.resize(particles_in_flight);
  write_track.resize(particles_in_flight);
  seeds.resize(N_STREAMS * particles_in_flight);
  stream.resize(particles_in_flight);
  secondary_bank.resize(max_secondary_particles * particles_in_flight);
  current_work.resize(particles_in_flight);
  flux_derivs.resize(model::tally_derivs.size() * particles_in_flight);
  filter_matches.resize(model::tally_filters.size() * particles_in_flight);
  nu_bank.resize(max_fissions * particles_in_flight);
  keff_tally_absorption.resize(particles_in_flight);
  keff_tally_collision.resize(particles_in_flight);
  keff_tally_tracklength.resize(particles_in_flight);
  keff_tally_leakage.resize(particles_in_flight);
  trace.resize(particles_in_flight);
  collision_distance.resize(particles_in_flight);
  n_event.resize(particles_in_flight);
#ifdef DAGMC
  history.resize(particles_in_flight);
  last_dir.resize(particles_in_flight);
#endif
  n_progeny.resize(particles_in_flight);
}

namespace soa {
std::vector<NuclideMicroXS> neutron_xs;
std::vector<ElementMicroXS> photon_xs;
std::vector<MacroXS> macro_xs;
std::vector<int64_t> id;
std::vector<ParticleType> type;
std::vector<int> n_coord;
std::vector<int> cell_instance;
std::vector<LocalCoord> coord;
std::vector<int> n_coord_last;
std::vector<int> cell_last;
std::vector<double> E;
std::vector<double> E_last;
std::vector<int> g;
std::vector<int> g_last;
std::vector<double> wgt;
std::vector<double> mu;
std::vector<char> alive;
std::vector<Position> r_last_current;
std::vector<Position> r_last;
std::vector<Direction> u_last;
std::vector<double> wgt_last;
std::vector<double> wgt_absorb;
std::vector<char> fission;
std::vector<TallyEvent> event;
std::vector<int> event_nuclide;
std::vector<int> event_mt;
std::vector<int> delayed_group;
std::vector<int> n_bank;
std::vector<int> n_bank_second;
std::vector<double> wgt_bank;
std::vector<int> n_delayed_bank; // MAX_DELAYED_GROUPS pitch
std::vector<int> surface;
std::vector<int> cell_born;
std::vector<int> material;
std::vector<int> material_last;
std::vector<BoundaryInfo> boundary;
std::vector<double> sqrtkT;
std::vector<double> sqrtkT_last;
std::vector<int> n_collision;
std::vector<char> write_track;
std::vector<uint64_t> seeds; // N_STREAMS pitch
std::vector<int> stream;
std::vector<ParticleBank> secondary_bank; // max_secondary_particles pitch
std::vector<int64_t> current_work;
std::vector<double> flux_derivs;         // tally_derivs.size() pitch
std::vector<FilterMatch> filter_matches; // tally_filters.size() pitch
std::vector<NuBank> nu_bank;             // max_fissions pitch
std::vector<double> keff_tally_absorption;
std::vector<double> keff_tally_collision;
std::vector<double> keff_tally_tracklength;
std::vector<double> keff_tally_leakage;
std::vector<char> trace;
std::vector<double> collision_distance;
std::vector<int> n_event;
#ifdef DAGMC
std::vector<moab::DagMC::RayHistory> history;
std::vector<Direction> last_dir;
#endif
std::vector<int64_t> n_progeny;
} // namespace soa

} // namespace openmc
