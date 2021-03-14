#include "openmc/cuda/calculate_xs.h"
#include "openmc/geometry.h" // find_cell
#include "openmc/search.h"

namespace openmc {
namespace gpu {

__constant__ unique_ptr<Material>* materials;
__constant__ unique_ptr<Nuclide>* nuclides;
__constant__ Particle* particles;
__constant__ NuclideMicroXS* micros;
__constant__ double energy_min_neutron;
__constant__ double energy_max_neutron;
__constant__ double log_spacing;
__constant__ unsigned number_nuclides;
__constant__ bool need_depletion_rx;

__managed__ unsigned managed_calculate_fuel_queue_index;
__managed__ unsigned managed_calculate_nonfuel_queue_index;

__global__ void process_calculate_xs_events_device(
  EventQueueItem* __restrict__ queue, unsigned queue_size)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  if (tid >= queue_size)
    return;
  Particle* __restrict__ p = particles + queue[tid].idx;
  auto const E = __ldg(&queue[tid].E);
  auto const mat_idx = __ldg(&queue[tid].material);

  // Store pre-collision particle properties
  p->wgt_last_ = p->wgt_;
  p->E_last_ = E;
  p->u_last_ = p->u();
  p->r_last_ = p->r();

  // Reset event variables
  p->event_ = TallyEvent::KILL;
  p->event_nuclide_ = NUCLIDE_NONE;
  p->event_mt_ = REACTION_NONE;

  p->macro_xs_.total = 0.0;
  p->macro_xs_.absorption = 0.0;
  p->macro_xs_.fission = 0.0;
  p->macro_xs_.nu_fission = 0.0;

  // Skip void material
  if (mat_idx == -1)
    return;

  Material const& m = *materials[mat_idx];

  unsigned i_log_union = std::log(E / energy_min_neutron) / log_spacing;

  // Add contribution from each nuclide in material
  auto const n_nuclides = m.nuclide_.size();
  for (int i = 0; i < n_nuclides; ++i) {
    auto const& i_nuclide = m.nuclide_[i];
    auto* __restrict__ micro {
      &micros[number_nuclides * queue[tid].idx + i_nuclide]};

    if (E != micro->last_E || p->sqrtkT_ != micro->last_sqrtkT) {
      auto const& nuclide = *nuclides[i_nuclide];
      micro->elastic = CACHE_INVALID;
      micro->thermal = 0.0;
      micro->thermal_elastic = 0.0;

      // Find the appropriate temperature index. why would someone use
      // nearest?
      double kT = p->sqrtkT_ * p->sqrtkT_;
      double f;

      int i_temp = -1;

      // Find temperatures that bound the actual temperature
      for (i_temp = 0; i_temp < nuclide.kTs_.size() - 1; ++i_temp) {
        if (nuclide.kTs_[i_temp] <= kT && kT < nuclide.kTs_[i_temp + 1])
          break;
      }

      // Randomly sample between temperature i and i+1
      f = (kT - nuclide.kTs_[i_temp]) /
          (nuclide.kTs_[i_temp + 1] - nuclide.kTs_[i_temp]);
      if (f > prn(p->seeds_))
        ++i_temp;

      const auto& grid {nuclide.grid_[i_temp]};
      int i_grid;
      if (E < grid.energy.front()) {
        i_grid = 0;
      } else if (E > grid.energy.back()) {
        i_grid = grid.energy.size() - 2;
      } else {
        // Determine bounding indices based on which equal log-spaced
        // interval the energy is in
        int i_low = __ldg(&grid.grid_index[i_log_union]);
        int i_high = __ldg(&grid.grid_index[i_log_union + 1]) + 1;

        // Perform binary search over reduced range
        i_grid = i_low + lower_bound_index_linear(
                           &grid.energy[i_low], &grid.energy[i_high], E);
      }
      const auto& xs_left {nuclide.xs_[i_temp][i_grid]};
      const auto& xs_right {nuclide.xs_[i_temp][i_grid + 1]};
      // check for rare case where two energy points are the same
      if (grid.energy[i_grid] == grid.energy[i_grid + 1])
        ++i_grid;

      // calculate interpolation factor
      f = (E - grid.energy[i_grid]) /
          (grid.energy[i_grid + 1] - grid.energy[i_grid]);

      micro->index_temp = i_temp;
      micro->index_grid = i_grid;
      micro->interp_factor = f;

      // Calculate all microscopic cross sections
      micro->total = (1.0 - f) * xs_left.total + f * xs_right.total;
      micro->absorption =
        (1.0 - f) * xs_left.absorption + f * xs_right.absorption;

      if (nuclide.fissionable_) {
        // Calculate microscopic nuclide total cross section
        micro->fission = (1.0 - f) * xs_left.fission + f * xs_right.fission;

        // Calculate microscopic nuclide nu-fission cross section
        micro->nu_fission =
          (1.0 - f) * xs_left.nu_fission + f * xs_right.nu_fission;
      } else {
        micro->fission = 0.0;
        micro->nu_fission = 0.0;
      }

      // Calculate microscopic nuclide photon production cross section
      micro->photon_prod =
        (1.0 - f) * xs_left.photon_production + f * xs_right.photon_production;

      micro->index_sab = C_NONE;
      micro->sab_frac = 0.0;
      micro->last_E = E;
      micro->last_sqrtkT = p->sqrtkT_;
    }

    double const& atom_density = m.atom_density_[i];
    p->macro_xs_.total += atom_density * micro->total;
    p->macro_xs_.absorption += atom_density * micro->absorption;
    p->macro_xs_.fission += atom_density * micro->fission;
    p->macro_xs_.nu_fission += atom_density * micro->nu_fission;
  }
}

} // namespace gpu
} // namespace openmc
