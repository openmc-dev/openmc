#include "openmc/random_ray/iteration.h"

namespace openmc {
  
  int openmc_run_random_ray(void);

  void update_neutron_source(double k_eff)
  {
    double inverse_k_eff = 1.0 / k_eff;
    int negroups = data::mg.num_energy_groups_;

    #pragma omp parallel for
    for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
      int material = random_ray::materials[sr]; 
      for (int energy_group_out = 0; energy_group_out < negroups; energy_group_out++) {
        float Sigma_t = data::mg.macro_xs_[material].get_xs(MgxsType::TOTAL, energy_group_out, NULL, NULL, NULL);
        float scatter_source = 0.0f;
        float fission_source = 0.0f;
        for (int energy_group_in = 0; energy_group_in < negroups; energy_group_in++) {
          float scalar_flux = random_ray::scalar_flux_old[sr * negroups + energy_group_in];
          float Sigma_s = data::mg.macro_xs_[material].get_xs(MgxsType::NU_SCATTER, energy_group_in, &energy_group_out, NULL, NULL);
          float nu_Sigma_f = data::mg.macro_xs_[material].get_xs(MgxsType::NU_FISSION, energy_group_in, NULL, NULL, NULL);
          float Chi = data::mg.macro_xs_[material].get_xs(MgxsType::CHI_PROMPT, energy_group_in, &energy_group_out, NULL, NULL);
          scatter_source += Sigma_s    * scalar_flux;
          fission_source += nu_Sigma_f * scalar_flux * Chi;
        }
        fission_source *= inverse_k_eff;
        float new_isotropic_source = (scatter_source + fission_source)  / Sigma_t;
        random_ray::source[sr * negroups + energy_group_out] = new_isotropic_source;
      }
    }
  }
  
  void normalize_scalar_flux_and_volumes(double total_active_distance_per_iteration, int iter);

  void add_source_to_scalar_flux(double total_active_distance_per_iteration, int iter);

  double compute_k_eff(double k_eff_old);
  
  void tally_fission_rates(void);
  
  double calculate_miss_rate(void);


} // namespace openmc
