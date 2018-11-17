#ifndef OPENMC_PHYSICS_H
#define OPENMC_PHYSICS_H

#include "openmc/bank.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/reaction.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

//! Sample a nuclide and reaction and then calls the appropriate routine
extern "C" void collision(Particle* p);

//! Samples an incident neutron reaction
extern "C" void sample_neutron_reaction(Particle* p);

//! Samples an element based on the macroscopic cross sections for each nuclide
//! within a material and then samples a reaction for that element and calls the
//! appropriate routine to process the physics.
extern "C" void sample_photon_reaction(Particle* p);

//! Terminates the particle and either deposits all energy locally
//! (electron_treatment = ELECTRON_LED) or creates secondary bremsstrahlung
//! photons from electron deflections with charged particles (electron_treatment
//! = ELECTRON_TTB).
extern "C" void sample_electron_reaction(Particle* p);

//! Terminates the particle and either deposits all energy locally
//! (electron_treatment = ELECTRON_LED) or creates secondary bremsstrahlung
//! photons from electron deflections with charged particles (electron_treatment
//! = ELECTRON_TTB). Two annihilation photons of energy MASS_ELECTRON_EV (0.511
//! MeV) are created and travel in opposite directions.
extern "C" void sample_positron_reaction(Particle* p);

extern "C" void sample_nuclide(const Particle* p, int mt, int* i_nuclide, int* i_nuc_mat);

//! Determine the average total, prompt, and delayed neutrons produced from
//! fission and creates appropriate bank sites.
extern "C" void create_fission_sites(Particle* p, int i_nuclide, int i_rx,
  Bank* bank_array, int64_t* bank_size, int64_t bank_capacity);

// void sample_element(Particle* p);

extern "C" int sample_fission(int i_nuclide, double E);

// void sample_photon_product(int i_nuclide, double E, int* i_rx, int* i_product);

extern "C" void absorption(Particle* p, int i_nuclide);

extern "C" void scatter(Particle*, int i_nuclide, int i_nuc_mat);

// void elastic_scatter(int i_nuclide, const Reaction& rx, double kT, double* E,
//   Direction* u, double* mu_lab, double* wgt);

// void sab_scatter(int i_nuclide, int i_sab, double* E, Direction* u, double* mu);

// void sample_target_velocity(int i_nuclide, Direction* v_target, double E, Direction u,
//   Direction v_neut, double* wgt, double xs_eff, double kT);

// void sample_cxs_target_velocity(int i_nuclide, Direction* v_target, double E, Direction u,
//   double kT);

extern "C" void sample_fission_neutron(int i_nuclide, int i_rx, double E_in, Bank* site);

// void inelastic_scatter(int i_nuclide, const Reaction& rx, Particle* p);

extern "C" void sample_secondary_photons(Particle* p, int i_nuclide);

} // namespace openmc

#endif // OPENMC_PHYSICS_H
