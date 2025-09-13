#ifndef OPENMC_TALLIES_TALLY_SCORING_H
#define OPENMC_TALLIES_TALLY_SCORING_H

#include "openmc/particle.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"

namespace openmc {

//==============================================================================
//! An iterator over all combinations of a tally's matching filter bins.
//
//! This iterator handles two distinct tasks.  First, it maps the N-dimensional
//! space created by the indices of N filters onto a 1D sequence.  In other
//! words, it provides a single number that uniquely identifies a combination of
//! bins for many filters.  Second, it handles the task of finding each valid
//! combination of filter bins given that each filter can have 1 or 2 or many
//! bins that are valid for the current tally event.
//==============================================================================

class FilterBinIter {
public:
  //! Construct an iterator over bins that match a given particle's state.
  FilterBinIter(const Tally& tally, Particle& p);

  //! Construct an iterator over all filter bin combinations.
  //
  //! \param end if true, the returned iterator indicates the end of a loop.
  FilterBinIter(
    const Tally& tally, bool end, vector<FilterMatch>* particle_filter_matches);

  bool operator==(const FilterBinIter& other) const
  {
    return index_ == other.index_;
  }

  bool operator!=(const FilterBinIter& other) const
  {
    return !(*this == other);
  }

  FilterBinIter& operator++();

  int index_ {1};
  double weight_ {1.};

  vector<FilterMatch>& filter_matches_;

private:
  void compute_index_weight();

  const Tally& tally_;
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Score tallies using a 1 / Sigma_t estimate of the flux.
//
//! This is triggered after every collision.  It is invalid for tallies that
//! require post-collison information because it can score reactions that didn't
//! actually occur, and we don't a priori know what the outcome will be for
//! reactions that we didn't sample.  It is assumed the material is not void
//! since collisions do not occur in voids.
//
//! \param p The particle being tracked
void score_collision_tally(Particle& p);

//! Score tallies using the point estimator.
//! This is triggered after every collision.
//! \param p The particle being tracked
void score_point_tally(Particle& p);
//! Retrieve the position and exclusion sphere (R0) of a detector.
//!
//! This function fetches the spatial position [x, y, z] and the exclusion
//! sphere radius (R0) of the point detector for a specified tally index. It
//! ensures that the tally contains exactly three spatial coordinates and an R0
//! value.
//!
//! \param det_pos Array to store the detector position [x, y, z] and exclusion
//! sphere radius R0. \param i_tally Index of the tally for which the detector
//! position and R0 are required.
void get_det_pos(double (&det_pos)[4], int i_tally);
//! Score to point tally using data from a source site.
//!
//! This function calculates and scores flux for active point tallies
//! based on the properties of a given source site. It accounts for the source's
//! spatial position, directional distribution, and type to compute the
//! appropriate contribution to the tally. The function handles various source
//! angle distributions, including isotropic and polar-azimuthal types.
//!
//! \param src Pointer to the source site containing particle and source
//! information. \throws RuntimeError If an unexpected angle distribution type
//! is encountered.
void score_point_tally_from_source(const SourceSite* src);
//! Calculate the total mean free path (MFP) for a ghost particle along a
//! specified distance.
//!
//! This function computes the total mean free path (MFP) for a ghost particle
//! traveling from a collision point to a detector located at a given distance.
//! The direction of travel is towards the point detector. The cross-section
//! used in the calculation is computed based on the ghost particle's energy,
//! which remains constant throughout the traversal. The function accounts for
//! shielding by iteratively advancing the particle to boundaries and summing
//! contributions to the MFP along the path.
//!
//! \param ghost_particle The particle whose MFP is being calculated, with its
//! direction towards the point detector. \param total_distance The distance
//! [cm] from the collision point to the detector. \return The total mean free
//! path for the particle over the specified distance. \throws RuntimeError If
//! the particle encounters unexpected behavior during boundary crossing.
double get_MFP(Particle& ghost_particle, double total_distance);
//! Calculate the probability density function (PDF) for elastic and inelastic
//! scattering to a point detector.
//!
//! This function computes the PDF for a particle scattering (elastic or
//! inelastic) to a detector at a specified position. All computations are
//! performed fully relativistically, including energy, momentum, and
//! transformations between the center-of-mass (CM) frame and the lab frame. The
//! Jacobian used in the calculations corresponds to a Lorentz transformation,
//! ensuring accurate mapping between frames of reference.
//!
//! Ghost particles are created to represent possible scattering trajectories.
//! Each ghost particle is initialized with the correct energy and direction as
//! if the particle had scattered towards the point detector.
//!
//! In the case of inelastic scattering, the post-collision energy of the
//! outgoing neutron is set in the CM frame, leading to an effective mass
//! (\(m_4\)) for the residual nucleus. For inelastic collisions,
//! \(m_4 \neq m_2\), as \(m_4\) reflects the excitation state of the residual
//! nucleus.

//! \param det_pos The position of the point detector in [x, y, z, R0].
//! \param p The incident particle being analyzed for scattering.
//! \param mu_cm A vector to store the cosine of the scattering angles in the CM
//! frame. \param Js A vector to store the Jacobians for transforming the
//! distributions to the lab frame using Lorentz transformations. \param
//! ghost_particles A vector to store the ghost particles representing
//! scattering trajectories.
//!                        Each ghost particle is initialized with the correct
//!                        energy and direction towards the point detector.
//! \param E3k_cm_given The center-of-mass kinetic energy of one outgoing
//! particle, if specified.
//!                      If set to -1, it is calculated based on kinematics.
void get_pdf_to_point_elastic(double det_pos[4], Particle& p,
  std::vector<double>& mu_cm, std::vector<double>& Js,
  std::vector<Particle>& ghost_particles, double E3k_cm_given = -1);
//! Score a ghost particle contribution to a point tally.
//!
//! This function calculates and scores the contribution of a ghost particle
//! to a specified point tally. The scoring process includes computing the flux
//! at the detector location.
//! It accounts for cases where the particle path is fully or partially within
//! the exclusion sphere (\(R_0\)) or outside it. The scoring respects particle
//! type, ensuring only neutrons are considered.
//!
//! \param ghost_p The ghost particle being scored.
//! \param pdf_lab The probability density function value in the lab frame for
//! the particle direction. \param i_tally The index of the tally to which the
//! ghost particle contribution is scored.
//!
//! \note Only continuous energy mode is currently supported; multi-group mode
//! is not implemented.
void score_ghost_particle(Particle& ghost_p, double pdf_lab, int i_tally);
//! Perform a Lorentz boost of a four-vector.
//!
//! This function boosts a four-vector \( B \) (in the lab frame) into the rest
//! frame defined by another four-vector \( A \). The result of the boost is
//! stored in \( X \), representing \( B \) in the rest frame of \( A \).
//!
//! \param A The four-vector defining the rest frame [E, px, py, pz].
//! \param B The four-vector to be boosted [E, px, py, pz].
//! \param X The resulting four-vector after the Lorentz boost [E, px, py, pz]
//! (rest frame).
void boostf(double A[4], double B[4], double X[4]);

//! Analog tallies are triggered at every collision, not every event.
//
//! \param p The particle being tracked
void score_analog_tally_ce(Particle& p);

//! Score tallies based on a simple count of events (for multigroup).
//
void score_general_ce_nonanalog(Particle& p, int i_tally, int start_index,
  int filter_index, double filter_weight, int i_nuclide, double atom_density,
  double flux);

//! Analog tallies are triggered at every collision, not every event.
//
//! \param p The particle being tracked
void score_analog_tally_mg(Particle& p);

//! Score tallies using a tracklength estimate of the flux.
//
//! This is triggered at every event (surface crossing, lattice crossing, or
//! collision) and thus cannot be done for tallies that require post-collision
//! information.
//
//! \param p The particle being tracked
//! \param distance The distance in [cm] traveled by the particle
void score_tracklength_tally(Particle& p, double distance);

//! Score time filtered tallies using a tracklength estimate of the flux.
//
//! This is triggered at every event (surface crossing, lattice crossing, or
//! collision) and thus cannot be done for tallies that require post-collision
//! information.
//
//! \param p The particle being tracked
//! \param total_distance The distance in [cm] traveled by the particle
void score_timed_tracklength_tally(Particle& p, double total_distance);

//! Score surface or mesh-surface tallies for particle currents.
//
//! \param p The particle being tracked
//! \param tallies A vector of the indices of the tallies to score to
void score_surface_tally(Particle& p, const vector<int>& tallies);

//! Score the pulse-height tally
//! This is triggered at the end of every particle history
//
//! \param p The particle being tracked
//! \param tallies A vector of the indices of the tallies to score to
void score_pulse_height_tally(Particle& p, const vector<int>& tallies);

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_SCORING_H
