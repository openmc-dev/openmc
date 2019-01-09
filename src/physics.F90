module physics

  use constants
  use error,                  only: fatal_error, warning, write_message
  use material_header,        only: Material, materials
  use math
  use message_passing
  use nuclide_header
  use particle_header
  use photon_header
  use photon_physics,         only: rayleigh_scatter, compton_scatter, &
                                    atomic_relaxation, pair_production, &
                                    thick_target_bremsstrahlung
  use physics_common
  use random_lcg,             only: prn
  use settings
  use simulation_header
  use tally_header

  implicit none

  interface
    subroutine collision(p) bind(C)
      import Particle
      type(Particle), intent(inout) :: p
    end subroutine
  end interface

contains

!===============================================================================
! SAMPLE_PHOTON_REACTION samples an element based on the macroscopic cross
! sections for each nuclide within a material and then samples a reaction for
! that element and calls the appropriate routine to process the physics.
!===============================================================================

  subroutine sample_photon_reaction(p) bind(C)
    type(Particle), intent(inout) :: p

    integer :: i_shell      ! index in subshells
    integer :: i_grid       ! index on energy grid
    integer :: i_element    ! index in nuclides array
    integer :: i_start      ! threshold index
    real(8) :: prob         ! cumulative probability
    real(8) :: cutoff       ! sampled total cross section
    real(8) :: f            ! interpolation factor
    real(8) :: xs           ! photoionization cross section
    real(8) :: r            ! random number
    real(8) :: prob_after
    real(8) :: alpha        ! photon energy divided by electron rest mass
    real(8) :: alpha_out    ! outgoing photon energy over electron rest mass
    real(8) :: mu           ! scattering cosine
    real(8) :: mu_electron  ! electron scattering cosine
    real(8) :: mu_positron  ! positron scattering cosine
    real(8) :: phi          ! azimuthal angle
    real(8) :: uvw(3)       ! new direction
    real(8) :: rel_vel      ! relative velocity of electron
    real(8) :: e_b          ! binding energy of electron
    real(8) :: E_electron   ! electron energy
    real(8) :: E_positron   ! positron energy

    ! Kill photon if below energy cutoff -- an extra check is made here because
    ! photons with energy below the cutoff may have been produced by neutrons
    ! reactions or atomic relaxation
    if (p % E < energy_cutoff(PHOTON)) then
      p % E = ZERO
      p % alive = .false.
      return
    end if

    ! Sample element within material
    i_element = sample_element(p)
    p % event_nuclide = i_element

    ! Calculate photon energy over electron rest mass equivalent
    alpha = p % E/MASS_ELECTRON_EV

    ! For tallying purposes, this routine might be called directly. In that
    ! case, we need to sample a reaction via the cutoff variable
    prob = ZERO
    cutoff = prn() * micro_photon_xs(i_element) % total

    associate (elm => elements(i_element))
      ! Coherent (Rayleigh) scattering
      prob = prob + micro_photon_xs(i_element) % coherent
      if (prob > cutoff) then
        call rayleigh_scatter(elm, alpha, mu)
        p % coord(1) % uvw = rotate_angle(p % coord(1) % uvw, mu)
        p % event_MT = COHERENT
        return
      end if

      ! Incoherent (Compton) scattering
      prob = prob + micro_photon_xs(i_element) % incoherent
      if (prob > cutoff) then
        call compton_scatter(elm, alpha, alpha_out, mu, i_shell, .true.)

        ! Determine binding energy of shell. The binding energy is zero if
        ! doppler broadening is not used.
        if (i_shell == 0) then
          e_b = ZERO
        else
          e_b = elm % binding_energy(i_shell)
        end if

        ! Create Compton electron
        E_electron = (alpha - alpha_out)*MASS_ELECTRON_EV - e_b
        mu_electron = (alpha - alpha_out*mu) &
             / sqrt(alpha**2 + alpha_out**2 - TWO*alpha*alpha_out*mu)
        phi = TWO*PI*prn()
        uvw = rotate_angle(p % coord(1) % uvw, mu_electron, phi)
        call particle_create_secondary(p, uvw, E_electron, ELECTRON, .true._C_BOOL)

        ! TODO: Compton subshell data does not match atomic relaxation data
        ! Allow electrons to fill orbital and produce auger electrons
        ! and fluorescent photons
        if (i_shell > 0) then
          call atomic_relaxation(p, elm, i_shell)
        end if

        phi = phi + PI
        p % E = alpha_out*MASS_ELECTRON_EV
        p % coord(1) % uvw = rotate_angle(p % coord(1) % uvw, mu, phi)
        p % event_MT = INCOHERENT
        return
      end if

      ! Photoelectric effect
      prob_after = prob + micro_photon_xs(i_element) % photoelectric
      if (prob_after > cutoff) then
        do i_shell = 1, size(elm % shells)
          ! Get grid index and interpolation factor
          i_grid = micro_photon_xs(i_element) % index_grid
          f      = micro_photon_xs(i_element) % interp_factor

          ! Check threshold of reaction
          i_start = elm % shells(i_shell) % threshold
          if (i_grid <= i_start) cycle

          ! Evaluation subshell photoionization cross section
          xs = exp(elm % shells(i_shell) % cross_section(i_grid - i_start) + &
               f*(elm % shells(i_shell) % cross_section(i_grid + 1 - i_start) - &
               elm % shells(i_shell) % cross_section(i_grid - i_start)))

          prob = prob + xs
          if (prob > cutoff) then
            E_electron = p % E - elm % shells(i_shell) % binding_energy

            ! Sample mu using non-relativistic Sauter distribution.
            ! See Eqns 3.19 and 3.20 in "Implementing a photon physics
            ! model in Serpent 2" by Toni Kaltiaisenaho
            SAMPLE_MU: do
              r = prn()
              if (FOUR * (ONE - r) * r >= prn()) then
                rel_vel = sqrt(E_electron * (E_electron + TWO * MASS_ELECTRON_EV))&
                     / (E_electron + MASS_ELECTRON_EV)
                mu = (TWO * r + rel_vel - ONE) / &
                     (TWO * rel_vel * r - rel_vel + ONE)
                exit SAMPLE_MU
              end if
            end do SAMPLE_MU

            phi = TWO*PI*prn()
            uvw(1) = mu
            uvw(2) = sqrt(ONE - mu*mu)*cos(phi)
            uvw(3) = sqrt(ONE - mu*mu)*sin(phi)

            ! Create secondary electron
            call particle_create_secondary(p, uvw, E_electron, ELECTRON, &
                 run_CE=.true._C_BOOL)

            ! Allow electrons to fill orbital and produce auger electrons
            ! and fluorescent photons
            call atomic_relaxation(p, elm, i_shell)
            p % event_MT = 533 + elm % shells(i_shell) % index_subshell
            p % alive = .false.
            p % E = ZERO

            return
          end if
        end do
      end if
      prob = prob_after

      ! Pair production
      prob = prob + micro_photon_xs(i_element) % pair_production
      if (prob > cutoff) then
        call pair_production(elm, alpha, E_electron, E_positron, mu_electron, &
             mu_positron)

        ! Create secondary electron
        uvw = rotate_angle(p % coord(1) % uvw, mu_electron)
        call particle_create_secondary(p, uvw, E_electron, ELECTRON, .true._C_BOOL)

        ! Create secondary positron
        uvw = rotate_angle(p % coord(1) % uvw, mu_positron)
        call particle_create_secondary(p, uvw, E_positron, POSITRON, .true._C_BOOL)

        p % event_MT = PAIR_PROD
        p % alive = .false.
        p % E = ZERO
      end if

    end associate

  end subroutine sample_photon_reaction

!===============================================================================
! SAMPLE_ELEMENT
!===============================================================================

  function sample_element(p) result(i_element)
    type(Particle), intent(in) :: p
    integer                    :: i_element

    integer :: i
    real(8) :: prob
    real(8) :: cutoff
    real(8) :: atom_density ! atom density of nuclide in atom/b-cm
    real(8) :: sigma        ! microscopic total xs for nuclide

    associate (mat => materials(p % material))
      ! Sample cumulative distribution function
      cutoff = prn() * material_xs % total

      i = 0
      prob = ZERO
      do while (prob < cutoff)
        i = i + 1

        ! Check to make sure that a nuclide was sampled
        if (i > mat % n_nuclides) then
          call particle_write_restart(p)
          call fatal_error("Did not sample any element during collision.")
        end if

        ! Find atom density
        i_element    = mat % element(i)
        atom_density = mat % atom_density(i)

        ! Determine microscopic cross section
        sigma = atom_density * micro_photon_xs(i_element) % total

        ! Increment probability to compare to cutoff
        prob = prob + sigma
      end do
    end associate

  end function sample_element

end module physics
