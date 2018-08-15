module physics

  use algorithm,              only: binary_search
  use constants
  use endf,                   only: reaction_name
  use error,                  only: fatal_error, warning, write_message
  use material_header,        only: Material, materials
  use math
  use mesh_header,            only: meshes
  use message_passing
  use nuclide_header
  use particle_header
  use photon_header
  use photon_physics,         only: rayleigh_scatter, compton_scatter, &
                                    atomic_relaxation, pair_production, &
                                    thick_target_bremsstrahlung
  use physics_common
  use random_lcg,             only: prn, advance_prn_seed, prn_set_stream
  use reaction_header,        only: Reaction
  use sab_header,             only: sab_tables
  use settings
  use simulation_header
  use string,                 only: to_str
  use tally_header

  implicit none

contains

!===============================================================================
! COLLISION samples a nuclide and reaction and then calls the appropriate
! routine for that reaction
!===============================================================================

  subroutine collision(p)

    type(Particle), intent(inout) :: p

    ! Add to collision counter for particle
    p % n_collision = p % n_collision + 1

    ! Sample reaction for the material the particle is in
    if (p % type == NEUTRON) then
      call sample_neutron_reaction(p)
    else if (p % type == PHOTON) then
      call sample_photon_reaction(p)
    else if (p % type == ELECTRON) then
      call sample_electron_reaction(p)
    else if (p % type == POSITRON) then
      call sample_positron_reaction(p)
    end if

    ! Kill particle if energy falls below cutoff
    if (p % E < energy_cutoff(p % type)) then
      p % alive = .false.
      p % wgt = ZERO
      p % last_wgt = ZERO
    end if

    ! Display information about collision
    if (verbosity >= 10 .or. trace) then
      if (p % type == NEUTRON) then
        call write_message("    " // trim(reaction_name(p % event_MT)) &
             &// " with " // trim(adjustl(nuclides(p % event_nuclide) % name)) &
             &// ". Energy = " // trim(to_str(p % E)) // " eV.")
      else
        call write_message("    " // trim(reaction_name(p % event_MT)) &
             &// " with " // trim(adjustl(elements(p % event_nuclide) % name)) &
             &// ". Energy = " // trim(to_str(p % E)) // " eV.")
      end if
    end if

  end subroutine collision

!===============================================================================
! SAMPLE_NEUTRON_REACTION samples a nuclide based on the macroscopic cross
! sections for each nuclide within a material and then samples a reaction for
! that nuclide and calls the appropriate routine to process the physics. Note
! that there is special logic when suvival biasing is turned on since fission
! and disappearance are treated implicitly.
!===============================================================================

  subroutine sample_neutron_reaction(p)

    type(Particle), intent(inout) :: p

    integer :: i_nuclide    ! index in nuclides array
    integer :: i_nuc_mat    ! index in material's nuclides array
    integer :: i_reaction   ! index in nuc % reactions array
    type(Nuclide), pointer :: nuc

    call sample_nuclide(p, 'total  ', i_nuclide, i_nuc_mat)

    ! Get pointer to table
    nuc => nuclides(i_nuclide)

    ! Save which nuclide particle had collision with
    p % event_nuclide = i_nuclide

    ! Create fission bank sites. Note that while a fission reaction is sampled,
    ! it never actually "happens", i.e. the weight of the particle does not
    ! change when sampling fission sites. The following block handles all
    ! absorption (including fission)

    if (nuc % fissionable) then
      if (run_mode == MODE_EIGENVALUE) then
        call sample_fission(i_nuclide, p % E, i_reaction)
        call create_fission_sites(p, i_nuclide, i_reaction, fission_bank, n_bank)
      elseif (run_mode == MODE_FIXEDSOURCE .and. create_fission_neutrons) then
        call sample_fission(i_nuclide, p % E, i_reaction)
        call create_fission_sites(p, i_nuclide, i_reaction, &
             p % secondary_bank, p % n_secondary)
      end if
    end if

    ! Create secondary photons
    if (photon_transport) then
      call prn_set_stream(STREAM_PHOTON)
      call sample_secondary_photons(p, i_nuclide)
      call prn_set_stream(STREAM_TRACKING)
    end if

    ! If survival biasing is being used, the following subroutine adjusts the
    ! weight of the particle. Otherwise, it checks to see if absorption occurs

    if (micro_xs(i_nuclide) % absorption > ZERO) then
      call absorption(p, i_nuclide)
    else
      p % absorb_wgt = ZERO
    end if
    if (.not. p % alive) return

    ! Sample a scattering reaction and determine the secondary energy of the
    ! exiting neutron
    call scatter(p, i_nuclide, i_nuc_mat)

    ! Advance URR seed stream 'N' times after energy changes
    if (p % E /= p % last_E) then
      call prn_set_stream(STREAM_URR_PTABLE)
      call advance_prn_seed(size(nuclides, kind=8))
      call prn_set_stream(STREAM_TRACKING)
    end if

    ! Play russian roulette if survival biasing is turned on
    if (survival_biasing) then
      call russian_roulette(p)
      if (.not. p % alive) return
    end if

  end subroutine sample_neutron_reaction

!===============================================================================
! SAMPLE_PHOTON_REACTION samples an element based on the macroscopic cross
! sections for each nuclide within a material and then samples a reaction for
! that element and calls the appropriate routine to process the physics.
!===============================================================================

  subroutine sample_photon_reaction(p)
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
! SAMPLE_ELECTRON_REACTION terminates the particle and either deposits all
! energy locally (electron_treatment = ELECTRON_LED) or creates secondary
! bremsstrahlung photons from electron deflections with charged particles
! (electron_treatment = ELECTRON_TTB).
!===============================================================================

  subroutine sample_electron_reaction(p)
    type(Particle), intent(inout) :: p

    real(8) :: E_lost  ! energy lost to bremsstrahlung photons

    ! TODO: create reaction types

    if (electron_treatment == ELECTRON_TTB) then
      call thick_target_bremsstrahlung(p, E_lost)
    end if

    p % E = ZERO
    p % alive = .false.

  end subroutine sample_electron_reaction

!===============================================================================
! SAMPLE_POSITRON_REACTION terminates the particle and either deposits all
! energy locally (electron_treatment = ELECTRON_LED) or creates secondary
! bremsstrahlung photons from electron deflections with charged particles
! (electron_treatment = ELECTRON_TTB). Two annihilation photons of energy
! MASS_ELECTRON_EV (0.511 MeV) are created and travel in opposite directions.
!===============================================================================

  subroutine sample_positron_reaction(p)
    type(Particle), intent(inout) :: p

    real(8) :: mu     ! scattering cosine
    real(8) :: phi    ! azimuthal angle
    real(8) :: uvw(3) ! new direction

    real(8) :: E_lost  ! energy lost to bremsstrahlung photons

    ! TODO: create reaction types

    if (electron_treatment == ELECTRON_TTB) then
      call thick_target_bremsstrahlung(p, E_lost)
    end if

    ! Sample angle isotropically
    mu = TWO*prn() - ONE
    phi = TWO*PI*prn()
    uvw(1) = mu
    uvw(2) = sqrt(ONE - mu*mu)*cos(phi)
    uvw(3) = sqrt(ONE - mu*mu)*sin(phi)

    ! Create annihilation photon pair traveling in opposite directions
    call particle_create_secondary(p, uvw, MASS_ELECTRON_EV, PHOTON, .true._C_BOOL)
    call particle_create_secondary(p, -uvw, MASS_ELECTRON_EV, PHOTON, .true._C_BOOL)

    p % E = ZERO
    p % alive = .false.

  end subroutine sample_positron_reaction

!===============================================================================
! SAMPLE_NUCLIDE
!===============================================================================

  subroutine sample_nuclide(p, base, i_nuclide, i_nuc_mat)

    type(Particle), intent(in) :: p
    character(7),   intent(in) :: base      ! which reaction to sample based on
    integer, intent(out)       :: i_nuclide
    integer, intent(out)       :: i_nuc_mat

    real(8) :: prob
    real(8) :: cutoff
    real(8) :: atom_density ! atom density of nuclide in atom/b-cm
    real(8) :: sigma        ! microscopic total xs for nuclide
    type(Material), pointer :: mat

    ! Get pointer to current material
    mat => materials(p % material)

    ! Sample cumulative distribution function
    select case (base)
    case ('total')
      cutoff = prn() * material_xs % total
    case ('scatter')
      cutoff = prn() * (material_xs % total - material_xs % absorption)
    case ('fission')
      cutoff = prn() * material_xs % fission
    end select

    i_nuc_mat = 0
    prob = ZERO
    do while (prob < cutoff)
      i_nuc_mat = i_nuc_mat + 1

      ! Check to make sure that a nuclide was sampled
      if (i_nuc_mat > mat % n_nuclides) then
        call particle_write_restart(p)
        call fatal_error("Did not sample any nuclide during collision.")
      end if

      ! Find atom density
      i_nuclide    = mat % nuclide(i_nuc_mat)
      atom_density = mat % atom_density(i_nuc_mat)

      ! Determine microscopic cross section
      select case (base)
      case ('total')
        sigma = atom_density * micro_xs(i_nuclide) % total
      case ('scatter')
        sigma = atom_density * (micro_xs(i_nuclide) % total - &
             micro_xs(i_nuclide) % absorption)
      case ('fission')
        sigma = atom_density * micro_xs(i_nuclide) % fission
      end select

      ! Increment probability to compare to cutoff
      prob = prob + sigma
    end do

  end subroutine sample_nuclide

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

!===============================================================================
! SAMPLE_FISSION
!===============================================================================

  subroutine sample_fission(i_nuclide, E, i_reaction)
    integer, intent(in)  :: i_nuclide  ! index in nuclides array
    real(8), intent(in)  :: E          ! incident neutron energy
    integer, intent(out) :: i_reaction ! index in nuc % reactions array

    integer :: i
    integer :: i_grid
    integer :: i_temp
    integer :: threshold
    real(8) :: f
    real(8) :: prob
    real(8) :: cutoff
    type(Nuclide), pointer :: nuc

    ! Get pointer to nuclide
    nuc => nuclides(i_nuclide)

    ! If we're in the URR, by default use the first fission reaction. We also
    ! default to the first reaction if we know that there are no partial fission
    ! reactions
    if (micro_xs(i_nuclide) % use_ptable .or. &
         .not. nuc % has_partial_fission) then
      i_reaction = nuc % index_fission(1)
      return
    end if

    ! Check to see if we are in a windowed multipole range.  WMP only supports
    ! the first fission reaction.
    if (nuc % mp_present) then
      if (E >= nuc % multipole % start_E .and. &
           E <= nuc % multipole % end_E) then
        i_reaction = nuc % index_fission(1)
        return
      end if
    end if

    ! Get grid index and interpolatoin factor and sample fission cdf
    i_temp = micro_xs(i_nuclide) % index_temp
    i_grid = micro_xs(i_nuclide) % index_grid
    f      = micro_xs(i_nuclide) % interp_factor
    cutoff = prn() * micro_xs(i_nuclide) % fission
    prob   = ZERO

    ! Loop through each partial fission reaction type

    FISSION_REACTION_LOOP: do i = 1, nuc % n_fission
      i_reaction = nuc % index_fission(i)

      associate (rx => nuc % reactions(i_reaction))
        ! if energy is below threshold for this reaction, skip it
        threshold = rx % xs_threshold(i_temp)
        if (i_grid < threshold) cycle

        ! add to cumulative probability
        prob = prob + ((ONE - f) * rx % xs(i_temp, i_grid - threshold + 1) &
             + f*(rx % xs(i_temp, i_grid - threshold + 2)))
      end associate

      ! Create fission bank sites if fission occurs
      if (prob > cutoff) exit FISSION_REACTION_LOOP
    end do FISSION_REACTION_LOOP

  end subroutine sample_fission

!===============================================================================
! SAMPLE_PHOTON_PRODUCT
!===============================================================================

  subroutine sample_photon_product(i_nuclide, E, i_reaction, i_product)
    integer, intent(in)  :: i_nuclide  ! index in nuclides array
    real(8), intent(in)  :: E          ! energy of neutron
    integer, intent(out) :: i_reaction ! index in nuc % reactions array
    integer, intent(out) :: i_product  ! index in reaction % products array

    integer :: i_grid
    integer :: i_temp
    integer :: threshold
    integer :: last_valid_reaction
    integer :: last_valid_product
    real(8) :: f
    real(8) :: prob
    real(8) :: cutoff
    real(8) :: yield

    ! Get pointer to nuclide
    associate (nuc => nuclides(i_nuclide))

      ! Get grid index and interpolation factor and sample photon production cdf
      i_temp = micro_xs(i_nuclide) % index_temp
      i_grid = micro_xs(i_nuclide) % index_grid
      f      = micro_xs(i_nuclide) % interp_factor
      cutoff = prn() * micro_xs(i_nuclide) % photon_prod
      prob   = ZERO

      ! Loop through each reaction type
      REACTION_LOOP: do i_reaction = 1, size(nuc % reactions)
        associate (rx => nuc % reactions(i_reaction))
          threshold = rx % xs_threshold(i_temp)

          ! if energy is below threshold for this reaction, skip it
          if (i_grid < threshold) cycle

          do i_product = 1, rx % products_size()
            if (rx % product_particle(i_product) == PHOTON) then
              ! add to cumulative probability
              yield = rx % product_yield(i_product, E)
              prob = prob + ((ONE - f) * rx % xs(i_temp, i_grid - threshold + 1) &
                   + f*(rx % xs(i_temp, i_grid - threshold + 2))) * yield

              if (prob > cutoff) return
              last_valid_reaction = i_reaction
              last_valid_product = i_product
            end if
          end do
        end associate
      end do REACTION_LOOP
    end associate

    i_reaction = last_valid_reaction
    i_product = last_valid_product

  end subroutine sample_photon_product

!===============================================================================
! ABSORPTION
!===============================================================================

  subroutine absorption(p, i_nuclide)
    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: i_nuclide

    if (survival_biasing) then
      ! Determine weight absorbed in survival biasing
      p % absorb_wgt = p % wgt * micro_xs(i_nuclide) % absorption / &
           micro_xs(i_nuclide) % total

      ! Adjust weight of particle by probability of absorption
      p % wgt = p % wgt - p % absorb_wgt
      p % last_wgt = p % wgt

      ! Score implicit absorption estimate of keff
      if (run_mode == MODE_EIGENVALUE) then
        global_tally_absorption = global_tally_absorption + p % absorb_wgt * &
             micro_xs(i_nuclide) % nu_fission / micro_xs(i_nuclide) % absorption
      end if
    else
      ! See if disappearance reaction happens
      if (micro_xs(i_nuclide) % absorption > &
           prn() * micro_xs(i_nuclide) % total) then
        ! Score absorption estimate of keff
        if (run_mode == MODE_EIGENVALUE) then
          global_tally_absorption = global_tally_absorption + p % wgt * &
               micro_xs(i_nuclide) % nu_fission / micro_xs(i_nuclide) % absorption
        end if

        p % alive = .false.
        p % event = EVENT_ABSORB
        p % event_MT = N_DISAPPEAR
      end if
    end if

  end subroutine absorption

!===============================================================================
! SCATTER
!===============================================================================

  subroutine scatter(p, i_nuclide, i_nuc_mat)
    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: i_nuclide
    integer,        intent(in)    :: i_nuc_mat

    integer :: i
    integer :: j
    integer :: i_temp
    integer :: i_grid
    integer :: threshold
    real(8) :: f
    real(8) :: prob
    real(8) :: cutoff
    real(8) :: uvw_new(3) ! outgoing uvw for iso-in-lab scattering
    real(8) :: uvw_old(3) ! incoming uvw for iso-in-lab scattering
    real(8) :: phi        ! azimuthal angle for iso-in-lab scattering
    real(8) :: kT         ! temperature in eV
    logical :: sampled    ! whether or not a reaction type has been sampled
    type(Nuclide),  pointer :: nuc

    ! copy incoming direction
    uvw_old(:) = p % coord(1) % uvw

    ! Get pointer to nuclide and grid index/interpolation factor
    nuc    => nuclides(i_nuclide)
    i_temp =  micro_xs(i_nuclide) % index_temp
    i_grid =  micro_xs(i_nuclide) % index_grid
    f      =  micro_xs(i_nuclide) % interp_factor

    ! For tallying purposes, this routine might be called directly. In that
    ! case, we need to sample a reaction via the cutoff variable
    cutoff = prn() * (micro_xs(i_nuclide) % total - &
         micro_xs(i_nuclide) % absorption)
    sampled = .false.

    ! Calculate elastic cross section if it wasn't precalculated
    if (micro_xs(i_nuclide) % elastic == CACHE_INVALID) then
      call nuc % calculate_elastic_xs(micro_xs(i_nuclide))
    end if

    prob = micro_xs(i_nuclide) % elastic - micro_xs(i_nuclide) % thermal
    if (prob > cutoff) then
      ! =======================================================================
      ! NON-S(A,B) ELASTIC SCATTERING

      ! Determine temperature
      if (nuc % mp_present) then
        kT = p % sqrtkT**2
      else
        kT = nuc % kTs(micro_xs(i_nuclide) % index_temp)
      end if

      ! Perform collision physics for elastic scattering
      call elastic_scatter(i_nuclide, nuc % reactions(1), kT, p % E, &
                           p % coord(1) % uvw, p % mu, p % wgt)

      p % event_MT = ELASTIC
      sampled = .true.
    end if

    prob = micro_xs(i_nuclide) % elastic
    if (prob > cutoff .and. .not. sampled) then
      ! =======================================================================
      ! S(A,B) SCATTERING

      call sab_scatter(i_nuclide, micro_xs(i_nuclide) % index_sab, p % E, &
                       p % coord(1) % uvw, p % mu)

      p % event_MT = ELASTIC
      sampled = .true.
    end if

    if (.not. sampled) then
      ! =======================================================================
      ! INELASTIC SCATTERING

      j = 0
      do while (prob < cutoff)
        j = j + 1
        i = nuc % index_inelastic_scatter(j)

        ! Check to make sure inelastic scattering reaction sampled
        if (i > size(nuc % reactions)) then
          call particle_write_restart(p)
          call fatal_error("Did not sample any reaction for nuclide " &
               &// trim(nuc % name))
        end if

        associate (rx => nuc % reactions(i))
          ! if energy is below threshold for this reaction, skip it
          threshold = rx % xs_threshold(i_temp)
          if (i_grid < threshold) cycle

          ! add to cumulative probability
          prob = prob + ((ONE - f)*rx % xs(i_temp, i_grid - threshold + 1) &
               + f*(rx % xs(i_temp, i_grid - threshold + 2)))
        end associate
      end do

      ! Perform collision physics for inelastic scattering
      call inelastic_scatter(nuc, nuc%reactions(i), p)
      p % event_MT = nuc % reactions(i) % MT

    end if

    ! Set event component
    p % event = EVENT_SCATTER

    ! Sample new outgoing angle for isotropic-in-lab scattering
    associate (mat => materials(p % material))
      if (mat % has_isotropic_nuclides) then
        if (materials(p % material) % p0(i_nuc_mat)) then
          ! Sample isotropic-in-lab outgoing direction
          uvw_new(1) = TWO * prn() - ONE
          phi = TWO * PI * prn()
          uvw_new(2) = cos(phi) * sqrt(ONE - uvw_new(1)*uvw_new(1))
          uvw_new(3) = sin(phi) * sqrt(ONE - uvw_new(1)*uvw_new(1))
          p % mu = dot_product(uvw_old, uvw_new)

          ! Change direction of particle
          p % coord(1) % uvw = uvw_new
        end if
      end if
    end associate

  end subroutine scatter

!===============================================================================
! ELASTIC_SCATTER treats the elastic scattering of a neutron with a
! target.
!===============================================================================

  subroutine elastic_scatter(i_nuclide, rxn, kT, E, uvw, mu_lab, wgt)
    integer, intent(in)     :: i_nuclide
    type(Reaction), intent(in) :: rxn
    real(8), intent(in)     :: kT      ! temperature in eV
    real(8), intent(inout)  :: E
    real(8), intent(inout)  :: uvw(3)
    real(8), intent(out)    :: mu_lab
    real(8), intent(inout)  :: wgt

    real(8) :: awr       ! atomic weight ratio of target
    real(8) :: mu_cm     ! cosine of polar angle in center-of-mass
    real(8) :: vel       ! magnitude of velocity
    real(8) :: v_n(3)    ! velocity of neutron
    real(8) :: v_cm(3)   ! velocity of center-of-mass
    real(8) :: v_t(3)    ! velocity of target nucleus
    real(8) :: uvw_cm(3) ! directional cosines in center-of-mass
    type(Nuclide), pointer :: nuc

    ! get pointer to nuclide
    nuc => nuclides(i_nuclide)

    vel = sqrt(E)
    awr = nuc % awr

    ! Neutron velocity in LAB
    v_n = vel * uvw

    ! Sample velocity of target nucleus
    if (.not. micro_xs(i_nuclide) % use_ptable) then
      call sample_target_velocity(nuc, v_t, E, uvw, v_n, wgt, &
           micro_xs(i_nuclide) % elastic, kT)
    else
      v_t = ZERO
    end if

    ! Velocity of center-of-mass
    v_cm = (v_n + awr*v_t)/(awr + ONE)

    ! Transform to CM frame
    v_n = v_n - v_cm

    ! Find speed of neutron in CM
    vel = sqrt(dot_product(v_n, v_n))

    ! Sample scattering angle
    mu_cm = rxn % sample_elastic_mu(E)

    ! Determine direction cosines in CM
    uvw_cm = v_n/vel

    ! Rotate neutron velocity vector to new angle -- note that the speed of the
    ! neutron in CM does not change in elastic scattering. However, the speed
    ! will change when we convert back to LAB
    v_n = vel * rotate_angle(uvw_cm, mu_cm)

    ! Transform back to LAB frame
    v_n = v_n + v_cm

    E = dot_product(v_n, v_n)
    vel = sqrt(E)

    ! compute cosine of scattering angle in LAB frame by taking dot product of
    ! neutron's pre- and post-collision angle
    mu_lab = dot_product(uvw, v_n) / vel

    ! Set energy and direction of particle in LAB frame
    uvw = v_n / vel

    ! Because of floating-point roundoff, it may be possible for mu_lab to be
    ! outside of the range [-1,1). In these cases, we just set mu_lab to exactly
    ! -1 or 1

    if (abs(mu_lab) > ONE) mu_lab = sign(ONE,mu_lab)

  end subroutine elastic_scatter

!===============================================================================
! SAB_SCATTER performs thermal scattering of a particle with a bound scatterer
! according to a specified S(a,b) table.
!===============================================================================

  subroutine sab_scatter(i_nuclide, i_sab, E, uvw, mu)
    integer, intent(in)     :: i_nuclide ! index in micro_xs
    integer, intent(in)     :: i_sab     ! index in sab_tables
    real(8), intent(inout)  :: E         ! incoming/outgoing energy
    real(8), intent(inout)  :: uvw(3)    ! directional cosines
    real(8), intent(out)    :: mu        ! scattering cosine

    real(C_DOUBLE) :: E_out
    type(C_PTR) :: ptr

    ! Sample from C++ side
    ptr = C_LOC(micro_xs(i_nuclide))
    call sab_tables(i_sab) % sample(ptr, E, E_out, mu)

    ! Set energy to outgoing, change direction of particle
    E = E_out
    uvw = rotate_angle(uvw, mu)
  end subroutine sab_scatter

!===============================================================================
! SAMPLE_TARGET_VELOCITY samples the target velocity. The constant cross section
! free gas model is the default method. Methods for correctly accounting
! for the energy dependence of cross sections in treating resonance elastic
! scattering such as the DBRC, WCM, and a new, accelerated scheme are also
! implemented here.
!===============================================================================

  subroutine sample_target_velocity(nuc, v_target, E, uvw, v_neut, wgt, xs_eff, kT)
    type(Nuclide), intent(in) :: nuc ! target nuclide at temperature T
    real(8), intent(out)   :: v_target(3) ! target velocity
    real(8), intent(in)    :: E           ! particle energy
    real(8), intent(in)    :: uvw(3)      ! direction cosines
    real(8), intent(in)    :: v_neut(3)   ! neutron velocity
    real(8), intent(inout) :: wgt         ! particle weight
    real(8), intent(in)    :: xs_eff      ! effective elastic xs at temperature T
    real(8), intent(in)    :: kT          ! equilibrium temperature of target in eV

    real(8) :: awr     ! target/neutron mass ratio
    real(8) :: E_rel   ! trial relative energy
    real(8) :: xs_0K   ! 0K xs at E_rel
    real(8) :: wcf     ! weight correction factor
    real(8) :: E_red   ! reduced energy (same as used by Cullen in SIGMA1)
    real(8) :: E_low   ! lowest practical relative energy
    real(8) :: E_up    ! highest practical relative energy
    real(8) :: E_t     ! trial target energy
    real(8) :: xs_max  ! max 0K xs over practical relative energies
    real(8) :: xs_low  ! 0K xs at lowest practical relative energy
    real(8) :: xs_up   ! 0K xs at highest practical relative energy
    real(8) :: m       ! slope for interpolation
    real(8) :: R       ! rejection criterion for DBRC / target speed
    real(8) :: cdf_low ! xs cdf at lowest practical relative energy
    real(8) :: cdf_up  ! xs cdf at highest practical relative energy
    real(8) :: cdf_rel ! trial xs cdf value
    real(8) :: mu      ! cosine between neutron and target velocities

    integer :: i_E_low ! 0K index to lowest practical relative energy
    integer :: i_E_up  ! 0K index to highest practical relative energy
    integer :: i_E_rel ! index to trial relative energy
    integer :: n_grid  ! number of energies on 0K grid

    integer :: sampling_method ! method of target velocity sampling

    awr = nuc % awr

    ! check if nuclide is a resonant scatterer
    if (nuc % resonant) then

      ! sampling method to use
      sampling_method = res_scat_method

      ! upper resonance scattering energy bound (target is at rest above this E)
      if (E > res_scat_energy_max) then
        v_target = ZERO
        return

      ! lower resonance scattering energy bound (should be no resonances below)
      else if (E < res_scat_energy_min) then
        sampling_method = RES_SCAT_CXS
      end if

    ! otherwise, use free gas model
    else
      if (E >= FREE_GAS_THRESHOLD * kT .and. awr > ONE) then
        v_target = ZERO
        return
      else
        sampling_method = RES_SCAT_CXS
      end if
    end if

    ! use appropriate target velocity sampling method
    select case (sampling_method)
    case (RES_SCAT_CXS)

      ! sample target velocity with the constant cross section (cxs) approx.
      call sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)

    case (RES_SCAT_WCM)

      ! sample target velocity with the constant cross section (cxs) approx.
      call sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)

      ! adjust weight as prescribed by the weight correction method (wcm)
      E_rel = dot_product((v_neut - v_target), (v_neut - v_target))
      xs_0K = elastic_xs_0K(E_rel, nuc)
      wcf = xs_0K / xs_eff
      wgt = wcf * wgt

    case (RES_SCAT_DBRC, RES_SCAT_ARES)
      E_red = sqrt(awr * E / kT)
      E_low = max(ZERO, E_red - FOUR)**2 * kT / awr
      E_up  = (E_red + FOUR)**2 * kT / awr

      ! find lower and upper energy bound indices
      ! lower index
      n_grid = size(nuc % energy_0K)
      if (E_low < nuc % energy_0K(1)) then
        i_E_low = 1
      elseif (E_low > nuc % energy_0K(n_grid)) then
        i_E_low = n_grid - 1
      else
        i_E_low = binary_search(nuc % energy_0K, n_grid, E_low)
      end if

      ! upper index
      if (E_up < nuc % energy_0K(1)) then
        i_E_up = 1
      elseif (E_up > nuc % energy_0K(n_grid)) then
        i_E_up = n_grid - 1
      else
        i_E_up = binary_search(nuc % energy_0K, n_grid, E_up)
      end if

      if (i_E_up == i_E_low) then
        ! Handle degenerate case -- if the upper/lower bounds occur for the same
        ! index, then using cxs is probably a good approximation
        call sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)

      else
        if (sampling_method == RES_SCAT_DBRC) then
          ! interpolate xs since we're not exactly at the energy indices
          xs_low = nuc % elastic_0K(i_E_low)
          m = (nuc % elastic_0K(i_E_low + 1) - xs_low) &
               / (nuc % energy_0K(i_E_low + 1) - nuc % energy_0K(i_E_low))
          xs_low = xs_low + m * (E_low - nuc % energy_0K(i_E_low))
          xs_up = nuc % elastic_0K(i_E_up)
          m = (nuc % elastic_0K(i_E_up + 1) - xs_up) &
               / (nuc % energy_0K(i_E_up + 1) - nuc % energy_0K(i_E_up))
          xs_up = xs_up + m * (E_up - nuc % energy_0K(i_E_up))

          ! get max 0K xs value over range of practical relative energies
          xs_max = max(xs_low, &
               maxval(nuc % elastic_0K(i_E_low + 1 : i_E_up)), xs_up)

          DBRC_REJECT_LOOP: do
            TARGET_ENERGY_LOOP: do
              ! sample target velocity with the constant cross section (cxs) approx.
              call sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)
              E_rel = dot_product((v_neut - v_target), (v_neut - v_target))
              if (E_rel < E_up) exit TARGET_ENERGY_LOOP
            end do TARGET_ENERGY_LOOP

            ! perform Doppler broadening rejection correction (dbrc)
            xs_0K = elastic_xs_0K(E_rel, nuc)
            R = xs_0K / xs_max
            if (prn() < R) exit DBRC_REJECT_LOOP
          end do DBRC_REJECT_LOOP

        elseif (sampling_method == RES_SCAT_ARES) then
          ! interpolate xs CDF since we're not exactly at the energy indices
          ! cdf value at lower bound attainable energy
          m = (nuc % xs_cdf(i_E_low) - nuc % xs_cdf(i_E_low - 1)) &
               / (nuc % energy_0K(i_E_low + 1) - nuc % energy_0K(i_E_low))
          cdf_low = nuc % xs_cdf(i_E_low - 1) &
               + m * (E_low - nuc % energy_0K(i_E_low))
          if (E_low <= nuc % energy_0K(1)) cdf_low = ZERO

          ! cdf value at upper bound attainable energy
          m = (nuc % xs_cdf(i_E_up) - nuc % xs_cdf(i_E_up - 1)) &
               / (nuc % energy_0K(i_E_up + 1) - nuc % energy_0K(i_E_up))
          cdf_up = nuc % xs_cdf(i_E_up - 1) &
               + m * (E_up - nuc % energy_0K(i_E_up))

          ARES_REJECT_LOOP: do

            ! directly sample Maxwellian
            E_t = -kT * log(prn())

            ! sample a relative energy using the xs cdf
            cdf_rel = cdf_low + prn() * (cdf_up - cdf_low)
            i_E_rel = binary_search(nuc % xs_cdf(i_E_low-1:i_E_up), &
                 i_E_up - i_E_low + 2, cdf_rel)
            E_rel = nuc % energy_0K(i_E_low + i_E_rel - 1)
            m = (nuc % xs_cdf(i_E_low + i_E_rel - 1) &
                 - nuc % xs_cdf(i_E_low + i_E_rel - 2)) &
                 / (nuc % energy_0K(i_E_low + i_E_rel) &
                 -  nuc % energy_0K(i_E_low + i_E_rel - 1))
            E_rel = E_rel + (cdf_rel - nuc % xs_cdf(i_E_low + i_E_rel - 2)) / m

            ! perform rejection sampling on cosine between
            ! neutron and target velocities
            mu = (E_t + awr * (E - E_rel)) / (TWO * sqrt(awr * E * E_t))

            if (abs(mu) < ONE) then
              ! set and accept target velocity
              E_t = E_t / awr
              v_target = sqrt(E_t) * rotate_angle(uvw, mu)
              exit ARES_REJECT_LOOP
            end if
          end do ARES_REJECT_LOOP
        end if
      end if
    end select

  end subroutine sample_target_velocity

!===============================================================================
! SAMPLE_CXS_TARGET_VELOCITY samples a target velocity based on the free gas
! scattering formulation, used by most Monte Carlo codes, in which cross section
! is assumed to be constant in energy. Excellent documentation for this method
! can be found in FRA-TM-123.
!===============================================================================

  subroutine sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)
    type(Nuclide), intent(in) :: nuc ! target nuclide at temperature
    real(8), intent(out)         :: v_target(3)
    real(8), intent(in)          :: E
    real(8), intent(in)          :: uvw(3)
    real(8), intent(in)          :: kT      ! equilibrium temperature of target in eV

    real(8) :: awr         ! target/neutron mass ratio
    real(8) :: alpha       ! probability of sampling f2 over f1
    real(8) :: mu          ! cosine of angle between neutron and target vel
    real(8) :: r1, r2      ! pseudo-random numbers
    real(8) :: c           ! cosine used in maxwell sampling
    real(8) :: accept_prob ! probability of accepting combination of vt and mu
    real(8) :: beta_vn     ! beta * speed of neutron
    real(8) :: beta_vt     ! beta * speed of target
    real(8) :: beta_vt_sq  ! (beta * speed of target)^2
    real(8) :: vt          ! speed of target

    awr = nuc % awr

    beta_vn = sqrt(awr * E / kT)
    alpha = ONE/(ONE + sqrt(pi)*beta_vn/TWO)

    do
      ! Sample two random numbers
      r1 = prn()
      r2 = prn()

      if (prn() < alpha) then
        ! With probability alpha, we sample the distribution p(y) =
        ! y*e^(-y). This can be done with sampling scheme C45 frmo the Monte
        ! Carlo sampler

        beta_vt_sq = -log(r1*r2)

      else
        ! With probability 1-alpha, we sample the distribution p(y) = y^2 *
        ! e^(-y^2). This can be done with sampling scheme C61 from the Monte
        ! Carlo sampler

        c = cos(PI/TWO * prn())
        beta_vt_sq = -log(r1) - log(r2)*c*c
      end if

      ! Determine beta * vt
      beta_vt = sqrt(beta_vt_sq)

      ! Sample cosine of angle between neutron and target velocity
      mu = TWO*prn() - ONE

      ! Determine rejection probability
      accept_prob = sqrt(beta_vn*beta_vn + beta_vt_sq - 2*beta_vn*beta_vt*mu) &
           /(beta_vn + beta_vt)

      ! Perform rejection sampling on vt and mu
      if (prn() < accept_prob) exit
    end do

    ! Determine speed of target nucleus
    vt = sqrt(beta_vt_sq*kT/awr)

    ! Determine velocity vector of target nucleus based on neutron's velocity
    ! and the sampled angle between them
    v_target = vt * rotate_angle(uvw, mu)

  end subroutine sample_cxs_target_velocity

!===============================================================================
! CREATE_FISSION_SITES determines the average total, prompt, and delayed
! neutrons produced from fission and creates appropriate bank sites.
!===============================================================================

  subroutine create_fission_sites(p, i_nuclide, i_reaction, bank_array, size_bank)
    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: i_nuclide
    integer,        intent(in)    :: i_reaction
    type(Bank),     intent(inout) :: bank_array(:)
    integer(8),     intent(inout) :: size_bank

    integer :: nu_d(MAX_DELAYED_GROUPS) ! number of delayed neutrons born
    integer :: i                        ! loop index
    integer :: nu                       ! actual number of neutrons produced
    integer :: mesh_bin                 ! mesh bin for source site
    real(8) :: nu_t                     ! total nu
    real(8) :: weight                   ! weight adjustment for ufs method
    type(Nuclide),  pointer :: nuc

    ! Get pointers
    nuc => nuclides(i_nuclide)

    ! TODO: Heat generation from fission

    ! If uniform fission source weighting is turned on, we increase of decrease
    ! the expected number of fission sites produced

    if (ufs) then
      associate (m => meshes(index_ufs_mesh))
        ! Determine indices on ufs mesh for current location
        call m % get_bin(p % coord(1) % xyz, mesh_bin)
        if (mesh_bin == NO_BIN_FOUND) then
          call particle_write_restart(p)
          call fatal_error("Source site outside UFS mesh!")
        end if

        if (source_frac(1, mesh_bin) /= ZERO) then
          weight = m % volume_frac / source_frac(1, mesh_bin)
        else
          weight = ONE
        end if
      end associate
    else
      weight = ONE
    end if

    ! Determine expected number of neutrons produced
    nu_t = p % wgt / keff * weight * micro_xs(i_nuclide) % nu_fission / &
         micro_xs(i_nuclide) % total

    ! Sample number of neutrons produced
    if (prn() > nu_t - int(nu_t)) then
      nu = int(nu_t)
    else
      nu = int(nu_t) + 1
    end if

    ! Check for bank size getting hit. For fixed source calculations, this is a
    ! fatal error. For eigenvalue calculations, it just means that k-effective
    ! was too high for a single batch.
    if (size_bank + nu > size(bank_array)) then
      if (run_mode == MODE_FIXEDSOURCE) then
        call fatal_error("Secondary particle bank size limit reached. If you &
             &are running a subcritical multiplication problem, k-effective &
             &may be too close to one.")
      else
        if (master) call warning("Maximum number of sites in fission bank &
             &reached. This can result in irreproducible results using different &
             &numbers of processes/threads.")
      end if
    end if

    ! Bank source neutrons
    if (nu == 0 .or. size_bank == size(bank_array)) return

    ! Initialize counter of delayed neutrons encountered for each delayed group
    ! to zero.
    nu_d(:) = 0

    p % fission = .true. ! Fission neutrons will be banked
    do i = int(size_bank,4) + 1, int(min(size_bank + nu, int(size(bank_array),8)),4)
      ! Bank source neutrons by copying particle data
      bank_array(i) % xyz = p % coord(1) % xyz

      ! Set particle as neutron
      bank_array(i) % particle = NEUTRON

      ! Set weight of fission bank site
      bank_array(i) % wgt = ONE/weight

      ! Sample delayed group and angle/energy for fission reaction
      call sample_fission_neutron(nuc, nuc % reactions(i_reaction), &
           p % E, bank_array(i))

      ! Set delayed group on particle too
      p % delayed_group = bank_array(i) % delayed_group

      ! Increment the number of neutrons born delayed
      if (p % delayed_group > 0) then
        nu_d(p % delayed_group) = nu_d(p % delayed_group) + 1
      end if
    end do

    ! increment number of bank sites
    size_bank = min(size_bank + nu, int(size(bank_array),8))

    ! Store total and delayed weight banked for analog fission tallies
    p % n_bank   = nu
    p % wgt_bank = nu/weight
    p % n_delayed_bank(:) = nu_d(:)

  end subroutine create_fission_sites

!===============================================================================
! SAMPLE_FISSION_NEUTRON
!===============================================================================

  subroutine sample_fission_neutron(nuc, rxn, E_in, site)
    type(Nuclide),  intent(in)    :: nuc
    type(Reaction), intent(in)    :: rxn
    real(8),        intent(in)    :: E_in
    type(Bank),     intent(inout) :: site

    integer :: group        ! index on nu energy grid / precursor group
    integer :: n_sample     ! number of resamples
    real(8) :: nu_t         ! total nu
    real(8) :: nu_d         ! delayed nu
    real(8) :: beta         ! delayed neutron fraction
    real(8) :: xi           ! random number
    real(8) :: yield        ! delayed neutron precursor yield
    real(8) :: prob         ! cumulative probability
    real(8) :: mu           ! cosine of scattering angle
    real(8) :: phi          ! azimuthal angle

    ! Sample cosine of angle -- fission neutrons are always emitted
    ! isotropically. Sometimes in ACE data, fission reactions actually have
    ! an angular distribution listed, but for those that do, it's simply just
    ! a uniform distribution in mu
    mu = TWO * prn() - ONE

    ! Sample azimuthal angle uniformly in [0,2*pi)
    phi = TWO*PI*prn()
    site % uvw(1) = mu
    site % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
    site % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

    ! Determine total nu, delayed nu, and delayed neutron fraction
    nu_t = nuc % nu(E_in, EMISSION_TOTAL)
    nu_d = nuc % nu(E_in, EMISSION_DELAYED)
    beta = nu_d / nu_t

    if (prn() < beta) then
      ! ====================================================================
      ! DELAYED NEUTRON SAMPLED

      ! sampled delayed precursor group
      xi = prn()*nu_d
      prob = ZERO
      do group = 1, nuc % n_precursor

        ! determine delayed neutron precursor yield for group j
        yield = rxn % product_yield(1 + group, E_in)

        ! Check if this group is sampled
        prob = prob + yield
        if (xi < prob) exit
      end do

      ! if the sum of the probabilities is slightly less than one and the
      ! random number is greater, j will be greater than nuc %
      ! n_precursor -- check for this condition
      group = min(group, nuc % n_precursor)

      ! set the delayed group for the particle born from fission
      site % delayed_group = group

      n_sample = 0
      do
        ! sample from energy/angle distribution -- note that mu has already been
        ! sampled above and doesn't need to be resampled
        call rxn % product_sample(1 + group, E_in, site % E, mu)

        ! resample if energy is greater than maximum neutron energy
        if (site % E < energy_max(NEUTRON)) exit

        ! check for large number of resamples
        n_sample = n_sample + 1
        if (n_sample == MAX_SAMPLE) then
          ! call particle_write_restart(p)
          call fatal_error("Resampled energy distribution maximum number of " &
               // "times for nuclide " // nuc % name)
        end if
      end do

    else
      ! ====================================================================
      ! PROMPT NEUTRON SAMPLED

      ! set the delayed group for the particle born from fission to 0
      site % delayed_group = 0

      ! sample from prompt neutron energy distribution
      n_sample = 0
      do
        call rxn % product_sample(1, E_in, site % E, mu)

        ! resample if energy is greater than maximum neutron energy
        if (site % E < energy_max(NEUTRON)) exit

        ! check for large number of resamples
        n_sample = n_sample + 1
        if (n_sample == MAX_SAMPLE) then
          ! call particle_write_restart(p)
          call fatal_error("Resampled energy distribution maximum number of " &
               // "times for nuclide " // nuc % name)
        end if
      end do
    end if

  end subroutine sample_fission_neutron

!===============================================================================
! INELASTIC_SCATTER handles all reactions with a single secondary neutron (other
! than fission), i.e. level scattering, (n,np), (n,na), etc.
!===============================================================================

  subroutine inelastic_scatter(nuc, rxn, p)
    type(Nuclide), intent(in)    :: nuc
    type(Reaction),   intent(in)    :: rxn
    type(Particle),   intent(inout) :: p

    integer :: i      ! loop index
    real(8) :: E      ! energy in lab (incoming/outgoing)
    real(8) :: mu     ! cosine of scattering angle in lab
    real(8) :: A      ! atomic weight ratio of nuclide
    real(8) :: E_in   ! incoming energy
    real(8) :: E_cm   ! outgoing energy in center-of-mass
    real(8) :: yield  ! neutron yield

    ! copy energy of neutron
    E_in = p % E

    ! sample outgoing energy and scattering cosine
    call rxn % product_sample(1, E_in, E, mu)

    ! if scattering system is in center-of-mass, transfer cosine of scattering
    ! angle and outgoing energy from CM to LAB
    if (rxn % scatter_in_cm) then
      E_cm = E

      ! determine outgoing energy in lab
      A = nuc%awr
      E = E_cm + (E_in + TWO * mu * (A+ONE) * sqrt(E_in * E_cm)) &
           / ((A+ONE)*(A+ONE))

      ! determine outgoing angle in lab
      mu = mu * sqrt(E_cm/E) + ONE/(A+ONE) * sqrt(E_in/E)
    end if

    ! Because of floating-point roundoff, it may be possible for mu to be
    ! outside of the range [-1,1). In these cases, we just set mu to exactly -1
    ! or 1
    if (abs(mu) > ONE) mu = sign(ONE,mu)

    ! Set outgoing energy and scattering angle
    p % E = E
    p % mu = mu

    ! change direction of particle
    p % coord(1) % uvw = rotate_angle(p % coord(1) % uvw, mu)

    ! evaluate yield
    yield = rxn % product_yield(1, E_in)
    if (mod(yield, ONE) == ZERO) then
      ! If yield is integral, create exactly that many secondary particles
      do i = 1, nint(yield) - 1
        call particle_create_secondary(p, p % coord(1) % uvw, p % E, &
             NEUTRON, run_CE=.true._C_BOOL)
      end do
    else
      ! Otherwise, change weight of particle based on yield
      p % wgt = yield * p % wgt
    end if

  end subroutine inelastic_scatter

!===============================================================================
! SAMPLE_SECONDARY_PHOTONS
!===============================================================================

  subroutine sample_secondary_photons(p, i_nuclide)
    type(Particle), intent(inout) :: p
    integer, intent(in)           :: i_nuclide

    integer :: i_reaction   ! index in nuc % reactions array
    integer :: i_product    ! index in nuc % reactions % products array

    real(8) :: nu_t
    real(8) :: mu
    real(8) :: E
    real(8) :: uvw(3)
    integer :: nu
    integer :: i

    ! Sample the number of photons produced
    nu_t =  p % wgt * micro_xs(i_nuclide) % photon_prod / &
         micro_xs(i_nuclide) % total
    if (prn() > nu_t - int(nu_t)) then
      nu = int(nu_t)
    else
      nu = int(nu_t) + 1
    end if

    ! Sample each secondary photon
    do i = 1, nu

      ! Sample the reaction and product
      call sample_photon_product(i_nuclide, p % E, i_reaction, i_product)

      ! Sample the outgoing energy and angle
      call nuclides(i_nuclide) % reactions(i_reaction) % &
           product_sample(i_product, p % E, E, mu)

      ! Sample the new direction
      uvw = rotate_angle(p % coord(1) % uvw, mu)

      ! Create the secondary photon
      call particle_create_secondary(p, uvw, E, PHOTON, run_CE=.true._C_BOOL)
    end do

  end subroutine sample_secondary_photons

end module physics
