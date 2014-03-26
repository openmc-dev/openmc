module physics

  use ace_header,             only: Nuclide, Reaction, DistEnergy
  use constants
  use endf,                   only: reaction_name
  use error,                  only: fatal_error, warning
  use fission,                only: nu_total, nu_delayed
  use global
  use interpolation,          only: interpolate_tab1
  use material_header,        only: Material
  use math,                   only: maxwell_spectrum, watt_spectrum
  use mesh,                   only: get_mesh_indices
  use output,                 only: write_message
  use particle_header,        only: Particle
  use particle_restart_write, only: write_particle_restart
  use random_lcg,             only: prn
  use search,                 only: binary_search
  use string,                 only: to_str

  implicit none

! TODO: Figure out how to write particle restart files in sample_angle,
! sample_energy, etc.

contains

!===============================================================================
! COLLISION samples a nuclide and reaction and then calls the appropriate
! routine for that reaction
!===============================================================================

  subroutine collision(p)

    type(Particle), intent(inout) :: p

    ! Store pre-collision particle properties
    p % last_wgt = p % wgt
    p % last_E   = p % E

    ! Add to collision counter for particle
    p % n_collision = p % n_collision + 1

    ! Sample nuclide/reaction for the material the particle is in
    call sample_reaction(p)

    ! Display information about collision
    if (verbosity >= 10 .or. trace) then
      message = "    " // trim(reaction_name(p % event_MT)) // " with " // &
           trim(adjustl(nuclides(p % event_nuclide) % name)) // &
           ". Energy = " // trim(to_str(p % E * 1e6_8)) // " eV."
      call write_message()
    end if

    ! check for very low energy
    if (p % E < 1.0e-100_8) then
      p % alive = .false.
      message = "Killing neutron with extremely low energy"
      call warning()
    end if

  end subroutine collision

!===============================================================================
! SAMPLE_REACTION samples a nuclide based on the macroscopic cross sections for
! each nuclide within a material and then samples a reaction for that nuclide
! and calls the appropriate routine to process the physics. Note that there is
! special logic when suvival biasing is turned on since fission and
! disappearance are treated implicitly.
!===============================================================================

  subroutine sample_reaction(p)

    type(Particle), intent(inout) :: p

    integer :: i_nuclide    ! index in nuclides array
    integer :: i_reaction   ! index in nuc % reactions array
    type(Nuclide), pointer, save :: nuc => null()
!$omp threadprivate(nuc)

    i_nuclide = sample_nuclide(p, 'total  ')

    ! Get pointer to table
    nuc => nuclides(i_nuclide)

    ! Save which nuclide particle had collision with
    p % event_nuclide = i_nuclide

    ! Create fission bank sites. Note that while a fission reaction is sampled,
    ! it never actually "happens", i.e. the weight of the particle does not
    ! change when sampling fission sites. The following block handles all
    ! absorption (including fission)

    if (nuc % fissionable) then
      call sample_fission(i_nuclide, i_reaction)
      call create_fission_sites(p, i_nuclide, i_reaction)
    end if

    ! If survival biasing is being used, the following subroutine adjusts the
    ! weight of the particle. Otherwise, it checks to see if absorption occurs

    call absorption(p, i_nuclide)
    if (.not. p % alive) return

    ! Play russian roulette if survival biasing is turned on

    if (survival_biasing) then
      call russian_roulette(p)
      if (.not. p % alive) return
    end if

    ! Sample a scattering reaction and determine the secondary energy of the
    ! exiting neutron

    call scatter(p, i_nuclide)

  end subroutine sample_reaction

!===============================================================================
! SAMPLE_NUCLIDE
!===============================================================================

  function sample_nuclide(p, base) result(i_nuclide)

    type(Particle), intent(in) :: p
    character(7),   intent(in) :: base      ! which reaction to sample based on
    integer                    :: i_nuclide

    integer :: i
    real(8) :: prob
    real(8) :: cutoff
    real(8) :: atom_density ! atom density of nuclide in atom/b-cm
    real(8) :: sigma        ! microscopic total xs for nuclide
    type(Material), pointer, save :: mat => null()
!$omp threadprivate(mat)

    ! Get pointer to current material
    mat => materials(p % material)

    ! Sample cumulative distribution function
    select case (base)
    case ('total')
      cutoff = prn() * material_xs % total
    case ('scatter')
      cutoff = prn() * material_xs % total - material_xs % absorption
    case ('fission')
      cutoff = prn() * material_xs % fission
    end select

    i = 0
    prob = ZERO
    do while (prob < cutoff)
      i = i + 1

      ! Check to make sure that a nuclide was sampled
      if (i > mat % n_nuclides) then
        call write_particle_restart(p)
        message = "Did not sample any nuclide during collision."
        call fatal_error()
      end if

      ! Find atom density
      i_nuclide    = mat % nuclide(i)
      atom_density = mat % atom_density(i)

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

  end function sample_nuclide

!===============================================================================
! SAMPLE_FISSION
!===============================================================================

  subroutine sample_fission(i_nuclide, i_reaction)

    integer, intent(in)  :: i_nuclide  ! index in nuclides array
    integer, intent(out) :: i_reaction ! index in nuc % reactions array

    integer :: i
    integer :: i_grid
    real(8) :: f
    real(8) :: prob
    real(8) :: cutoff
    type(Nuclide),  pointer, save :: nuc => null()
    type(Reaction), pointer, save :: rxn => null()
!$omp threadprivate(nuc, rxn)

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

    ! Get grid index and interpolatoin factor and sample fission cdf
    i_grid = micro_xs(i_nuclide) % index_grid
    f      = micro_xs(i_nuclide) % interp_factor
    cutoff = prn() * micro_xs(i_nuclide) % fission
    prob   = ZERO

    ! Loop through each partial fission reaction type

    FISSION_REACTION_LOOP: do i = 1, nuc % n_fission
      i_reaction = nuc % index_fission(i)
      rxn => nuc % reactions(i_reaction)

      ! if energy is below threshold for this reaction, skip it
      if (i_grid < rxn % threshold) cycle

      ! add to cumulative probability
      prob = prob + ((ONE - f)*rxn%sigma(i_grid - rxn%threshold + 1) &
           + f*(rxn%sigma(i_grid - rxn%threshold + 2)))

      ! Create fission bank sites if fission occus
      if (prob > cutoff) exit FISSION_REACTION_LOOP
    end do FISSION_REACTION_LOOP

  end subroutine sample_fission

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
!$omp critical
      global_tallies(K_ABSORPTION) % value = &
           global_tallies(K_ABSORPTION) % value + p % absorb_wgt * &
           micro_xs(i_nuclide) % nu_fission / micro_xs(i_nuclide) % absorption
!$omp end critical

    else
      ! See if disappearance reaction happens
      if (micro_xs(i_nuclide) % absorption > &
           prn() * micro_xs(i_nuclide) % total) then
        ! Score absorption estimate of keff
!$omp critical
        global_tallies(K_ABSORPTION) % value = &
             global_tallies(K_ABSORPTION) % value + p % wgt * &
             micro_xs(i_nuclide) % nu_fission / micro_xs(i_nuclide) % absorption
!$omp end critical

        p % alive = .false.
        p % event = EVENT_ABSORB
        p % event_MT = N_DISAPPEAR
      end if
    end if

  end subroutine absorption

!===============================================================================
! RUSSIAN_ROULETTE
!===============================================================================

  subroutine russian_roulette(p)

    type(Particle), intent(inout) :: p

    if (p % wgt < weight_cutoff) then
      if (prn() < p % wgt / weight_survive) then
        p % wgt = weight_survive
        p % last_wgt = p % wgt
      else
        p % wgt = ZERO
        p % alive = .false.
      end if
    end if

  end subroutine russian_roulette

!===============================================================================
! SCATTER
!===============================================================================

  subroutine scatter(p, i_nuclide)

    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: i_nuclide

    integer :: i
    integer :: i_grid
    real(8) :: f
    real(8) :: prob
    real(8) :: cutoff
    type(Nuclide),  pointer, save :: nuc => null()
    type(Reaction), pointer, save :: rxn => null()
!$omp threadprivate(nuc, rxn)

    ! Get pointer to nuclide and grid index/interpolation factor
    nuc    => nuclides(i_nuclide)
    i_grid =  micro_xs(i_nuclide) % index_grid
    f      =  micro_xs(i_nuclide) % interp_factor

    ! For tallying purposes, this routine might be called directly. In that
    ! case, we need to sample a reaction via the cutoff variable
    prob = ZERO
    cutoff = prn() * (micro_xs(i_nuclide) % total - &
         micro_xs(i_nuclide) % absorption)

    prob = prob + micro_xs(i_nuclide) % elastic
    if (prob > cutoff) then
      ! =======================================================================
      ! ELASTIC SCATTERING

      if (micro_xs(i_nuclide) % index_sab /= NONE) then
        ! S(a,b) scattering
        call sab_scatter(i_nuclide, micro_xs(i_nuclide) % index_sab, &
             p % E, p % coord0 % uvw, p % mu)

      else
        ! get pointer to elastic scattering reaction
        rxn => nuc % reactions(1)

        ! Perform collision physics for elastic scattering
        call elastic_scatter(i_nuclide, rxn, &
             p % E, p % coord0 % uvw, p % mu)
      end if

      p % event_MT = ELASTIC

    else
      ! =======================================================================
      ! INELASTIC SCATTERING

      ! note that indexing starts from 2 since nuc % reactions(1) is elastic
      ! scattering
      i = 1
      do while (prob < cutoff)
        i = i + 1

        ! Check to make sure inelastic scattering reaction sampled
        if (i > nuc % n_reaction) then
          call write_particle_restart(p)
          message = "Did not sample any reaction for nuclide " // &
               trim(nuc % name)
          call fatal_error()
        end if

        rxn => nuc % reactions(i)

        ! Skip fission reactions
        if (rxn % MT == N_FISSION .or. rxn % MT == N_F .or. rxn % MT == N_NF &
             .or. rxn % MT == N_2NF .or. rxn % MT == N_3NF) cycle

        ! some materials have gas production cross sections with MT > 200 that
        ! are duplicates. Also MT=4 is total level inelastic scattering which
        ! should be skipped
        if (rxn % MT >= 200 .or. rxn % MT == N_LEVEL) cycle

        ! if energy is below threshold for this reaction, skip it
        if (i_grid < rxn % threshold) cycle

        ! add to cumulative probability
        prob = prob + ((ONE - f)*rxn%sigma(i_grid - rxn%threshold + 1) &
             + f*(rxn%sigma(i_grid - rxn%threshold + 2)))
      end do

      ! Perform collision physics for inelastics scattering
      call inelastic_scatter(nuc, rxn, p % E, p % coord0 % uvw, &
           p % mu, p % wgt)
      p % event_MT = rxn % MT

    end if

    ! Set event component
    p % event = EVENT_SCATTER

  end subroutine scatter

!===============================================================================
! ELASTIC_SCATTER treats the elastic scattering of a neutron with a
! target.
!===============================================================================

  subroutine elastic_scatter(i_nuclide, rxn, E, uvw, mu_lab)

    integer, intent(in)     :: i_nuclide
    type(Reaction), pointer :: rxn
    real(8), intent(inout)  :: E
    real(8), intent(inout)  :: uvw(3)
    real(8), intent(out)    :: mu_lab

    real(8) :: awr       ! atomic weight ratio of target
    real(8) :: mu_cm     ! cosine of polar angle in center-of-mass
    real(8) :: vel       ! magnitude of velocity
    real(8) :: v_n(3)    ! velocity of neutron
    real(8) :: v_cm(3)   ! velocity of center-of-mass
    real(8) :: v_t(3)    ! velocity of target nucleus
    real(8) :: uvw_cm(3) ! directional cosines in center-of-mass
    type(Nuclide), pointer, save :: nuc => null()
!$omp threadprivate(nuc)

    ! get pointer to nuclide
    nuc => nuclides(i_nuclide)

    vel = sqrt(E)
    awr = nuc % awr

    ! Neutron velocity in LAB
    v_n = vel * uvw

    ! Sample velocity of target nucleus
    if (.not. micro_xs(i_nuclide) % use_ptable) then
      call sample_target_velocity(nuc, v_t, E, uvw)
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
    mu_cm = sample_angle(rxn, E)

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

    integer :: i            ! incoming energy bin
    integer :: j            ! outgoing energy bin
    integer :: k            ! outgoing cosine bin
    integer :: n_energy_out ! number of outgoing energy bins
    real(8) :: f            ! interpolation factor
    real(8) :: r            ! used for skewed sampling & continuous
    real(8) :: E_ij         ! outgoing energy j for E_in(i)
    real(8) :: E_i1j        ! outgoing energy j for E_in(i+1)
    real(8) :: mu_ijk       ! outgoing cosine k for E_in(i) and E_out(j)
    real(8) :: mu_i1jk      ! outgoing cosine k for E_in(i+1) and E_out(j)
    real(8) :: prob         ! probability for sampling Bragg edge
    type(SAlphaBeta), pointer, save :: sab => null()
    ! Following are needed only for SAB_SECONDARY_CONT scattering
    integer :: l              ! sampled incoming E bin (is i or i + 1)
    real(8) :: E_i_1, E_i_J   ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_J ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_J       ! endpoints interpolated between i and i+1
    real(8) :: E_l_j, E_l_j1  ! adjacent E on outgoing grid l
    real(8) :: p_l_j, p_l_j1  ! adjacent p on outgoing grid l
    real(8) :: c_j, c_j1      ! cumulative probability
    real(8) :: frac           ! interpolation factor on outgoing energy
    real(8) :: r1             ! RNG for outgoing energy

!$omp threadprivate(sab)

    ! Get pointer to S(a,b) table
    sab => sab_tables(i_sab)

    ! Determine whether inelastic or elastic scattering will occur
    if (prn() < micro_xs(i_nuclide) % elastic_sab / &
         micro_xs(i_nuclide) % elastic) then
      ! elastic scattering

      ! Get index and interpolation factor for elastic grid
      if (E < sab % elastic_e_in(1)) then
        i = 1
        f = ZERO
      else
        i = binary_search(sab % elastic_e_in, sab % n_elastic_e_in, E)
        f = (E - sab%elastic_e_in(i)) / &
             (sab%elastic_e_in(i+1) - sab%elastic_e_in(i))
      end if

      ! Select treatment based on elastic mode
      if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
        ! With this treatment, we interpolate between two discrete cosines
        ! corresponding to neighboring incoming energies. This is used for
        ! data derived in the incoherent approximation

        ! Sample outgoing cosine bin
        k = 1 + int(prn() * sab % n_elastic_mu)

        ! Determine outgoing cosine corresponding to E_in(i) and E_in(i+1)
        mu_ijk  = sab % elastic_mu(k,i)
        mu_i1jk = sab % elastic_mu(k,i+1)

        ! Cosine of angle between incoming and outgoing neutron
        mu = (1 - f)*mu_ijk + f*mu_i1jk

      elseif (sab % elastic_mode == SAB_ELASTIC_EXACT) then
        ! This treatment is used for data derived in the coherent
        ! approximation, i.e. for crystalline structures that have Bragg
        ! edges.

        ! Sample a Bragg edge between 1 and i
        prob = prn() * sab % elastic_P(i+1)
        if (prob < sab % elastic_P(1)) then
          k = 1
        else
          k = binary_search(sab % elastic_P(1:i+1), i+1, prob)
        end if

        ! Characteristic scattering cosine for this Bragg edge
        mu = ONE - TWO*sab % elastic_e_in(k) / E

      end if

      ! Outgoing energy is same as incoming energy -- no need to do anything

    else
      ! Perform inelastic calculations

      ! Get index and interpolation factor for inelastic grid
      if (E < sab % inelastic_e_in(1)) then
        i = 1
        f = ZERO
      else
        i = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, E)
        f = (E - sab%inelastic_e_in(i)) / &
             (sab%inelastic_e_in(i+1) - sab%inelastic_e_in(i))
      end if

      ! Now that we have an incoming energy bin, we need to determine the
      ! outgoing energy bin. This will depend on the "secondary energy
      ! mode". If the mode is 0, then the outgoing energy bin is chosen from a
      ! set of equally-likely bins. If the mode is 1, then the first
      ! two and last two bins are skewed to have lower probabilities than the
      ! other bins (0.1 for the first and last bins and 0.4 for the second and
      ! second to last bins, relative to a normal bin probability of 1).
      ! Finally, if the mode is 2, then a continuous distribution (with
      ! accompanying PDF and CDF is utilized)

      if ((sab % secondary_mode == SAB_SECONDARY_EQUAL) .or. &
          (sab % secondary_mode == SAB_SECONDARY_SKEWED)) then
        if (sab % secondary_mode == SAB_SECONDARY_EQUAL) then
          ! All bins equally likely

          j = 1 + int(prn() * sab % n_inelastic_e_out)
        elseif (sab % secondary_mode == SAB_SECONDARY_SKEWED) then
          ! Distribution skewed away from edge points

          ! Determine number of outgoing energy and angle bins
          n_energy_out = sab % n_inelastic_e_out

          r = prn() * (n_energy_out - 3)
          if (r > ONE) then
            ! equally likely N-4 middle bins
            j = int(r) + 2
          elseif (r > 0.6) then
            ! second to last bin has relative probability of 0.4
            j = n_energy_out - 1
          elseif (r > 0.5) then
            ! last bin has relative probability of 0.1
            j = n_energy_out
          elseif (r > 0.1) then
            ! second bin has relative probability of 0.4
            j = 2
          else
            ! first bin has relative probability of 0.1
            j = 1
          end if
        end if

        ! Determine outgoing energy corresponding to E_in(i) and E_in(i+1)
        E_ij  = sab % inelastic_e_out(j,i)
        E_i1j = sab % inelastic_e_out(j,i+1)

        ! Outgoing energy
        E = (1 - f)*E_ij + f*E_i1j

        ! Sample outgoing cosine bin
        k = 1 + int(prn() * sab % n_inelastic_mu)

        ! Determine outgoing cosine corresponding to E_in(i) and E_in(i+1)
        mu_ijk  = sab % inelastic_mu(k,j,i)
        mu_i1jk = sab % inelastic_mu(k,j,i+1)

        ! Cosine of angle between incoming and outgoing neutron
        mu = (1 - f)*mu_ijk + f*mu_i1jk

      else if (sab % secondary_mode == SAB_SECONDARY_CONT) then
        ! Continuous secondary energy - this is to be similar to
        ! Law 61 interpolation on outgoing energy

        ! Sample between ith and (i+1)th bin
        r = prn()
        if (f > r) then
          l = i + 1
        else
          l = i
        end if

        ! Determine endpoints on grid i
        n_energy_out = sab % inelastic_data(i) % n_e_out
        E_i_1 = sab % inelastic_data(i) % e_out(1)
        E_i_J = sab % inelastic_data(i) % e_out(n_energy_out)

        ! Determine endpoints on grid i + 1
        n_energy_out = sab % inelastic_data(i + 1) % n_e_out
        E_i1_1 = sab % inelastic_data(i + 1) % e_out(1)
        E_i1_J = sab % inelastic_data(i + 1) % e_out(n_energy_out)

        E_1 = E_i_1 + f * (E_i1_1 - E_i_1)
        E_J = E_i_J + f * (E_i1_J - E_i_J)

        ! Determine outgoing energy bin
        ! (First reset n_energy_out to the right value)
        n_energy_out = sab % inelastic_data(l) % n_e_out
        r1 = prn()
        c_j = sab % inelastic_data(l) % e_out_cdf(1)
        do j = 1, n_energy_out - 1
          c_j1 = sab % inelastic_data(l) % e_out_cdf(j + 1)
          if (r1 < c_j1) exit
          c_j = c_j1
        end do

        ! check to make sure k is <= n_energy_out - 1
        j = min(j, n_energy_out - 1)

        ! Get the data to interpolate between
        E_l_j = sab % inelastic_data(l) % e_out(j)
        p_l_j = sab % inelastic_data(l) % e_out_pdf(j)

        ! Next part assumes linear-linear interpolation in standard
        E_l_j1 = sab % inelastic_data(l) % e_out(j + 1)
        p_l_j1 = sab % inelastic_data(l) % e_out_pdf(j + 1)

        ! Find secondary energy (variable E)
        frac = (p_l_j1 - p_l_j) / (E_l_j1 - E_l_j)
        if (frac == ZERO) then
          E = E_l_j + (r1 - c_j) / p_l_j
        else
          E = E_l_j + (sqrt(max(ZERO, p_l_j * p_l_j + &
                       TWO * frac * (r1 - c_j))) - p_l_j) / frac
        end if

        ! Now interpolate between incident energy bins i and i + 1
        if (l == i) then
          E = E_1 + (E - E_i_1) * (E_J - E_1) / (E_i_J - E_i_1)
        else
          E = E_1 + (E - E_i1_1) * (E_J - E_1) / (E_i1_J - E_i1_1)
        end if

        ! Find angular distribution for closest outgoing energy bin
        if (r1 - c_j < c_j1 - r1) then
          j = j
        else
          j = j + 1
        end if

        ! Sample outgoing cosine bin
        k = 1 + int(prn() * sab % n_inelastic_mu)

        ! Will use mu from the randomly chosen incoming and closest outgoing
        ! energy bins
        mu = sab % inelastic_data(l) % mu(k, j)

      else
        message = "Invalid secondary energy mode on S(a,b) table " // &
             trim(sab % name)
      end if  ! (inelastic secondary energy treatment)
    end if  ! (elastic or inelastic)

    ! change direction of particle
    uvw = rotate_angle(uvw, mu)

  end subroutine sab_scatter

!===============================================================================
! SAMPLE_TARGET_VELOCITY samples the target velocity based on the free gas
! scattering formulation used by most Monte Carlo codes. Excellent documentation
! for this method can be found in FRA-TM-123.
!===============================================================================

  subroutine sample_target_velocity(nuc, v_target, E, uvw)

    type(Nuclide),  pointer :: nuc
    real(8), intent(out)    :: v_target(3)
    real(8), intent(in)     :: E
    real(8), intent(in)     :: uvw(3)

    real(8) :: kT          ! equilibrium temperature of target in MeV
    real(8) :: alpha       ! probability of sampling f2 over f1
    real(8) :: mu          ! cosine of angle between neutron and target vel
    real(8) :: r1, r2      ! pseudo-random numbers
    real(8) :: c           ! cosine used in maxwell sampling
    real(8) :: accept_prob ! probability of accepting combination of vt and mu
    real(8) :: beta_vn     ! beta * speed of neutron
    real(8) :: beta_vt     ! beta * speed of target
    real(8) :: beta_vt_sq  ! (beta * speed of target)^2
    real(8) :: vt          ! speed of target

    ! Determine equilibrium temperature in MeV
    kT = nuc % kT

    ! Check if energy is above threshold
    if (E >= FREE_GAS_THRESHOLD * kT .and. nuc % awr > ONE) then
      v_target = ZERO
      return
    end if

    ! calculate beta
    beta_vn = sqrt(nuc%awr * E / kT)

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

    ! determine speed of target nucleus
    vt = sqrt(beta_vt_sq*kT/nuc % awr)

    ! determine velocity vector of target nucleus based on neutron's velocity
    ! and the sampled angle between them
    v_target = vt * rotate_angle(uvw, mu)

  end subroutine sample_target_velocity

!===============================================================================
! CREATE_FISSION_SITES determines the average total, prompt, and delayed
! neutrons produced from fission and creates appropriate bank sites.
!===============================================================================

  subroutine create_fission_sites(p, i_nuclide, i_reaction)

    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: i_nuclide
    integer,        intent(in)    :: i_reaction

    integer :: i            ! loop index
    integer :: nu           ! actual number of neutrons produced
    integer :: ijk(3)       ! indices in ufs mesh
    real(8) :: nu_t         ! total nu
    real(8) :: mu           ! fission neutron angular cosine
    real(8) :: phi          ! fission neutron azimuthal angle
    real(8) :: weight       ! weight adjustment for ufs method
    logical :: in_mesh      ! source site in ufs mesh?
    type(Nuclide),  pointer, save :: nuc => null()
    type(Reaction), pointer, save :: rxn => null()
!$omp threadprivate(nuc, rxn)

    ! Get pointers
    nuc => nuclides(i_nuclide)
    rxn => nuc % reactions(i_reaction)

    ! TODO: Heat generation from fission

    ! If uniform fission source weighting is turned on, we increase of decrease
    ! the expected number of fission sites produced

    if (ufs) then
      ! Determine indices on ufs mesh for current location
      call get_mesh_indices(ufs_mesh, p % coord0 % xyz, ijk, in_mesh)
      if (.not. in_mesh) then
        call write_particle_restart(p)
        message = "Source site outside UFS mesh!"
        call fatal_error()
      end if

      if (source_frac(1,ijk(1),ijk(2),ijk(3)) /= ZERO) then
        weight = ufs_mesh % volume_frac / source_frac(1,ijk(1),ijk(2),ijk(3))
      else
        weight = ONE
      end if
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

    ! Bank source neutrons
    if (nu == 0 .or. n_bank == 3*work) return
    p % fission = .true. ! Fission neutrons will be banked
    do i = int(n_bank,4) + 1, int(min(n_bank + nu, 3*work),4)
      ! Bank source neutrons by copying particle data
      fission_bank(i) % xyz = p % coord0 % xyz

      ! Set weight of fission bank site
      fission_bank(i) % wgt = ONE/weight

      ! Sample cosine of angle -- fission neutrons are always emitted
      ! isotropically. Sometimes in ACE data, fission reactions actually have
      ! an angular distribution listed, but for those that do, it's simply just
      ! a uniform distribution in mu
      mu = TWO * prn() - ONE

      ! Sample azimuthal angle uniformly in [0,2*pi)
      phi = TWO*PI*prn()
      fission_bank(i) % uvw(1) = mu
      fission_bank(i) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
      fission_bank(i) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

      ! Sample secondary energy distribution for fission reaction and set energy
      ! in fission bank
      fission_bank(i) % E = sample_fission_energy(nuc, rxn, p % E)
    end do

    ! increment number of bank sites
    n_bank = min(n_bank + nu, 3*work)

    ! Store total weight banked for analog fission tallies
    p % n_bank   = nu
    p % wgt_bank = nu/weight

  end subroutine create_fission_sites

!===============================================================================
! SAMPLE_FISSION_ENERGY
!===============================================================================

  function sample_fission_energy(nuc, rxn, E) result(E_out)

    type(Nuclide),  pointer :: nuc
    type(Reaction), pointer :: rxn
    real(8), intent(in)     :: E     ! incoming energy of neutron
    real(8)                 :: E_out ! outgoing energy of fission neutron

    integer :: j            ! index on nu energy grid / precursor group
    integer :: lc           ! index before start of energies/nu values
    integer :: NR           ! number of interpolation regions
    integer :: NE           ! number of energies tabulated
    integer :: law          ! energy distribution law
    integer :: n_sample     ! number of times resampling
    real(8) :: nu_t         ! total nu
    real(8) :: nu_d         ! delayed nu
    real(8) :: mu           ! fission neutron angular cosine
    real(8) :: beta         ! delayed neutron fraction
    real(8) :: xi           ! random number
    real(8) :: yield        ! delayed neutron precursor yield
    real(8) :: prob         ! cumulative probability
    type(DistEnergy), pointer, save :: edist => null()
!$omp threadprivate(edist)

    ! Determine total nu
    nu_t = nu_total(nuc, E)

    ! Determine delayed nu
    nu_d = nu_delayed(nuc, E)

    ! Determine delayed neutron fraction
    beta = nu_d / nu_t

    if (prn() < beta) then
      ! ====================================================================
      ! DELAYED NEUTRON SAMPLED

      ! sampled delayed precursor group
      xi = prn()
      lc = 1
      prob = ZERO
      do j = 1, nuc % n_precursor
        ! determine number of interpolation regions and energies
        NR = int(nuc % nu_d_precursor_data(lc + 1))
        NE = int(nuc % nu_d_precursor_data(lc + 2 + 2*NR))

        ! determine delayed neutron precursor yield for group j
        yield = interpolate_tab1(nuc % nu_d_precursor_data( &
             lc+1:lc+2+2*NR+2*NE), E)

        ! Check if this group is sampled
        prob = prob + yield
        if (xi < prob) exit

        ! advance pointer
        lc = lc + 2 + 2*NR + 2*NE + 1
      end do

      ! if the sum of the probabilities is slightly less than one and the
      ! random number is greater, j will be greater than nuc %
      ! n_precursor -- check for this condition
      j = min(j, nuc % n_precursor)

      ! select energy distribution for group j
      law = nuc % nu_d_edist(j) % law
      edist => nuc % nu_d_edist(j)

      ! sample from energy distribution
      n_sample = 0
      do
        if (law == 44 .or. law == 61) then
          call sample_energy(edist, E, E_out, mu)
        else
          call sample_energy(edist, E, E_out)
        end if

        ! resample if energy is >= 20 MeV
        if (E_out < 20) exit

        ! check for large number of resamples
        n_sample = n_sample + 1
        if (n_sample == MAX_SAMPLE) then
          ! call write_particle_restart(p)
          message = "Resampled energy distribution maximum number of " // &
               "times for nuclide " // nuc % name
          call fatal_error()
        end if
      end do

    else
      ! ====================================================================
      ! PROMPT NEUTRON SAMPLED

      ! sample from prompt neutron energy distribution
      law = rxn % edist % law
      n_sample = 0
      do
        if (law == 44 .or. law == 61) then
          call sample_energy(rxn%edist, E, E_out, prob)
        else
          call sample_energy(rxn%edist, E, E_out)
        end if

        ! resample if energy is >= 20 MeV
        if (E_out < 20) exit

        ! check for large number of resamples
        n_sample = n_sample + 1
        if (n_sample == MAX_SAMPLE) then
          ! call write_particle_restart(p)
          message = "Resampled energy distribution maximum number of " // &
               "times for nuclide " // nuc % name
          call fatal_error()
        end if
      end do

    end if

  end function sample_fission_energy

!===============================================================================
! INELASTIC_SCATTER handles all reactions with a single secondary neutron (other
! than fission), i.e. level scattering, (n,np), (n,na), etc.
!===============================================================================

  subroutine inelastic_scatter(nuc, rxn, E, uvw, mu, wgt)

    type(Nuclide),  pointer :: nuc
    type(Reaction), pointer :: rxn
    real(8), intent(inout)  :: E      ! energy in lab (incoming/outgoing)
    real(8), intent(inout)  :: uvw(3) ! directional cosines
    real(8), intent(out)    :: mu     ! cosine of scattering angle in lab
    real(8), intent(inout)  :: wgt    ! particle weight

    integer :: law         ! secondary energy distribution law
    real(8) :: A           ! atomic weight ratio of nuclide
    real(8) :: E_in        ! incoming energy
    real(8) :: E_cm        ! outgoing energy in center-of-mass
    real(8) :: Q           ! Q-value of reaction

    ! copy energy of neutron
    E_in = E

    ! determine A and Q
    A = nuc % awr
    Q = rxn % Q_value

    ! determine secondary energy distribution law
    law = rxn % edist % law

    ! sample scattering angle
    mu = sample_angle(rxn, E_in)

    ! sample outgoing energy
    if (law == 44 .or. law == 61) then
      call sample_energy(rxn%edist, E_in, E, mu)
    elseif (law == 66) then
      call sample_energy(rxn%edist, E_in, E, A=A, Q=Q)
    else
      call sample_energy(rxn%edist, E_in, E)
    end if

    ! if scattering system is in center-of-mass, transfer cosine of scattering
    ! angle and outgoing energy from CM to LAB
    if (rxn % scatter_in_cm) then
      E_cm = E

      ! determine outgoing energy in lab
      E = E_cm + (E_in + TWO * mu * (A+ONE) * sqrt(E_in * E_cm)) &
           / ((A+ONE)*(A+ONE))

      ! determine outgoing angle in lab
      mu = mu * sqrt(E_cm/E) + ONE/(A+ONE) * sqrt(E_in/E)
    end if

    ! change direction of particle
    uvw = rotate_angle(uvw, mu)

    ! change weight of particle based on multiplicity
    wgt = rxn % multiplicity * wgt

  end subroutine inelastic_scatter

!===============================================================================
! SAMPLE_ANGLE samples the cosine of the angle between incident and exiting
! particle directions either from 32 equiprobable bins or from a tabular
! distribution.
!===============================================================================

  function sample_angle(rxn, E) result(mu)

    type(Reaction), pointer    :: rxn ! reaction
    real(8),        intent(in) :: E   ! incoming energy

    real(8)        :: xi      ! random number on [0,1)
    integer        :: interp  ! type of interpolation
    integer        :: type    ! angular distribution type
    integer        :: i       ! incoming energy bin
    integer        :: n       ! number of incoming energy bins
    integer        :: lc      ! location in data array
    integer        :: NP      ! number of points in cos distribution
    integer        :: k       ! index on cosine grid
    real(8)        :: r       ! interpolation factor on incoming energy
    real(8)        :: frac    ! interpolation fraction on cosine
    real(8)        :: mu0     ! cosine in bin k
    real(8)        :: mu1     ! cosine in bin k+1
    real(8)        :: mu      ! final cosine sampled
    real(8)        :: c_k     ! cumulative frequency at k
    real(8)        :: c_k1    ! cumulative frequency at k+1
    real(8)        :: p0,p1   ! probability distribution

    ! check if reaction has angular distribution -- if not, sample outgoing
    ! angle isotropically
    if (.not. rxn % has_angle_dist) then
      mu = TWO * prn() - ONE
      return
    end if

    ! determine number of incoming energies
    n = rxn % adist % n_energy

    ! find energy bin and calculate interpolation factor -- if the energy is
    ! outside the range of the tabulated energies, choose the first or last bins
    if (E < rxn % adist % energy(1)) then
      i = 1
      r = ZERO
    elseif (E > rxn % adist % energy(n)) then
      i = n - 1
      r = ONE
    else
      i = binary_search(rxn % adist % energy, n, E)
      r = (E - rxn % adist % energy(i)) / &
           (rxn % adist % energy(i+1) - rxn % adist % energy(i))
    end if

    ! Sample between the ith and (i+1)th bin
    if (r > prn()) i = i + 1

    ! check whether this is a 32-equiprobable bin or a tabular distribution
    lc  = rxn % adist % location(i)
    type = rxn % adist % type(i)
    if (type == ANGLE_ISOTROPIC) then
      mu = TWO * prn() - ONE
    elseif (type == ANGLE_32_EQUI) then
      ! sample cosine bin
      xi = prn()
      k = 1 + int(32.0_8*xi)

      ! calculate cosine
      mu0 = rxn % adist % data(lc + k)
      mu1 = rxn % adist % data(lc + k+1)
      mu = mu0 + (32.0_8 * xi - k) * (mu1 - mu0)

    elseif (type == ANGLE_TABULAR) then
      interp = int(rxn % adist % data(lc + 1))
      NP     = int(rxn % adist % data(lc + 2))

      ! determine outgoing cosine bin
      xi = prn()
      lc = lc + 2
      c_k = rxn % adist % data(lc + 2*NP + 1)
      do k = 1, NP - 1
        c_k1 = rxn % adist % data(lc + 2*NP + k+1)
        if (xi < c_k1) exit
        c_k = c_k1
      end do

      ! check to make sure k is <= NP - 1
      k = min(k, NP - 1)

      p0  = rxn % adist % data(lc + NP + k)
      mu0 = rxn % adist % data(lc + k)
      if (interp == HISTOGRAM) then
        ! Histogram interpolation
        if (p0 > ZERO) then
          mu = mu0 + (xi - c_k)/p0
        else
          mu = mu0
        end if

      elseif (interp == LINEAR_LINEAR) then
        ! Linear-linear interpolation
        p1  = rxn % adist % data(lc + NP + k+1)
        mu1 = rxn % adist % data(lc + k+1)

        frac = (p1 - p0)/(mu1 - mu0)
        if (frac == ZERO) then
          mu = mu0 + (xi - c_k)/p0
        else
          mu = mu0 + (sqrt(max(ZERO, p0*p0 + 2*frac*(xi - c_k))) - p0)/frac
        end if
      else
        ! call write_particle_restart(p)
        message = "Unknown interpolation type: " // trim(to_str(interp))
        call fatal_error()
      end if

      ! Because of floating-point roundoff, it may be possible for mu to be
      ! outside of the range [-1,1). In these cases, we just set mu to exactly
      ! -1 or 1

      if (abs(mu) > ONE) mu = sign(ONE,mu)

    else
      ! call write_particle_restart(p)
      message = "Unknown angular distribution type: " // trim(to_str(type))
      call fatal_error()
    end if

  end function sample_angle

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  function rotate_angle(uvw0, mu) result(uvw)

    real(8), intent(in) :: uvw0(3) ! directional cosine
    real(8), intent(in) :: mu      ! cosine of angle in lab or CM
    real(8)             :: uvw(3)  ! rotated directional cosine

    real(8) :: phi    ! azimuthal angle
    real(8) :: sinphi ! sine of azimuthal angle
    real(8) :: cosphi ! cosine of azimuthal angle
    real(8) :: a      ! sqrt(1 - mu^2)
    real(8) :: b      ! sqrt(1 - w^2)
    real(8) :: u0     ! original cosine in x direction
    real(8) :: v0     ! original cosine in y direction
    real(8) :: w0     ! original cosine in z direction

    ! Copy original directional cosines
    u0 = uvw0(1)
    v0 = uvw0(2)
    w0 = uvw0(3)

    ! Sample azimuthal angle in [0,2pi)
    phi = TWO * PI * prn()

    ! Precompute factors to save flops
    sinphi = sin(phi)
    cosphi = cos(phi)
    a = sqrt(max(ZERO, ONE - mu*mu))
    b = sqrt(max(ZERO, ONE - w0*w0))

    ! Need to treat special case where sqrt(1 - w**2) is close to zero by
    ! expanding about the v component rather than the w component
    if (b > 1e-10) then
      uvw(1) = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
      uvw(2) = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
      uvw(3) = mu*w0 - a*b*cosphi
    else
      b = sqrt(ONE - v0*v0)
      uvw(1) = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
      uvw(2) = mu*v0 - a*b*cosphi
      uvw(3) = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end function rotate_angle

!===============================================================================
! SAMPLE_ENERGY samples an outgoing energy distribution, either for a secondary
! neutron from a collision or for a prompt/delayed fission neutron
!===============================================================================

  recursive subroutine sample_energy(edist, E_in, E_out, mu_out, A, Q)

    type(DistEnergy),  pointer       :: edist
    real(8), intent(in)              :: E_in   ! incoming energy of neutron
    real(8), intent(out)             :: E_out  ! outgoing energy
    real(8), intent(inout), optional :: mu_out ! outgoing cosine of angle
    real(8), intent(in),    optional :: A      ! mass number of nuclide
    real(8), intent(in),    optional :: Q      ! Q-value of reaction

    integer :: i           ! index on incoming energy grid
    integer :: k           ! sampled index on outgoing grid
    integer :: l           ! sampled index on incoming grid
    integer :: n_sample    ! number of rejections
    integer :: lc          ! dummy index
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies
    integer :: NET         ! number of outgoing energies
    integer :: INTTp       ! combination of INTT and ND
    integer :: INTT        ! 1 = histogram, 2 = linear-linear
    integer :: JJ          ! 1 = histogram, 2 = linear-linear
    integer :: ND          ! number of discrete lines
    integer :: NP          ! number of points in distribution

    real(8) :: p_valid     ! probability of law validity

    real(8) :: E_i_1, E_i_K   ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_K ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_K       ! endpoints interpolated between i and i+1

    real(8) :: E_l_k, E_l_k1  ! adjacent E on outgoing grid l
    real(8) :: p_l_k, p_l_k1  ! adjacent p on outgoing grid l
    real(8) :: c_k, c_k1      ! cumulative probability

    real(8) :: KM_A           ! Kalbach-Mann parameter R
    real(8) :: KM_R           ! Kalbach-Mann parameter R
    real(8) :: A_k, A_k1      ! Kalbach-Mann A on outgoing grid l
    real(8) :: R_k, R_k1      ! Kalbach-Mann R on outgoing grid l

    real(8) :: Watt_a, Watt_b ! Watt spectrum parameters

    real(8) :: mu_k    ! angular cosine in bin k
    real(8) :: mu_k1   ! angular cosine in bin k+1
    real(8) :: p_k     ! angular pdf in bin k
    real(8) :: p_k1    ! angular pdf in bin k+1

    real(8) :: r           ! interpolation factor on incoming energy
    real(8) :: frac        ! interpolation factor on outgoing energy
    real(8) :: U           ! restriction energy
    real(8) :: T           ! nuclear temperature

    real(8) :: Ap          ! total mass ratio for n-body dist
    integer :: n_bodies    ! number of bodies for n-body dist
    real(8) :: E_max       ! parameter for n-body dist
    real(8) :: x, y, v     ! intermediate variables for n-body dist
    real(8) :: r1, r2, r3, r4, r5, r6

    ! ==========================================================================
    ! SAMPLE ENERGY DISTRIBUTION IF THERE ARE MULTIPLE

    if (associated(edist % next)) then
      p_valid = interpolate_tab1(edist % p_valid, E_in)

      if (prn() > p_valid) then
        if (edist % law == 44 .or. edist % law == 61) then
          call sample_energy(edist%next, E_in, E_out, mu_out)
        elseif (edist % law == 66) then
          call sample_energy(edist%next, E_in, E_out, A=A, Q=Q)
        else
          call sample_energy(edist%next, E_in, E_out)
        end if
        return
      end if
    end if

    ! Determine which secondary energy distribution law to use
    select case (edist % law)
    case (1)
      ! =======================================================================
      ! TABULAR EQUIPROBABLE ENERGY BINS

      ! read number of interpolation regions, incoming energies, and outgoing
      ! energies
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))
      NET = int(edist % data(3 + 2*NR + NE))
      if (NR > 0) then
        ! call write_particle_restart(p)
        message = "Multiple interpolation regions not supported while &
             &attempting to sample equiprobable energy bins."
        call fatal_error()
      end if

      ! determine index on incoming energy grid and interpolation factor
      lc = 2 + 2*NR
      i = binary_search(edist % data(lc+1:lc+NE), NE, E_in)
      r = (E_in - edist%data(lc+i)) / &
           (edist%data(lc+i+1) - edist%data(lc+i))

      ! Sample outgoing energy bin
      r1 = prn()
      k = 1 + int(NET * r1)

      ! Determine E_1 and E_K
      lc    = 3 + 3*NR + NE + (i-1)*NET
      E_i_1 = edist % data(lc + 1)
      E_i_K = edist % data(lc + NET)

      lc     = 3 + 3*NR + NE + i*NET
      E_i1_1 = edist % data(lc + 1)
      E_i1_K = edist % data(lc + NET)

      E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
      E_K = E_i_K + r*(E_i1_K - E_i_K)

      ! Randomly select between the outgoing table for incoming energy E_i and
      ! E_(i+1)
      if (prn() < r) then
        l = i + 1
      else
        l = i
      end if

      ! Determine E_l_k and E_l_k+1
      lc     = 3 + 2*NR + NE + (l-1)*NET
      E_l_k  = edist % data(lc+k)
      E_l_k1 = edist % data(lc+k+1)

      ! Determine E' (denoted here as E_out)
      r2 = prn()
      E_out  = E_l_k + r2*(E_l_k1 - E_l_k)

      ! Now interpolate between incident energy bins i and i + 1
      if (l == i) then
        E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
      else
        E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
      end if

    case (3)
      ! =======================================================================
      ! INELASTIC LEVEL SCATTERING

      E_out = edist%data(2) * (E_in - edist%data(1))

    case (4)
      ! =======================================================================
      ! CONTINUOUS TABULAR DISTRIBUTION

      ! read number of interpolation regions and incoming energies
      NR  = int(edist % data(1))
      NE  = int(edist % data(2 + 2*NR))
      if (NR == 1) then
        message = "Assuming linear-linear interpolation when sampling &
             &continuous tabular distribution"
        call warning()
      else if (NR > 1) then
        ! call write_particle_restart(p)
        message = "Multiple interpolation regions not supported while &
             &attempting to sample continuous tabular distribution."
        call fatal_error()
      end if

      ! find energy bin and calculate interpolation factor -- if the energy is
      ! outside the range of the tabulated energies, choose the first or last
      ! bins
      lc = 2 + 2*NR
      if (E_in < edist % data(lc+1)) then
        i = 1
        r = ZERO
      elseif (E_in > edist % data(lc+NE)) then
        i = NE - 1
        r = ONE
      else
        i = binary_search(edist % data(lc+1:lc+NE), NE, E_in)
        r = (E_in - edist%data(lc+i)) / &
             (edist%data(lc+i+1) - edist%data(lc+i))
      end if

      ! Sample between the ith and (i+1)th bin
      r2 = prn()
      if (r > r2) then
        l = i + 1
      else
        l = i
      end if

      ! interpolation for energy E1 and EK
      lc    = int(edist%data(2 + 2*NR + NE + i))
      NP    = int(edist%data(lc + 2))
      E_i_1 = edist%data(lc + 2 + 1)
      E_i_K = edist%data(lc + 2 + NP)

      lc     = int(edist%data(2 + 2*NR + NE + i + 1))
      NP     = int(edist%data(lc + 2))
      E_i1_1 = edist%data(lc + 2 + 1)
      E_i1_K = edist%data(lc + 2 + NP)

      E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
      E_K = E_i_K + r*(E_i1_K - E_i_K)

      ! determine location of outgoing energies, pdf, cdf for E(l)
      lc = int(edist % data(2 + 2*NR + NE + l))

      ! determine type of interpolation and number of discrete lines
      INTTp = int(edist % data(lc + 1))
      NP    = int(edist % data(lc + 2))
      if (INTTp > 10) then
        INTT = mod(INTTp,10)
        ND = (INTTp - INTT)/10
      else
        INTT = INTTp
        ND = 0
      end if

      if (ND > 0) then
        ! discrete lines present
        ! call write_particle_restart(p)
        message = "Discrete lines in continuous tabular distributed not &
             &yet supported"
        call fatal_error()
      end if

      ! determine outgoing energy bin
      r1 = prn()
      lc = lc + 2 ! start of EOUT
      c_k = edist % data(lc + 2*NP + 1)
      do k = 1, NP - 1
        c_k1 = edist % data(lc + 2*NP + k+1)
        if (r1 < c_k1) exit
        c_k = c_k1
      end do

      ! check to make sure k is <= NP - 1
      k = min(k, NP - 1)

      E_l_k = edist % data(lc+k)
      p_l_k = edist % data(lc+NP+k)
      if (INTT == HISTOGRAM) then
        ! Histogram interpolation
        if (p_l_k > ZERO) then
          E_out = E_l_k + (r1 - c_k)/p_l_k
        else
          E_out = E_l_k
        end if

      elseif (INTT == LINEAR_LINEAR) then
        ! Linear-linear interpolation
        E_l_k1 = edist % data(lc+k+1)
        p_l_k1 = edist % data(lc+NP+k+1)

        frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
        if (frac == ZERO) then
          E_out = E_l_k + (r1 - c_k)/p_l_k
        else
          E_out = E_l_k + (sqrt(max(ZERO, p_l_k*p_l_k + &
               2*frac*(r1 - c_k))) - p_l_k)/frac
        end if
      else
        ! call write_particle_restart(p)
        message = "Unknown interpolation type: " // trim(to_str(INTT))
        call fatal_error()
      end if

      ! Now interpolate between incident energy bins i and i + 1
      if (l == i) then
        E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
      else
        E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
      end if

    case (5)
      ! =======================================================================
      ! GENERAL EVAPORATION SPECTRUM

    case (7)
      ! =======================================================================
      ! MAXWELL FISSION SPECTRUM

      ! read number of interpolation regions and incoming energies
      NR = int(edist % data(1))
      NE = int(edist % data(2 + 2*NR))

      ! determine nuclear temperature from tabulated function
      T = interpolate_tab1(edist % data, E_in)

      ! determine restriction energy
      lc = 2 + 2*NR + 2*NE
      U = edist % data(lc + 1)

      n_sample = 0
      do
        ! sample maxwell fission spectrum
        E_out = maxwell_spectrum(T)

        ! accept energy based on restriction energy
        if (E_out <= E_in - U) exit

        ! check for large number of rejections
        n_sample = n_sample + 1
        if (n_sample == MAX_SAMPLE) then
          ! call write_particle_restart(p)
          message = "Too many rejections on Maxwell fission spectrum."
          call fatal_error()
        end if
      end do

    case (9)
      ! =======================================================================
      ! EVAPORATION SPECTRUM

      ! read number of interpolation regions and incoming energies
      NR = int(edist % data(1))
      NE = int(edist % data(2 + 2*NR))

      ! determine nuclear temperature from tabulated function
      T = interpolate_tab1(edist % data, E_in)

      ! determine restriction energy
      lc = 2 + 2*NR + 2*NE
      U = edist % data(lc + 1)

      ! sample outgoing energy based on evaporation spectrum probability
      ! density function
      n_sample = 0
      do
        r1 = prn()
        r2 = prn()
        E_out = -T * log(r1*r2)
        if (E_out <= E_in - U) exit

        ! check for large number of rejections
        n_sample = n_sample + 1
        if (n_sample == MAX_SAMPLE) then
          ! call write_particle_restart(p)
          message = "Too many rejections on evaporation spectrum."
          call fatal_error()
        end if
      end do

    case (11)
      ! =======================================================================
      ! ENERGY-DEPENDENT WATT SPECTRUM

      ! read number of interpolation regions and incoming energies for
      ! parameter 'a'
      NR = int(edist % data(1))
      NE = int(edist % data(2 + 2*NR))

      ! determine Watt parameter 'a' from tabulated function
      Watt_a = interpolate_tab1(edist % data, E_in)

      ! determine Watt parameter 'b' from tabulated function
      lc = 2 + 2*(NR + NE)
      Watt_b = interpolate_tab1(edist % data, E_in, lc + 1)

      ! read number of interpolation regions and incoming energies for
      ! parameter 'a'
      NR = int(edist % data(lc + 1))
      NE = int(edist % data(lc + 2 + 2*NR))

      ! determine restriction energy
      lc = lc + 2 + 2*(NR + NE)
      U = edist % data(lc + 1)

      n_sample = 0
      do
        ! Sample energy-dependent Watt fission spectrum
        E_out = watt_spectrum(Watt_a, Watt_b)

        ! accept energy based on restriction energy
        if (E_out <= E_in - U) exit

        ! check for large number of rejections
        n_sample = n_sample + 1
        if (n_sample == MAX_SAMPLE) then
          ! call write_particle_restart(p)
          message = "Too many rejections on Watt spectrum."
          call fatal_error()
        end if
      end do

    case (44)
      ! =======================================================================
      ! KALBACH-MANN CORRELATED SCATTERING

      if (.not. present(mu_out)) then
        ! call write_particle_restart(p)
        message = "Law 44 called without giving mu_out as argument."
        call fatal_error()
      end if

      ! read number of interpolation regions and incoming energies
      NR = int(edist % data(1))
      NE = int(edist % data(2 + 2*NR))
      if (NR > 0) then
        ! call write_particle_restart(p)
        message = "Multiple interpolation regions not supported while &
             &attempting to sample Kalbach-Mann distribution."
        call fatal_error()
      end if

      ! find energy bin and calculate interpolation factor -- if the energy is
      ! outside the range of the tabulated energies, choose the first or last
      ! bins
      lc = 2 + 2*NR
      if (E_in < edist % data(lc+1)) then
        i = 1
        r = ZERO
      elseif (E_in > edist % data(lc+NE)) then
        i = NE - 1
        r = ONE
      else
        i = binary_search(edist % data(lc+1:lc+NE), NE, E_in)
        r = (E_in - edist%data(lc+i)) / &
             (edist%data(lc+i+1) - edist%data(lc+i))
      end if

      ! Sample between the ith and (i+1)th bin
      r2 = prn()
      if (r > r2) then
        l = i + 1
      else
        l = i
      end if

      ! determine endpoints on grid i
      lc    = int(edist%data(2+2*NR+NE + i)) ! start of LDAT for i
      NP    = int(edist%data(lc + 2))
      E_i_1 = edist%data(lc + 2 + 1)
      E_i_K = edist%data(lc + 2 + NP)

      ! determine endpoints on grid i+1
      lc     = int(edist%data(2+2*NR+NE + i+1)) ! start of LDAT for i+1
      NP     = int(edist%data(lc + 2))
      E_i1_1 = edist%data(lc + 2 + 1)
      E_i1_K = edist%data(lc + 2 + NP)

      E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
      E_K = E_i_K + r*(E_i1_K - E_i_K)

      ! determine location of outgoing energies, pdf, cdf for E(l)
      lc = int(edist % data(2 + 2*NR + NE + l))

      ! determine type of interpolation and number of discrete lines
      INTTp = int(edist % data(lc + 1))
      NP    = int(edist % data(lc + 2))
      if (INTTp > 10) then
        INTT = mod(INTTp,10)
        ND = (INTTp - INTT)/10
      else
        INTT = INTTp
        ND = 0
      end if

      if (ND > 0) then
        ! discrete lines present
        ! call write_particle_restart(p)
        message = "Discrete lines in continuous tabular distributed not &
             &yet supported"
        call fatal_error()
      end if

      ! determine outgoing energy bin
      r1 = prn()
      lc = lc + 2 ! start of EOUT
      c_k = edist % data(lc + 2*NP + 1)
      do k = 1, NP - 1
        c_k1 = edist % data(lc + 2*NP + k+1)
        if (r1 < c_k1) exit
        c_k = c_k1
      end do

      ! check to make sure k is <= NP - 1
      k = min(k, NP - 1)

      E_l_k = edist % data(lc+k)
      p_l_k = edist % data(lc+NP+k)
      if (INTT == HISTOGRAM) then
        ! Histogram interpolation
        if (p_l_k > ZERO) then
          E_out = E_l_k + (r1 - c_k)/p_l_k
        else
          E_out = E_l_k
        end if

        ! Determine Kalbach-Mann parameters
        KM_R = edist % data(lc + 3*NP + k)
        KM_A = edist % data(lc + 4*NP + k)

      elseif (INTT == LINEAR_LINEAR) then
        ! Linear-linear interpolation
        E_l_k1 = edist % data(lc+k+1)
        p_l_k1 = edist % data(lc+NP+k+1)

        ! Find E prime
        frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
        if (frac == ZERO) then
          E_out = E_l_k + (r1 - c_k)/p_l_k
        else
          E_out = E_l_k + (sqrt(max(ZERO, p_l_k*p_l_k + &
               2*frac*(r1 - c_k))) - p_l_k)/frac
        end if

        ! Determine Kalbach-Mann parameters
        R_k  = edist % data(lc + 3*NP + k)
        R_k1 = edist % data(lc + 3*NP + k+1)
        A_k  = edist % data(lc + 4*NP + k)
        A_k1 = edist % data(lc + 4*NP + k+1)

        KM_R = R_k + (R_k1 - R_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
        KM_A = A_k + (A_k1 - A_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
      else
        ! call write_particle_restart()
        message = "Unknown interpolation type: " // trim(to_str(INTT))
        call fatal_error()
      end if

      ! Now interpolate between incident energy bins i and i + 1
      if (l == i) then
        E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
      else
        E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
      end if

      ! Sampled correlated angle from Kalbach-Mann parameters
      r3 = prn()
      r4 = prn()
      if (r3 > KM_R) then
        T = (TWO*r4 - ONE) * sinh(KM_A)
        mu_out = log(T + sqrt(T*T + ONE))/KM_A
      else
        mu_out = log(r4*exp(KM_A) + (ONE - r4)*exp(-KM_A))/KM_A
      end if

    case (61)
      ! =======================================================================
      ! CORRELATED ENERGY AND ANGLE DISTRIBUTION

      if (.not. present(mu_out)) then
        ! call write_particle_restart()
        message = "Law 61 called without giving mu_out as argument."
        call fatal_error()
      end if

      ! read number of interpolation regions and incoming energies
      NR = int(edist % data(1))
      NE = int(edist % data(2 + 2*NR))
      if (NR > 0) then
        ! call write_particle_restart()
        message = "Multiple interpolation regions not supported while &
             &attempting to sample correlated energy-angle distribution."
        call fatal_error()
      end if

      ! find energy bin and calculate interpolation factor -- if the energy is
      ! outside the range of the tabulated energies, choose the first or last
      ! bins
      lc = 2 + 2*NR
      if (E_in < edist % data(lc+1)) then
        i = 1
        r = ZERO
      elseif (E_in > edist % data(lc+NE)) then
        i = NE - 1
        r = ONE
      else
        i = binary_search(edist % data(lc+1:lc+NE), NE, E_in)
        r = (E_in - edist%data(lc+i)) / &
             (edist%data(lc+i+1) - edist%data(lc+i))
      end if

      ! Sample between the ith and (i+1)th bin
      r2 = prn()
      if (r > r2) then
        l = i + 1
      else
        l = i
      end if

      ! determine endpoints on grid i
      lc    = int(edist%data(2+2*NR+NE + i)) ! start of LDAT for i
      NP    = int(edist%data(lc + 2))
      E_i_1 = edist%data(lc + 2 + 1)
      E_i_K = edist%data(lc + 2 + NP)

      ! determine endpoints on grid i+1
      lc     = int(edist%data(2+2*NR+NE + i+1)) ! start of LDAT for i+1
      NP     = int(edist%data(lc + 2))
      E_i1_1 = edist%data(lc + 2 + 1)
      E_i1_K = edist%data(lc + 2 + NP)

      E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
      E_K = E_i_K + r*(E_i1_K - E_i_K)

      ! determine location of outgoing energies, pdf, cdf for E(l)
      lc = int(edist % data(2 + 2*NR + NE + l))

      ! determine type of interpolation and number of discrete lines
      INTTp = int(edist % data(lc + 1))
      NP    = int(edist % data(lc + 2))
      if (INTTp > 10) then
        INTT = mod(INTTp,10)
        ND = (INTTp - INTT)/10
      else
        INTT = INTTp
        ND = 0
      end if

      if (ND > 0) then
        ! discrete lines present
        ! call write_particle_restart()
        message = "Discrete lines in continuous tabular distributed not &
             &yet supported"
        call fatal_error()
      end if

      ! determine outgoing energy bin
      r1 = prn()
      lc = lc + 2 ! start of EOUT
      c_k = edist % data(lc + 2*NP + 1)
      do k = 1, NP - 1
        c_k1 = edist % data(lc + 2*NP + k+1)
        if (r1 < c_k1) exit
        c_k = c_k1
      end do

      ! check to make sure k is <= NP - 1
      k = min(k, NP - 1)

      E_l_k = edist % data(lc+k)
      p_l_k = edist % data(lc+NP+k)
      if (INTT == HISTOGRAM) then
        ! Histogram interpolation
        if (p_l_k > ZERO) then
          E_out = E_l_k + (r1 - c_k)/p_l_k
        else
          E_out = E_l_k
        end if

      elseif (INTT == LINEAR_LINEAR) then
        ! Linear-linear interpolation
        E_l_k1 = edist % data(lc+k+1)
        p_l_k1 = edist % data(lc+NP+k+1)

        ! Find E prime
        frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
        if (frac == ZERO) then
          E_out = E_l_k + (r1 - c_k)/p_l_k
        else
          E_out = E_l_k + (sqrt(max(ZERO, p_l_k*p_l_k + &
               2*frac*(r1 - c_k))) - p_l_k)/frac
        end if
      else
        ! call write_particle_restart()
        message = "Unknown interpolation type: " // trim(to_str(INTT))
        call fatal_error()
      end if

      ! Now interpolate between incident energy bins i and i + 1
      if (l == i) then
        E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
      else
        E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
      end if

      ! Find correlated angular distribution for closest outgoing energy bin
      if (r1 - c_k < c_k1 - r1) then
        lc = int(edist % data(lc + 3*NP + k))
      else
        lc = int(edist % data(lc + 3*NP + k + 1))
      end if

      ! Check if angular distribution is isotropic
      if (lc == 0) then
        mu_out = TWO * prn() - ONE
        return
      end if

      ! interpolation type and number of points in angular distribution
      JJ = int(edist % data(lc + 1))
      NP = int(edist % data(lc + 2))

      ! determine outgoing cosine bin
      r3 = prn()
      lc = lc + 2
      c_k = edist % data(lc + 2*NP + 1)
      do k = 1, NP - 1
        c_k1 = edist % data(lc + 2*NP + k+1)
        if (r3 < c_k1) exit
        c_k = c_k1
      end do

      ! check to make sure k is <= NP - 1
      k = min(k, NP - 1)

      p_k  = edist % data(lc + NP + k)
      mu_k = edist % data(lc + k)
      if (JJ == HISTOGRAM) then
        ! Histogram interpolation
        if (p_k > ZERO) then
          mu_out = mu_k + (r3 - c_k)/p_k
        else
          mu_out = mu_k
        end if

      elseif (JJ == LINEAR_LINEAR) then
        ! Linear-linear interpolation
        p_k1  = edist % data(lc + NP + k+1)
        mu_k1 = edist % data(lc + k+1)

        frac = (p_k1 - p_k)/(mu_k1 - mu_k)
        if (frac == ZERO) then
          mu_out = mu_k + (r3 - c_k)/p_k
        else
          mu_out = mu_k + (sqrt(p_k*p_k + 2*frac*(r3 - c_k))-p_k)/frac
        end if
      else
        ! call write_particle_restart()
        message = "Unknown interpolation type: " // trim(to_str(JJ))
        call fatal_error()
      end if

    case (66)
      ! =======================================================================
      ! N-BODY PHASE SPACE DISTRIBUTION

      ! read number of bodies in phase space and total mass ratio
      n_bodies = int(edist % data(1))
      Ap       = edist % data(2)

      ! determine E_max parameter
      E_max = (Ap - ONE)/Ap * (A/(A+ONE)*E_in + Q)

      ! x is essentially a Maxwellian distribution
      x = maxwell_spectrum(ONE)

      select case (n_bodies)
      case (3)
        y = maxwell_spectrum(ONE)
      case (4)
        r1 = prn()
        r2 = prn()
        r3 = prn()
        y = -log(r1*r2*r3)
      case (5)
        r1 = prn()
        r2 = prn()
        r3 = prn()
        r4 = prn()
        r5 = prn()
        r6 = prn()
        y = -log(r1*r2*r3*r4) - log(r5) * cos(PI/TWO*r6)**2
      end select

      ! now determine v and E_out
      v = x/(x+y)
      E_out = E_max * v

    case (67)
      ! =======================================================================
      ! LABORATORY ENERGY-ANGLE LAW

    end select

  end subroutine sample_energy

end module physics
