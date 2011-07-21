module physics

  use global
  use geometry,    only: find_cell, dist_to_boundary, cross_surface, &
       &                 cross_lattice
  use types,       only: Particle
  use mcnp_random, only: rang
  use output,      only: error, warning, message, print_particle
  use search,      only: binary_search
  use endf,        only: reaction_name
  use score,       only: score_tally

  implicit none

contains

!=====================================================================
! TRANSPORT encompasses the main logic for moving a particle through
! geometry.
!=====================================================================

  subroutine transport(p)

    type(Particle), pointer :: p

    integer        :: surf           ! surface which particle is on
    integer        :: IE             ! index on energy grid
    real(8)        :: d_to_boundary  ! distance to nearest boundary
    real(8)        :: d_to_collision ! sampled distance to collision
    real(8)        :: distance       ! distance particle travels
    real(8)        :: Sigma          ! total cross-section
    real(8)        :: f              ! interpolation factor
    logical        :: found_cell     ! found cell which particle is in?
    logical        :: in_lattice     ! is surface crossing in lattice?
    character(250) :: msg            ! output/error message
    type(Universe), pointer :: univ

    if (p % cell == 0) then
       univ => universes(BASE_UNIVERSE)
       call find_cell(univ, p, found_cell)

       ! if particle couldn't be located, print error
       if (.not. found_cell) then
          write(msg, 100) "Could not locate cell for particle at: ", p % xyz
100       format (A,3ES11.3)
          call error(msg)
       end if
    end if

    if (verbosity >= 9) then
       msg = "Simulating Particle " // trim(int_to_str(p % uid))
       call message(msg)
    end if

    if (verbosity >= 10) then
       msg = "    Born in cell " // trim(int_to_str(cells(p%cell)%uid))
       call message(msg)
    end if

    ! find energy index, interpolation factor
    call find_energy_index(p)

    do while (p % alive)

       ! Copy energy index and interpolation factor
       IE = p % IE
       f  = p % interp
    
       ! Determine material total cross-section
       Sigma = (1-f)*cMaterial%total_xs(IE) + f*cMaterial%total_xs(IE+1)

       ! Find the distance to the nearest boundary
       call dist_to_boundary(p, d_to_boundary, surf, in_lattice)

       ! Sample a distance to collision
       d_to_collision = -log(rang()) / Sigma
       
       ! Select smaller of the two distances
       distance = min(d_to_boundary, d_to_collision)

       ! Advance particle
       p%xyz = p%xyz + distance * p%uvw
       p%xyz_local = p%xyz_local + distance * p%uvw

       if (d_to_collision > d_to_boundary) then
          p % cell = 0
          if (in_lattice) then
             p % surface = 0
             call cross_lattice(p)
          else
             p % surface = surf
             call cross_surface(p)
          end if
       else
          ! collision
          p % surface = 0
          call collision(p)
       end if
       
    end do

  end subroutine transport

!=====================================================================
! FIND_ENERGY_INDEX determines the index on the union energy grid and
! the interpolation factor for a particle at a certain energy
!=====================================================================

  subroutine find_energy_index(p)

    type(Particle), pointer :: p

    integer :: IE     ! index on union energy grid
    real(8) :: E      ! energy of particle
    real(8) :: interp ! interpolation factor

    ! copy particle's energy
    E = p % E

    ! if particle's energy is outside of energy grid range, set to
    ! first or last index. Otherwise, do a binary search through the
    ! union energy grid.
    if (E < e_grid(1)) then
       IE = 1
    elseif (E > e_grid(n_grid)) then
       IE = n_grid - 1
    else
       IE = binary_search(e_grid, n_grid, E)
    end if
    
    ! calculate the interpolation factor -- note this will be outside
    ! of [0,1) for a particle outside the energy range of the union
    ! grid
    interp = (E - e_grid(IE))/(e_grid(IE+1) - e_grid(IE))

    ! set particle attributes
    p % IE     = IE
    p % interp = interp
    
  end subroutine find_energy_index

!=====================================================================
! COLLISION samples a nuclide and reaction and then calls the
! appropriate routine for that reaction
!=====================================================================

  subroutine collision(p)

    type(Particle), pointer :: p

    integer :: i          ! index over isotopes in a material
    integer :: n_isotopes ! number of isotopes in material
    integer :: IE         ! index on nuclide energy grid
    real(8) :: f          ! interpolation factor
    real(8) :: sigma      ! microscopic total xs for nuclide
    real(8) :: total      ! total macroscopic xs for material
    real(8) :: density    ! atom density of material
    real(8) :: density_i  ! atom density of nuclide
    real(8) :: prob       ! cumulative probability
    real(8) :: r1         ! random number
    real(8) :: flux       ! collision estimator of flux
    real(8), allocatable :: Sigma_rxn(:) ! macroscopic xs for each nuclide
    character(250) :: msg ! output/error message
    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn

    density = cMaterial%atom_density

    ! calculate total cross-section for each nuclide at current energy
    ! in order to create discrete pdf for sampling nuclide
    n_isotopes = cMaterial%n_isotopes
    allocate(Sigma_rxn(n_isotopes))
    do i = 1, n_isotopes
       table => xs_continuous(cMaterial%table(i))

       ! determine nuclide atom density
       density_i = cMaterial%atom_percent(i) * density

       ! search nuclide energy grid
       IE = table%grid_index(p % IE)
       f = (p%E - table%energy(IE))/(table%energy(IE+1) - table%energy(IE))

       ! calculate nuclide macroscopic cross-section
       Sigma_rxn(i) = density_i*((ONE-f)*table%Sigma_t(IE) + & 
            & f*(table%Sigma_t(IE+1)))
    end do
    total = sum(Sigma_rxn)

    ! sample nuclide
    r1 = rang()
    prob = ZERO
    do i = 1, n_isotopes
       prob = prob + Sigma_rxn(i) / total
       if (r1 < prob) exit
    end do

    ! Get table, total xs, interpolation factor
    density_i = cMaterial%atom_percent(i)*density
    table => xs_continuous(cMaterial%table(i))
    sigma = Sigma_rxn(i) / density_i
    IE = table%grid_index(p % IE)
    f = (p%E - table%energy(IE))/(table%energy(IE+1) - table%energy(IE))

    ! free memory
    deallocate(Sigma_rxn)

    ! sample reaction channel
    r1 = rang()*sigma
    prob = ZERO
    do i = 1, table%n_reaction
       rxn => table%reactions(i)

       ! some materials have gas production cross sections with MT >
       ! 200 that are duplicates. Also MT=4 is total level inelastic
       ! scattering which should be skipped
       if (rxn%MT >= 200 .or. rxn%MT == 4) cycle

       ! if energy is below threshold for this reaction, skip it
       if (IE < rxn%IE) cycle

       ! add to cumulative probability
       prob = prob + ((ONE-f)*rxn%sigma(IE-rxn%IE+1) & 
            & + f*(rxn%sigma(IE-rxn%IE+2)))
       if (r1 < prob) exit
    end do

    if (verbosity >= 10) then
       msg = "    " // trim(reaction_name(rxn%MT)) // " with nuclide " // &
            & trim(table%name)
       call message(msg)
    end if

    ! call appropriate subroutine
    select case (rxn % MT)
    case (ELASTIC)
       call elastic_scatter(p, table, rxn)
    case (N_2ND, N_2N, N_3N, N_2NA, N_3NA, N_2N2A, N_4N, N_2NP, N_3NP)
       call n_xn(p, table, rxn)
    case (N_NA, N_N3A, N_NP, N_N2A, N_ND, N_NT, N_N3HE, &
         & N_NT2A, N_N2P, N_NPA, N_N1 : N_NC)
       call inelastic_scatter(p, table, rxn)
    case (FISSION, N_F, N_NF, N_2NF, N_3NF)
       call n_fission(p, table, rxn)
    case (N_GAMMA : N_DA)
       call n_absorption(p)
    case default
       msg = "Cannot simulate reaction with MT " // int_to_str(rxn%MT)
       call warning(msg)
    end select

    ! Add collision estimator tallies
    if (tallies_on) then
       flux = p%wgt / SIGMA
       call score_tally(p, flux)
    end if

    ! check for very low energy
    if (p % E < 1.0e-100_8) then
       p % alive = .false.
       ! msg = "Killing neutron with extremely low energy"
       ! call warning(msg)
    end if

    ! find energy index, interpolation factor
    call find_energy_index(p)

  end subroutine collision

!=====================================================================
! ELASTIC_SCATTER treats the elastic scattering of a neutron with a
! target. Currently this assumes target-at-rest kinematics --
! obviously will need to be fixed
!=====================================================================

  subroutine elastic_scatter(p, table, rxn)

    type(Particle),      pointer :: p
    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn

    real(8) :: awr ! atomic weight ratio of target
    real(8) :: phi ! azimuthal angle
    real(8) :: mu  ! cosine of polar angle
    real(8) :: vx  ! velocity of neutron in x-direction
    real(8) :: vy  ! velocity of neutron in y-direction
    real(8) :: vz  ! velocity of neutron in z-direction
    real(8) :: vcx ! velocity of CM in x-direction
    real(8) :: vcy ! velocity of CM in y-direction
    real(8) :: vcz ! velocity of CM in z-direction
    real(8) :: vel ! magnitude of velocity
    real(8) :: u   ! x-direction
    real(8) :: v   ! y-direction
    real(8) :: w   ! z-direction
    real(8) :: E   ! energy

    vel = sqrt(p % E)
    awr = table % awr

    vx = vel * p%uvw(1)
    vy = vel * p%uvw(2)
    vz = vel * p%uvw(3)

    vcx = vx/(awr + ONE)
    vcy = vy/(awr + ONE)
    vcz = vz/(awr + ONE)

    ! Transform to CM frame
    vx = vx - vcx
    vy = vy - vcy
    vz = vz - vcz

    vel = sqrt(vx*vx + vy*vy + vz*vz)

    ! Sample scattering angle
    mu = sample_angle(rxn, p % E)

    ! Determine direction cosines in CM
    u = vx/vel
    v = vy/vel
    w = vz/vel

    ! Change direction cosines according to mu
    call rotate_angle(u, v, w, mu)

    vx = u*vel
    vy = v*vel
    vz = w*vel

    ! Transform back to LAB frame
    vx = vx + vcx
    vy = vy + vcy
    vz = vz + vcz

    E = vx*vx + vy*vy + vz*vz
    vel = sqrt(E)

    p % E = E
    p % uvw(1) = vx/vel
    p % uvw(2) = vy/vel
    p % uvw(3) = vz/vel

  end subroutine elastic_scatter

!=====================================================================
! N_FISSION determines the average total, prompt, and delayed neutrons
! produced from fission and creates appropriate bank sites. This
! routine will not work with implicit absorption, namely sampling of
! the number of neutrons!
!=====================================================================

  subroutine n_fission(p, table, rxn)

    type(Particle),      pointer :: p
    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn

    integer :: i           ! loop index
    integer :: j           ! index on nu energy grid / precursor group
    integer :: k           ! index on precursor yield grid
    integer :: loc         ! index before start of energies/nu values
    integer :: NC          ! number of polynomial coefficients
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies tabulated
    integer :: nu          ! actual number of neutrons produced
    integer :: law         ! energy distribution law
    real(8) :: E           ! incoming energy of neutron
    real(8) :: E_out       ! outgoing energy of fission neutron
    real(8) :: c           ! polynomial coefficient
    real(8) :: f           ! interpolation factor
    real(8) :: nu_total    ! total nu
    real(8) :: nu_prompt   ! prompt nu
    real(8) :: nu_delay    ! delayed nu
    real(8) :: mu          ! fission neutron angular cosine
    real(8) :: phi         ! fission neutron azimuthal angle
    real(8) :: beta        ! delayed neutron fraction
    real(8) :: xi          ! random number
    real(8) :: yield       ! delayed neutron precursor yield
    real(8) :: prob        ! cumulative probability
    character(250) :: msg  ! error message

    ! copy energy of neutron
    E = p % E

    ! ================================================================
    ! DETERMINE TOTAL NU
    if (table % nu_t_type == NU_NONE) then
       msg = "No neutron emission data for table: " // table % name
       call error(msg)
    elseif (table % nu_t_type == NU_POLYNOMIAL) then
       ! determine number of coefficients
       NC = int(table % nu_t_data(1))

       ! sum up polynomial in energy
       nu_total = ZERO
       do i = 0, NC - 1
          c = table % nu_t_data(i+2)
          nu_total = nu_total + c * E**i
       end do
    elseif (table % nu_t_type == NU_TABULAR) then
       ! determine number of interpolation regions -- as far as I can
       ! tell, no nu data has multiple interpolation
       ! regions. Furthermore, it seems all are lin-lin.
       NR = int(table % nu_t_data(1))
       if (NR /= 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to determine total nu."
          call error(msg)
       end if

       ! determine number of energies
       loc = 2 + 2*NR
       NE = int(table % nu_t_data(loc))

       ! do binary search over tabuled energies to determine
       ! appropriate index and interpolation factor
       if (E < table % nu_t_data(loc+1)) then
          j = 1
          f = ZERO
       elseif (E > table % nu_t_data(loc+NE)) then
          j = NE - 1
          f = ONE
       else
          j = binary_search(table % nu_t_data(loc+1), NE, E)
          f = (E - table % nu_t_data(loc+j)) / &
               & (table % nu_t_data(loc+j+1) - table % nu_t_data(loc+j))
       end if

       ! determine nu total
       loc = loc + NE
       nu_total = table % nu_t_data(loc+j) + f * & 
            & (table % nu_t_data(loc+j+1) - table % nu_t_data(loc+j))
    end if
          
    ! ================================================================
    ! DETERMINE PROMPT NU
    if (table % nu_p_type == NU_NONE) then
       ! since no prompt or delayed data is present, this means all
       ! neutron emission is prompt
       nu_prompt = nu_total
    elseif (table % nu_p_type == NU_POLYNOMIAL) then
       ! determine number of coefficients
       NC = int(table % nu_p_data(1))

       ! sum up polynomial in energy
       nu_prompt = ZERO
       do i = 0, NC - 1
          c = table % nu_p_data(i+2)
          nu_prompt = nu_prompt + c * E**i
       end do
    elseif (table % nu_p_type == NU_TABULAR) then
       ! determine number of interpolation regions
       NR = int(table % nu_p_data(1))
       if (NR /= 0) then
          msg = "Multiple interpolation regions not supported while & 
               &attempting to determine prompt nu."
          call error(msg)
       end if

       ! determine number of energies
       loc = 2 + 2*NR
       NE = int(table % nu_p_data(loc))

       ! do binary search over tabuled energies to determine
       ! appropriate index and interpolation factor
       if (E < table % nu_p_data(loc+1)) then
          j = 1
          f = ZERO
       elseif (E > table % nu_p_data(loc+NE)) then
          j = NE - 1
          f = ONE
       else
          j = binary_search(table % nu_p_data(loc+1), NE, E)
          f = (E - table % nu_p_data(loc+j)) / &
               & (table % nu_p_data(loc+j+1) - table % nu_p_data(loc+j))
       end if

       ! determine nu total
       loc = loc + NE
       nu_prompt = table % nu_p_data(loc+j) + f * & 
            & (table % nu_p_data(loc+j+1) - table % nu_p_data(loc+j))
    end if
       
    ! ================================================================
    ! DETERMINE DELAYED NU
    if (table % nu_d_type == NU_NONE) then
       nu_delay = ZERO
    elseif (table % nu_d_type == NU_TABULAR) then
       ! determine number of interpolation regions
       NR = int(table % nu_d_data(1))
       if (NR /= 0) then
          msg = "Multiple interpolation regions not supported while & 
               &attempting to determine delayed nu."
          call error(msg)
       end if

       ! determine number of energies
       loc = 2 + 2*NR
       NE = int(table % nu_d_data(loc))

       ! do binary search over tabuled energies to determine
       ! appropriate index and interpolation factor
       if (E < table % nu_d_data(loc+1)) then
          j = 1
          f = ZERO
       elseif (E > table % nu_d_data(loc+NE)) then
          j = NE - 1
          f = ONE
       else
          j = binary_search(table % nu_d_data(loc+1), NE, E)
          f = (E - table % nu_d_data(loc+j)) / &
               & (table % nu_d_data(loc+j+1) - table % nu_d_data(loc+j))
       end if

       ! determine nu total
       loc = loc + NE
       nu_delay = table % nu_d_data(loc+j) + f * & 
            & (table % nu_d_data(loc+j+1) - table % nu_d_data(loc+j))
    end if

    beta = nu_delay / nu_total

    ! TODO: Heat generation from fission

    ! Sample number of neutrons produced
    nu_total = p%wgt / keff * nu_total
    if (rang() > nu_total - int(nu_total)) then
       nu = int(nu_total)
    else
       nu = int(nu_total) + 1
    end if

    ! Bank source neutrons
    if (nu == 0 .or. n_bank == 3*n_particles) return
    do i = n_bank + 1, min(n_bank + nu, 3*n_particles)
       ! Bank source neutrons by copying particle data
       fission_bank(i) % uid = p % uid
       fission_bank(i) % xyz = p % xyz

       ! sample cosine of angle
       mu = sample_angle(rxn, E)

       ! sample between delayed and prompt neutrons
       if (rang() < beta) then
          ! ==========================================================
          ! DELAYED NEUTRON SAMPLED

          ! sampled delayed precursor group
          xi = rang()
          loc = 1
          prob = ZERO
          do j = 1, table % n_precursor
             ! determine number of interpolation regions and energies
             NR  = table % nu_d_precursor_data(loc + 1)
             NE  = table % nu_d_precursor_data(loc + 2 + 2*NR)
             if (NR > 0) then
                msg = "Multiple interpolation regions not supported while & 
                     &sampling delayed neutron precursor yield."
                call error(msg)
             end if

             ! interpolate on energy grid
             loc = loc + 2 + 2*NR
             if (E < table%nu_d_precursor_data(loc+1)) then
                k = 1
                f = ZERO
             elseif (E > table%nu_d_precursor_data(loc+NE)) then
                k = NE - 1
                f = ONE
             else
                k = binary_search(table%nu_d_precursor_data(loc+1), NE, E)
                f = (E - table%nu_d_precursor_data(loc+k)) / & 
                     & (table%nu_d_precursor_data(loc+k+1) - &
                     & table%nu_d_precursor_data(loc+k))
             end if

             ! determine delayed neutron precursor yield for group j
             loc = loc + NE
             yield = table%nu_d_precursor_data(loc+k) + f * &
                  (table%nu_d_precursor_data(loc+k+1) - &
                  & table%nu_d_precursor_data(loc+k))
             prob = prob + yield
             if (xi < prob) exit

             ! advance pointer
             loc = loc + NE + 1
          end do

          ! sample from energy distribution for group j
          law = table % nu_d_edist(j) % law
          do
             if (law == 44 .or. law == 61) then
                call sample_energy(table%nu_d_edist(j), E, E_out, mu)
             else
                call sample_energy(table%nu_d_edist(j), E, E_out)
             end if
             ! resample if energy is >= 20 MeV
             if (E_out < 20) exit
          end do

       else
          ! ==========================================================
          ! PROMPT NEUTRON SAMPLED

          ! sample from prompt neutron energy distribution
          law = rxn % edist % law
          do
             if (law == 44 .or. law == 61) then
                call sample_energy(rxn%edist, E, E_out, prob)
             else
                call sample_energy(rxn%edist, E, E_out)
             end if
             ! resample if energy is >= 20 MeV
             if (E_out < 20) exit
          end do

       end if

       ! Sample azimuthal angle uniformly in [0,2*pi)
       phi = TWO*PI*rang()
       fission_bank(i) % uvw(1) = mu
       fission_bank(i) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
       fission_bank(i) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

       ! set energy of fission neutron
       fission_bank(i) % E = E_out
    end do

    ! increment number of bank sites
    n_bank = min(n_bank + nu, 3*n_particles)

    ! kill original neutron
    p % alive = .false.

  end subroutine n_fission

!=====================================================================
! INELASTIC_SCATTER handles all reactions with a single secondary
! neutron (other than fission), i.e. level scattering, (n,np), (n,na),
! etc.
!=====================================================================

  subroutine inelastic_scatter(p, table, rxn)

    type(Particle),      pointer :: p
    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn

    integer :: law
    real(8) :: A          ! atomic weight ratio of nuclide
    real(8) :: E_in       ! incoming energy
    real(8) :: mu         ! cosine of scattering angle
    real(8) :: E          ! outgoing energy in laboratory
    real(8) :: E_cm       ! outgoing energy in center-of-mass
    real(8) :: u,v,w      ! direction cosines
    character(250) :: msg ! error message
    
    ! copy energy of neutron
    E_in = p % E

    ! determine A
    A = table % awr

    ! determine secondary energy distribution law
    law = rxn % edist % law

    ! sample scattering angle
    mu = sample_angle(rxn, E_in)

    ! sample outgoing energy
    if (law == 44 .or. law == 61) then
       call sample_energy(rxn%edist, E_in, E, mu)
    else
       call sample_energy(rxn%edist, E_in, E)
    end if

    ! if scattering system is in center-of-mass, transfer cosine of
    ! scattering angle and outgoing energy from CM to LAB
    if (rxn % TY < 0) then
       E_cm = E

       ! determine outgoing energy in lab
       E = E_cm + (E_in + TWO * mu * (A+ONE) * sqrt(E_in * E_cm)) & 
            & / ((A+ONE)*(A+ONE))

       ! determine outgoing angle in lab
       mu = mu * sqrt(E_cm/E) + ONE/(A+ONE) * sqrt(E_in/E)
    end if

    ! copy directional cosines
    u = p % uvw(1)
    v = p % uvw(2)
    w = p % uvw(3)

    ! change direction of particle
    call rotate_angle(u, v, w, mu)
    p % uvw = (/ u, v, w /)

    ! change energy of particle
    p % E = E

  end subroutine inelastic_scatter

!=====================================================================
! N_XN
!=====================================================================

  subroutine n_xn(p, table, rxn)

    type(Particle),      pointer :: p
    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn

    integer :: i           ! loop index
    integer :: n_secondary ! number of secondary particles
    integer :: law         ! secondary energy distribution law
    real(8) :: A           ! atomic weight ratio of nuclide
    real(8) :: E_in        ! incoming energy
    real(8) :: mu          ! cosine of scattering angle
    real(8) :: E           ! outgoing energy in laboratory
    real(8) :: E_cm        ! outgoing energy in center-of-mass
    real(8) :: u,v,w       ! direction cosines
    character(250) :: msg  ! error message
    
    ! copy energy of neutron
    E_in = p % E

    ! determine A
    A = table % awr
    
    ! determine secondary energy distribution law
    law = rxn % edist % law

    ! sample scattering angle
    mu = sample_angle(rxn, E_in)

    ! sample outgoing energy
    if (law == 44 .or. law == 61) then
       call sample_energy(rxn%edist, E_in, E, mu)
    else
       call sample_energy(rxn%edist, E_in, E)
    end if

    ! if scattering system is in center-of-mass, transfer cosine of
    ! scattering angle and outgoing energy from CM to LAB
    if (rxn % TY < 0) then
       E_cm = E

       ! determine outgoing energy in lab
       E = E_cm + (E_in + TWO * mu * (A+ONE) * sqrt(E_in * E_cm)) & 
            & / ((A+ONE)*(A+ONE))

       ! determine outgoing angle in lab
       mu = mu * sqrt(E_cm/E) + ONE/(A+ONE) * sqrt(E_in/E)
    end if

    ! copy directional cosines
    u = p % uvw(1)
    v = p % uvw(2)
    w = p % uvw(3)

    ! change direction of particle
    call rotate_angle(u, v, w, mu)
    p % uvw = (/ u, v, w /)

    ! change energy of particle
    p % E = E

    ! produce secondary particles in source bank
    n_secondary = abs(rxn % TY) - 1
    do i = 1, n_secondary
       ! TODO: create particle in source bank
    end do

  end subroutine n_xn

!=====================================================================
! N_ABSORPTION handles all absorbing reactions, i.e. (n,gamma), (n,p),
! (n,a), etc.
!=====================================================================

  subroutine n_absorption(p)

    type(Particle), pointer :: p

    integer        :: cell_num ! user-specified cell number
    character(250) :: msg      ! output/error message

    p % alive = .false.
    if (verbosity >= 10) then
       cell_num = cells(p % cell)%uid
       msg = "    Absorbed in cell " // trim(int_to_str(cell_num))
       call message(msg)
    end if

  end subroutine n_absorption

!=====================================================================
! SAMPLE_ANGLE samples the cosine of the angle between incident and
! exiting particle directions either from 32 equiprobable bins or from
! a tabular distribution.
!=====================================================================

  function sample_angle(rxn, E) result(mu)

    type(AceReaction), pointer    :: rxn ! reaction
    real(8),           intent(in) :: E   ! incoming energy

    real(8)        :: xi      ! random number on [0,1)
    integer        :: interp  ! type of interpolation
    integer        :: type    ! angular distribution type
    integer        :: i       ! incoming energy bin
    integer        :: n       ! number of incoming energy bins
    integer        :: loc     ! location in data array
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
    character(250) :: msg     ! error message

    ! check if reaction has angular distribution -- if not, sample
    ! outgoing angle isotropically
    if (.not. rxn % has_angle_dist) then
       mu = TWO * rang() - ONE
       return
    end if

    ! determine number of incoming energies
    n = rxn % adist % n_energy

    ! find energy bin and calculate interpolation factor -- if the
    ! energy is outside the range of the tabulated energies, choose
    ! the first or last bins
    if (E < rxn % adist % energy(1)) then
       i = 1
       r = ZERO
    elseif (E > rxn % adist % energy(n)) then
       i = n - 1
       r = ONE
    else
       i = binary_search(rxn % adist % energy, n, E)
       r = (E - rxn % adist % energy(i)) / & 
            & (rxn % adist % energy(i+1) - rxn % adist % energy(i))
    end if

    ! Sample between the ith and (i+1)th bin
    if (r > rang()) i = i + 1

    ! check whether this is a 32-equiprobable bin or a tabular
    ! distribution
    loc  = rxn % adist % location(i)
    type = rxn % adist % type(i)
    if (type == ANGLE_ISOTROPIC) then
       mu = TWO * rang() - ONE
    elseif (type == ANGLE_32_EQUI) then
       ! sample cosine bin
       xi = rang()
       k = 1 + int(32.0_8*xi)

       ! calculate cosine
       mu0 = rxn % adist % data(loc + k)
       mu1 = rxn % adist % data(loc + k+1)
       mu = mu0 + (32.0_8 * xi - k) * (mu1 - mu0)

    elseif (type == ANGLE_TABULAR) then
       interp = rxn % adist % data(loc + 1)
       NP     = rxn % adist % data(loc + 2)

       ! determine outgoing cosine bin
       xi = rang()
       loc = loc + 2
       c_k = rxn % adist % data(loc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = rxn % adist % data(loc + 2*NP + k+1)
          if (xi < c_k1) exit
          c_k = c_k1
       end do

       p0  = rxn % adist % data(loc + NP + k)
       mu0 = rxn % adist % data(loc + k)
       if (interp == HISTOGRAM) then
          ! Histogram interpolation
          mu = mu0 + (xi - c_k)/p0

       elseif (interp == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          p1  = rxn % adist % data(loc + NP + k+1)
          mu1 = rxn % adist % data(loc + k+1)

          frac = (p1 - p0)/(mu1 - mu0)
          if (frac == ZERO) then
             mu = mu0 + (xi - c_k)/p0
          else
             mu = mu0 + (sqrt(p0*p0 + 2*frac*(xi - c_k))-p0)/frac
          end if
       else
          msg = "Unknown interpolation type: " // trim(int_to_str(interp))
          call error(msg)
       end if

       if (abs(mu) > ONE) then
          msg = "Sampled cosine of angle outside [-1, 1)."
          call warning(msg)

          mu = sign(ONE,mu)
       end if
         
    else
       msg = "Unknown angular distribution type: " // trim(int_to_str(type))
       call error(msg)
    end if
    
  end function sample_angle

!=====================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose
! cosine is mu and through an azimuthal angle sampled uniformly. Note
! that this is done with direct sampling rather than rejection as is
! done in MCNP and SERPENT.
!=====================================================================

  subroutine rotate_angle(u, v, w, mu)

!    type(Particle), pointer :: p
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v
    real(8), intent(inout) :: w
    real(8), intent(in)    :: mu ! cosine of angle in lab

    real(8) :: phi, sinphi, cosphi
    real(8) :: a,b
    real(8) :: u0, v0, w0

    ! Copy original directional cosines
    u0 = u
    v0 = v
    w0 = w

    ! Sample azimuthal angle in [0,2pi)
    phi = TWO * PI * rang()

    ! Precompute factors to save flops
    sinphi = sin(phi)
    cosphi = cos(phi)
    a = sqrt(ONE - mu*mu)
    b = sqrt(ONE - w0*w0)

    ! Need to treat special case where sqrt(1 - w**2) is close to zero
    ! by expanding about the v component rather than the w component
    if (b > 1e-10) then
       u = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
       v = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
       w = mu*w0 - a*b*cosphi
    else
       b = sqrt(ONE - v0*v0)
       u = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
       v = mu*v0 - a*b*cosphi
       w = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end subroutine rotate_angle
    
!=====================================================================
! SAMPLE_ENERGY
!=====================================================================

  subroutine sample_energy(edist, E_in, E_out, mu_out)

    type(AceDistEnergy), intent(inout):: edist
    real(8), intent(in)               :: E_in
    real(8), intent(out)              :: E_out
    real(8), intent(inout), optional  :: mu_out

    integer :: i           ! index on incoming energy grid
    integer :: k           ! sampled index on outgoing grid
    integer :: l           ! sampled index on incoming grid
    integer :: loc         ! dummy index
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies
    integer :: NET         ! number of outgoing energies
    integer :: INTTp       ! combination of INTT and ND
    integer :: INTT        ! 1 = histogram, 2 = linear-linear
    integer :: JJ          ! 1 = histogram, 2 = linear-linear
    integer :: ND          ! number of discrete lines
    integer :: NP          ! number of points in distribution

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

    real(8) :: E_cm
    real(8) :: xi1, xi2, xi3, xi4
    real(8) :: r           ! interpolation factor on incoming energy
    real(8) :: frac        ! interpolation factor on outgoing energy
    real(8) :: U           ! restriction energy
    real(8) :: T           ! nuclear temperature
    character(250) :: msg  ! error message

    ! TODO: If there are multiple scattering laws, sample scattering
    ! law

    ! Check for multiple interpolation regions
    if (edist % n_interp > 0) then
       msg = "Multiple interpolation regions not supported while &
            &attempting to sample secondary energy distribution."
       call error(msg)
    end if
       
    ! Determine which secondary energy distribution law to use
    select case (edist % law)
    case (1)
       ! =============================================================
       ! TABULAR EQUIPROBABLE ENERGY BINS

       ! read number of interpolation regions, incoming energies, and
       ! outgoing energies
       NR  = edist % data(1)
       NE  = edist % data(2 + 2*NR)
       NET = edist % data(3 + 2*NR + NE)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample equiprobable energy bins."
          call error(msg)
       end if

       ! determine index on incoming energy grid and interpolation
       ! factor
       loc = 2 + 2*NR
       i = binary_search(edist % data(loc+1), NE, E_in)
       r = (E_in - edist%data(loc+i)) / &
            & (edist%data(loc+i+1) - edist%data(loc+i))

       ! Sample outgoing energy bin
       xi1 = rang()
       k = 1 + int(NET * xi1)

       ! Randomly select between the outgoing table for incoming
       ! energy E_i and E_(i+1)
       if (rang() < r) then
          l = i + 1
       else
          l = i
       end if

       loc    = 3 + 2*NR + NE + (l-1)*NET
       E_l_k  = edist % data(loc+k)
       E_l_k1 = edist % data(loc+k+1)
       xi2 = rang()
       E_out  = E_l_k + xi2*(E_l_k1 - E_l_k)

       ! TODO: Add scaled interpolation

    case (3)
       ! =============================================================
       ! INELASTIC LEVEL SCATTERING

       E_cm = edist%data(2) * (E_in - edist%data(1))
       
       E_out = E_cm

    case (4)
       ! =============================================================
       ! CONTINUOUS TABULAR DISTRIBUTION

       ! read number of interpolation regions and incoming energies 
       NR  = edist % data(1)
       NE  = edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample continuous tabular distribution."
          call error(msg)
       end if

       ! find energy bin and calculate interpolation factor -- if the
       ! energy is outside the range of the tabulated energies, choose
       ! the first or last bins
       loc = 2 + 2*NR
       if (E_in < edist % data(loc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(loc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(loc+1), NE, E_in)
          r = (E_in - edist%data(loc+i)) / & 
               & (edist%data(loc+i+1) - edist%data(loc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       xi2 = rang()
       if (r > xi2) then
          l = i + 1
       else
          l = i
       end if

       ! interpolation for energy E1 and EK
       loc   = edist%data(2 + 2*NR + NE + i)
       NP    = edist%data(loc + 2)
       E_i_1 = edist%data(loc + 2 + 1)
       E_i_K = edist%data(loc + 2 + NP)

       loc    = edist%data(2 + 2*NR + NE + i + 1)
       NP     = edist%data(loc + 2)
       E_i1_1 = edist%data(loc + 2 + 1)
       E_i1_K = edist%data(loc + 2 + NP)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       loc = edist % data(2 + 2*NR + NE + l)

       ! determine type of interpolation and number of discrete lines
       INTTp = edist % data(loc + 1)
       NP    = edist % data(loc + 2)
       if (INTTp > 10) then
          INTT = mod(INTTp,10)
          ND = (INTTp - INTT)/10
       else
          INTT = INTTp
          ND = 0
       end if

       if (ND > 0) then
          ! discrete lines present
          msg = "Discrete lines in continuous tabular distributed not &
               &yet supported"
          call error(msg)
       end if

       ! determine outgoing energy bin
       xi1 = rang()
       loc = loc + 2 ! start of EOUT
       c_k = edist % data(loc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(loc + 2*NP + k+1)
          if (xi1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = edist % data(loc+k)
       p_l_k = edist % data(loc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (xi1 - c_k)/p_l_k

       elseif (INTT == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          E_l_k1 = edist % data(loc+k+1)
          p_l_k1 = edist % data(loc+NP+k+1)

          frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
          if (frac == ZERO) then
             E_out = E_l_k + (xi1 - c_k)/p_l_k
          else
             E_out = E_l_k + (sqrt(p_l_k*p_l_k + 2*frac*(xi1 - c_k)) - & 
                  & p_l_k)/frac
          end if
       else
          msg = "Unknown interpolation type: " // trim(int_to_str(INTT))
          call error(msg)
       end if

       ! Now interpolate between incident energy bins i and i + 1
       if (l == i) then
          E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
       else
          E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
       end if

    case (5)
       ! =============================================================
       ! GENERAL EVAPORATION SPECTRUM

    case (7)
       ! =============================================================
       ! MAXWELL FISSION SPECTRUM

       ! read number of interpolation regions and incoming energies 
       NR  = edist % data(1)
       NE  = edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Maxwell fission spectrum."
          call error(msg)
       end if

       ! find incident energy bin and calculate interpolation factor
       loc = 2 + 2*NR
       if (E_in < edist % data(loc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(loc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(loc+1), NE, E_in)
          r = (E_in - edist%data(loc+i)) / & 
               & (edist%data(loc+i+1) - edist%data(loc+i))
       end if

       ! determine nuclear temperature from tabulated function
       loc = loc + NE
       T = edist%data(loc+i) + r * &
            & (edist%data(loc+i+1) - edist%data(loc+i))
       
       ! sample maxwell fission spectrum
       E_out = maxwell_spectrum(T)
       
    case (9)
       ! =============================================================
       ! EVAPORATION SPECTRUM

       ! read number of interpolation regions and incoming energies 
       NR  = edist % data(1)
       NE  = edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample evaporation spectrum."
          call error(msg)
       end if

       ! find energy bin and calculate interpolation factor -- if the
       ! energy is outside the range of the tabulated energies, choose
       ! the first or last bins
       loc = 2 + 2*NR
       if (E_in < edist % data(loc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(loc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(loc+1), NE, E_in)
          r = (E_in - edist%data(loc+i)) / & 
               & (edist%data(loc+i+1) - edist%data(loc+i))
       end if

       ! determine nuclear temperature from tabulated function
       loc = loc + NE
       T = edist%data(loc+i) + r * &
            & (edist%data(loc+i+1) - edist%data(loc+i))

       ! determine restriction energy
       loc = loc + NE
       U = edist % data(loc + 1)

       ! sample outgoing energy based on evaporation spectrum
       ! probability density function
       do
          xi1 = rang()
          xi2 = rang()
          E_out = -T * log(xi1*xi2)
          if (E_out <= E_in - U) exit
       end do
       
    case (11)
       ! =============================================================
       ! ENERGY-DEPENDENT WATT SPECTRUM

       ! read number of interpolation regions and incoming energies
       ! for parameter 'a'
       NR  = edist % data(1)
       NE  = edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Watt fission spectrum."
          call error(msg)
       end if

       ! find incident energy bin and calculate interpolation factor
       loc = 2 + 2*NR
       if (E_in < edist % data(loc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(loc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(loc+1), NE, E_in)
          r = (E_in - edist%data(loc+i)) / & 
               & (edist%data(loc+i+1) - edist%data(loc+i))
       end if

       ! determine Watt parameter 'a' from tabulated function
       loc = loc + NE
       Watt_a = edist%data(loc+i) + r * &
            & (edist%data(loc+i+1) - edist%data(loc+i))

       ! read number of interpolation regions and incoming energies
       ! for parameter 'b'
       loc = loc + NE
       NR  = edist % data(loc + 1)
       NE  = edist % data(loc + 2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Watt fission spectrum."
          call error(msg)
       end if

       ! find incident energy bin and calculate interpolation factor
       loc = loc + 2 + 2*NR
       if (E_in < edist % data(loc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(loc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(loc+1), NE, E_in)
          r = (E_in - edist%data(loc+i)) / & 
               & (edist%data(loc+i+1) - edist%data(loc+i))
       end if

       ! determine Watt parameter 'b' from tabulated function
       loc = loc + NE
       Watt_b = edist%data(loc+i) + r * &
            & (edist%data(loc+i+1) - edist%data(loc+i))

       ! Sample energy-dependent Watt fission spectrum
       E_out = watt_spectrum(Watt_a, Watt_b)

    case (44)
       ! =============================================================
       ! KALBACH-MANN CORRELATED SCATTERING

       if (.not. present(mu_out)) then
          msg = "Law 44 called without giving mu_out as argument."
          call error(msg)
       end if

       ! read number of interpolation regions and incoming energies 
       NR  = edist % data(1)
       NE  = edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Kalbach-Mann distribution."
          call error(msg)
       end if

       ! find energy bin and calculate interpolation factor -- if the
       ! energy is outside the range of the tabulated energies, choose
       ! the first or last bins
       loc = 2 + 2*NR
       if (E_in < edist % data(loc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(loc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(loc+1), NE, E_in)
          r = (E_in - edist%data(loc+i)) / & 
               & (edist%data(loc+i+1) - edist%data(loc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       xi2 = rang()
       if (r > xi2) then
          l = i + 1
       else
          l = i
       end if

       ! determine endpoints on grid i
       loc   = edist%data(2+2*NR+NE + i) ! start of LDAT for i
       NP    = edist%data(loc + 2)
       E_i_1 = edist%data(loc + 2 + 1)
       E_i_K = edist%data(loc + 2 + NP)

       ! determine endpoints on grid i+1
       loc    = edist%data(2+2*NR+NE + i+1) ! start of LDAT for i+1
       NP     = edist%data(loc + 2)
       E_i1_1 = edist%data(loc + 2 + 1)
       E_i1_K = edist%data(loc + 2 + NP)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       loc = edist % data(2 + 2*NR + NE + l)

       ! determine type of interpolation and number of discrete lines
       INTTp = edist % data(loc + 1)
       NP    = edist % data(loc + 2)
       if (INTTp > 10) then
          INTT = mod(INTTp,10)
          ND = (INTTp - INTT)/10
       else
          INTT = INTTp
          ND = 0
       end if

       if (ND > 0) then
          ! discrete lines present
          msg = "Discrete lines in continuous tabular distributed not &
               &yet supported"
          call error(msg)
       end if

       ! determine outgoing energy bin
       xi1 = rang()
       loc = loc + 2 ! start of EOUT
       c_k = edist % data(loc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(loc + 2*NP + k+1)
          if (xi1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = edist % data(loc+k)
       p_l_k = edist % data(loc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (xi1 - c_k)/p_l_k

          ! Determine Kalbach-Mann parameters
          KM_R = edist % data(loc + 3*NP + k)
          KM_A = edist % data(loc + 4*NP + k)

       elseif (INTT == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          E_l_k1 = edist % data(loc+k+1)
          p_l_k1 = edist % data(loc+NP+k+1)

          ! Find E prime
          frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
          if (frac == ZERO) then
             E_out = E_l_k + (xi1 - c_k)/p_l_k
          else
             E_out = E_l_k + (sqrt(p_l_k*p_l_k + 2*frac*(xi1 - c_k)) - & 
                  & p_l_k)/frac
          end if

          ! Determine Kalbach-Mann parameters
          R_k  = edist % data(loc + 3*NP + k)
          R_k1 = edist % data(loc + 3*NP + k+1)
          A_k  = edist % data(loc + 4*NP + k)
          A_k1 = edist % data(loc + 4*NP + k+1)
          
          KM_R = R_k + (R_k1 - R_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
          KM_A = A_k + (A_k1 - A_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
       else
          msg = "Unknown interpolation type: " // trim(int_to_str(INTT))
          call error(msg)
       end if

       ! Now interpolate between incident energy bins i and i + 1
       if (l == i) then
          E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
       else
          E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
       end if

       ! Sampled correlated angle from Kalbach-Mann parameters
       xi3 = rang()
       xi4 = rang()
       T = (TWO*xi4 - ONE) * sinh(KM_A)
       if (xi3 > KM_R) then
          mu_out = log(T + sqrt(T*T + ONE))/KM_A
       else
          mu_out = log(xi4*exp(KM_A) + (ONE - xi4)*exp(-KM_A))/KM_A
       end if

    case (61)
       ! =============================================================
       ! CORRELATED ENERGY AND ANGLE DISTRIBUTION

       if (.not. present(mu_out)) then
          msg = "Law 44 called without giving mu_out as argument."
          call error(msg)
       end if

       ! read number of interpolation regions and incoming energies 
       NR  = edist % data(1)
       NE  = edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample correlated energy-angle distribution."
          call error(msg)
       end if

       ! find energy bin and calculate interpolation factor -- if the
       ! energy is outside the range of the tabulated energies, choose
       ! the first or last bins
       loc = 2 + 2*NR
       if (E_in < edist % data(loc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(loc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(loc+1), NE, E_in)
          r = (E_in - edist%data(loc+i)) / & 
               & (edist%data(loc+i+1) - edist%data(loc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       xi2 = rang()
       if (r > xi2) then
          l = i + 1
       else
          l = i
       end if

       ! determine endpoints on grid i
       loc   = edist%data(2+2*NR+NE + i) ! start of LDAT for i
       NP    = edist%data(loc + 2)
       E_i_1 = edist%data(loc + 2 + 1)
       E_i_K = edist%data(loc + 2 + NP)

       ! determine endpoints on grid i+1
       loc    = edist%data(2+2*NR+NE + i+1) ! start of LDAT for i+1
       NP     = edist%data(loc + 2)
       E_i1_1 = edist%data(loc + 2 + 1)
       E_i1_K = edist%data(loc + 2 + NP)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       loc = edist % data(2 + 2*NR + NE + l)

       ! determine type of interpolation and number of discrete lines
       INTTp = edist % data(loc + 1)
       NP    = edist % data(loc + 2)
       if (INTTp > 10) then
          INTT = mod(INTTp,10)
          ND = (INTTp - INTT)/10
       else
          INTT = INTTp
          ND = 0
       end if

       if (ND > 0) then
          ! discrete lines present
          msg = "Discrete lines in continuous tabular distributed not &
               &yet supported"
          call error(msg)
       end if

       ! determine outgoing energy bin
       xi1 = rang()
       loc = loc + 2 ! start of EOUT
       c_k = edist % data(loc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(loc + 2*NP + k+1)
          if (xi1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = edist % data(loc+k)
       p_l_k = edist % data(loc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (xi1 - c_k)/p_l_k

       elseif (INTT == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          E_l_k1 = edist % data(loc+k+1)
          p_l_k1 = edist % data(loc+NP+k+1)

          ! Find E prime
          frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
          if (frac == ZERO) then
             E_out = E_l_k + (xi1 - c_k)/p_l_k
          else
             E_out = E_l_k + (sqrt(p_l_k*p_l_k + 2*frac*(xi1 - c_k)) - & 
                  & p_l_k)/frac
          end if
       else
          msg = "Unknown interpolation type: " // trim(int_to_str(INTT))
          call error(msg)
       end if

       ! Now interpolate between incident energy bins i and i + 1
       if (l == i) then
          E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
       else
          E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
       end if

       ! Find location of correlated angular distribution
       loc = edist % data(loc+3*NP+k)

       ! Check if angular distribution is isotropic
       if (loc == 0) then
          mu_out = TWO * rang() - ONE
          return
       end if

       ! interpolation type and number of points in angular distribution
       JJ = edist % data(loc + 1)
       NP = edist % data(loc + 2)

       ! determine outgoing cosine bin
       xi3 = rang()
       loc = loc + 2
       c_k = edist % data(loc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(loc + 2*NP + k+1)
          if (xi3 < c_k1) exit
          c_k = c_k1
       end do

       p_k  = edist % data(loc + NP + k)
       mu_k = edist % data(loc + k)
       if (JJ == HISTOGRAM) then
          ! Histogram interpolation
          mu_out = mu_k + (xi3 - c_k)/p_k

       elseif (JJ == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          p_k1  = edist % data(loc + NP + k+1)
          mu_k1 = edist % data(loc + k+1)

          frac = (p_k1 - p_k)/(mu_k1 - mu_k)
          if (frac == ZERO) then
             mu_out = mu_k + (xi3 - c_k)/p_k
          else
             mu_out = mu_k + (sqrt(p_k*p_k + 2*frac*(xi3 - c_k))-p_k)/frac
          end if
       else
          msg = "Unknown interpolation type: " // trim(int_to_str(JJ))
          call error(msg)
       end if

    case (66)
       ! =============================================================
       ! N-BODY PHASE SPACE DISTRIBUTION

    case (67)
       ! =============================================================
       ! LABORATORY ENERGY-ANGLE LAW

    end select
    
  end subroutine sample_energy

!=====================================================================
! MAXWELL_SPECTRUM samples an energy from the Maxwell fission
! distribution based on a rejection sampling scheme. This is described
! in the MCNP manual volume I -- need to verify formula
!=====================================================================

  function maxwell_spectrum(T) result(E_out)

    real(8), intent(in)  :: T     ! tabulated function of incoming E
    real(8)              :: E_out ! sampled energy

    real(8) :: r1, r2, r3, r4  ! random numbers
    real(8) :: d               ! r1^2 + r2^2

    ! TODO: Change this to direct sampling

    r1 = rang()
    do
       r2 = rang()
       d = r1*r1 + r2*r2
       if (d < ONE) exit
       r1 = r2
    end do

    r3 = rang()
    r4 = rang()
    E_out = -T*(r1*r1 * log(r3) / d + log(r4))

  end function maxwell_spectrum

!=====================================================================
! WATT_SPECTRUM samples the outgoing energy from a Watt
! energy-dependent fission spectrum. Although fitted parameters exist
! for many nuclides, generally the continuous tabular distributions
! (LAW 4) should be used in lieu of the Watt spectrum
!=====================================================================

  function watt_spectrum(a, b) result(E_out)

    real(8), intent(in) :: a     ! Watt parameter a
    real(8), intent(in) :: b     ! Watt parameter b
    real(8)             :: E_out ! energy of emitted neutron

    real(8) :: g      ! from sampling scheme R28
    real(8) :: r1, r2 ! random numbers

    g = sqrt((1 + a*b/8)**2 - 1) + (1 + a*b/8)
    do
       r1 = log(rang())
       r2 = log(rang())
       E_out = -a*g*r1
       if (((1 - g)*(1 - r1) - r2)**2 < b*E_out) exit
    end do

  end function watt_spectrum

!=====================================================================
! WIGNER samples a Wigner distribution of energy level spacings. Note
! that this scheme is C50 in the Monte Carlo Sampler from Los Alamos
! (LA-9721-MS).
!=====================================================================

  function wigner(D_avg) result (D)

    real(8), intent(in) :: D_avg ! average level spacing
    real(8)             :: D     ! sampled level spacing

    real(8) :: c

    c = -4.0*D_avg*D_avg/PI * log(rang())
    D = sqrt(c)

  end function wigner

!=====================================================================
! CHI_SQUARED samples a chi-squared distribution with n degrees of
! freedom. The distribution of resonance widths in the unresolved
! region is given by a chi-squared distribution. For the special case
! of n=1, this is a Porter-Thomas distribution. This currently
! hardwires cases for n = 1 to 4, but could be made more general by
! treating all n odd with rule C64 and all n even with rule C45.
!=====================================================================

  function chi_squared(n, G_avg) result(G)

    real(8), intent(in)           :: n     ! number of degrees of freedom
    real(8), intent(in), optional :: G_avg ! average resonance width

    real(8) :: G          ! sampled random variable (or resonance width)
    real(8) :: x, y       ! dummy variables
    real(8) :: r1, r2, r3 ! psuedorandom numbers

    select case (n)
    case (1)
       ! For this case, we get p(x) = c*x^(-1/2)*exp(-x). This can be
       ! sampled with rule C64
       r1 = rang()
       r2 = rang()
       y = cos(PI/2*r2)
       x = -2*log(r1)*y*y
       G = x * G_avg
    case (2)
       ! This case corresponds to a simple exponential distribution
       ! p(x) = exp(-x) whose CDF can be directly inverted.
       r1 = rang()
       x = -log(r1)
    case (3)
       ! For this case, we get p(x) = c*x^(1/2)*exp(-x). This can be
       ! sampled with rule C64
       r1 = rang()
       r2 = rang()
       r3 = rang()
       y = cos(PI/2*r2)
       x = -2.0/3.0*(log(r1)*y*y + log(r3))
    case (4)
       ! For this case, we get p(x) = c*x*exp(-x). This can be sampled
       ! with rule C45.
       r1 = rang()
       r2 = rang()
       x = -0.5*log(r1*r2)
    end select

    ! If sampling a chi-squared distribution for a resonance width and
    ! the average resonance width has been given, return the sampled
    ! resonance width.
    if (present(G_avg)) then
       G = x * G_avg
    else
       G = x
    end if

  end function chi_squared

end module physics
