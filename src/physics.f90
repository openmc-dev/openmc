module physics

  use global
  use geometry,    only: find_cell, dist_to_boundary, cross_boundary
  use types,       only: Particle
  use mcnp_random, only: rang
  use output,      only: error, message
  use search,      only: binary_search
  use endf,        only: reaction_name

  implicit none

contains

!=====================================================================
! TRANSPORT encompasses the main logic for moving a particle through
! geometry.
!=====================================================================

  subroutine transport(p)

    type(Particle), pointer :: p

    character(250) :: msg
    integer :: i, surf
    logical :: alive
    logical :: found_cell
 
    real(8) :: d_to_boundary
    real(8) :: d_to_collision
    real(8) :: distance
    real(8) :: Sigma          ! total cross-section
    real(8) :: f              ! interpolation factor
    integer :: IE             ! index on energy grid
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

    if (verbosity >= 10) then
       msg = "=== Particle " // trim(int_to_str(p % uid)) // " ==="
       call message(msg, 10)

       i = cells(p % cell)%uid
       msg = "    Born in cell " // trim(int_to_str(i))
       call message(msg, 10)
    end if

    ! sample energy from Watt fission energy spectrum for U-235
    p % E = watt_spectrum(0.988_8, 2.249_8)

    ! find energy index, interpolation factor
    call find_energy_index(p)
    IE = p % IE
    f  = p % interp

    ! Determine material total cross-section
    Sigma = f*cMaterial%total_xs(IE) + (1-f)*cMaterial%total_xs(IE+1)

    do while (p % alive)

       ! Find the distance to the nearest boundary
       call dist_to_boundary(p, d_to_boundary, surf)

       ! Sample a distance to collision
       d_to_collision = -log(rang()) / 1.0 ! Sigma
       
       ! Select smaller of the two distances
       distance = min(d_to_boundary, d_to_collision)

       ! Advance particle
       p%xyz = p%xyz + distance * p%uvw

       ! Add pathlength tallies

       if (d_to_collision > d_to_boundary) then
          p % surface = surf
          p % cell = 0
          call cross_boundary(p)
       else
          ! collision
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

    integer :: IE
    real(8) :: E
    real(8) :: interp

    ! copy particle's energy
    E = p % E

    ! if particle's energy is outside of energy grid range, set to
    ! first or last index. Otherwise, do a binary search through the
    ! union energy grid.
    if (E < e_grid(1)) then
       IE = 1
    elseif (E > e_grid(n_grid)) then
       IE = n_grid
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
! COLLISION
!=====================================================================

  subroutine collision(p)

    type(Particle), pointer :: p

    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn
    real(8) :: r1
    character(250) :: msg
    integer :: i,j
    integer :: n_isotopes
    integer :: IE
    real(8) :: f, Sigma, total
    real(8) :: density, density_i
    real(8) :: prob
    real(8), allocatable :: Sigma_t(:)

    ! tallies

    density = cMaterial%atom_density

    ! calculate total cross-section for each nuclide at current energy
    ! in order to create discrete pdf for sampling nuclide
    n_isotopes = cMaterial%n_isotopes
    allocate(Sigma_t(n_isotopes))
    do i = 1, n_isotopes
       table => xs_continuous(cMaterial%table(i))
       density_i = cMaterial%atom_percent(i)*density

       ! search nuclide energy grid
       IE = table%grid_index(p % IE)
       f = (p%E - table%energy(IE))/(table%energy(IE+1) - table%energy(IE))

       Sigma = density_i*(f*table%sigma_t(IE) + (1-f)*(table%sigma_t(IE+1)))
       Sigma_t(i) = Sigma
    end do

    ! sample nuclide
    r1 = rang()
    prob = ZERO
    total = sum(Sigma_t)
    do i = 1, n_isotopes
       prob = prob + Sigma_t(i) / total
       if (r1 < prob) exit
    end do

    ! Get table, total xs, interpolation factor
    table => xs_continuous(cMaterial%table(i))
    Sigma = Sigma_t(i)
    IE = table%grid_index(p % IE)
    f = (p%E - table%energy(IE))/(table%energy(IE+1) - table%energy(IE))
    density = cMaterial%atom_percent(i)*density

    ! free memory
    deallocate(Sigma_t)

    ! sample reaction channel
    r1 = rang()*Sigma
    prob = ZERO
    do i = 1, table%n_reaction
       rxn => table%reactions(i)
       if (rxn%MT >= 200) cycle
       if (IE < rxn%IE) cycle
       prob = prob + density * (f*rxn%sigma(IE-rxn%IE+1) & 
            & + (1-f)*(rxn%sigma(IE-rxn%IE+2)))
       if (r1 < prob) exit
    end do
    if (verbosity >= 10) then
       msg = "    " // trim(reaction_name(rxn%MT)) // " with nuclide " // &
            & trim(table%name)
       call message(msg, 10)
    end if

    ! call appropriate subroutine
    select case (rxn%MT)
    case (2)
       call elastic_scatter(p, table%awr)
    case (18:21, 38)
       call n_fission(p, table, rxn)
    case (51:90)
       call level_scatter(p, table, rxn)
    case (102)
       call n_gamma(p)
    case default
       call elastic_scatter(p, table%awr)
    end select
    
  end subroutine collision

!=====================================================================
! ELASTIC_SCATTER
!=====================================================================

  subroutine elastic_scatter(p, awr)

    type(Particle), pointer :: p
    real(8), intent(in) :: awr

    real(8) :: phi ! azimuthal angle
    real(8) :: mu  ! cosine of polar angle
    real(8) :: vx, vy, vz
    real(8) :: vcx, vcy ,vcz
    real(8) :: vel
    real(8) :: u, v, w
    real(8) :: E
    integer :: IE

    vel = sqrt(p % E)

    vx = vel*p % uvw(1)
    vy = vel*p % uvw(2)
    vz = vel*p % uvw(3)


    vcx = vx/(awr + ONE)
    vcy = vy/(awr + ONE)
    vcz = vz/(awr + ONE)

    ! Transform to CM frame
    vx = vx - vcx
    vy = vy - vcy
    vz = vz - vcz

    vel = sqrt(vx*vx + vy*vy + vz*vz)

    ! Select isotropic direcion -- this is only valid for s-wave
    ! scattering
    phi = TWO*PI*rang()
    mu = TWO*rang() - ONE
    u = mu
    v = sqrt(ONE - mu**2) * cos(phi)
    w = sqrt(ONE - mu**2) * sin(phi)

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

    ! find energy index, interpolation factor
    call find_energy_index(p)

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
    integer :: j           ! index on nu energy grid
    integer :: loc         ! index before start of energies/nu values
    integer :: NC          ! number of polynomial coefficients
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies tabulated
    real(8) :: E           ! incoming energy of neutron
    real(8) :: E_out       ! outgoing energy of fission neutron
    real(8) :: mu_out      ! outgoing angle of fission neutron (if needed)
    real(8) :: c           ! polynomial coefficient
    real(8) :: f           ! interpolation factor
    real(8) :: nu_total    ! total nu
    real(8) :: nu_prompt   ! prompt nu
    real(8) :: nu_delay    ! delayed nu
    integer :: nu          ! actual number of neutrons produced
    real(8) :: mu          ! fission neutron angular cosine
    real(8) :: phi         ! fission neutron azimuthal angle
    real(8) :: beta        ! delayed neutron fraction
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
       j = binary_search(table % nu_t_data(loc+1), NE, E)
       f = (E - table % nu_t_data(loc+j)) / &
            & (table % nu_t_data(loc+j+1) - table % nu_t_data(loc+j))

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
       j = binary_search(table % nu_p_data(loc+1), NE, E)
       f = (E - table % nu_p_data(loc+j)) / &
            & (table % nu_p_data(loc+j+1) - table % nu_p_data(loc+j))

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
       j = binary_search(table % nu_d_data(loc+1), NE, E)
       f = (E - table % nu_d_data(loc+j)) / &
            & (table % nu_d_data(loc+j+1) - table % nu_d_data(loc+j))

       ! determine nu total
       loc = loc + NE
       nu_delay = table % nu_d_data(loc+j) + f * & 
            & (table % nu_d_data(loc+j+1) - table % nu_d_data(loc+j))
    end if

    beta = nu_delay / nu_total

    ! TODO: Heat generation from fission

    ! Sample number of neutrons produced
    nu_total = p % wgt * nu_total
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
       if (rxn % has_angle_dist) then
          mu = sample_angle(rxn, E)
       else
          mu = TWO * rang() - ONE
       end if

       ! Sample azimuthal angle uniformly in [0,2*pi)
       phi = TWO*PI*rang()
       fission_bank(i) % uvw(1) = mu
       fission_bank(i) % uvw(2) = sqrt(1. - mu**2) * cos(phi)
       fission_bank(i) % uvw(3) = sqrt(1. - mu**2) * sin(phi)

       ! determine energy of fission neutron
       call sample_energy(rxn, E, E_out, mu_out)
       fission_bank(i) % E = E_out
    end do

    ! increment number of bank sites
    n_bank = min(n_bank + nu, 3*n_particles)

    ! kill original neutron
    p % alive = .false.

  end subroutine n_fission

!=====================================================================
! LEVEL_SCATTER
!=====================================================================

  subroutine level_scatter(p, table, rxn)

    type(Particle),      pointer :: p
    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn

    real(8) :: A          ! atomic weight ratio of nuclide
    real(8) :: E_in       ! incoming energy
    real(8) :: mu_cm      ! cosine of scattering angle in center-of-mass
    real(8) :: mu_lab     ! cosine of scattering angle in laboratory
    real(8) :: E_cm       ! outgoing energy in center-of-mass
    real(8) :: E_lab      ! outgoing energy in laboratory
    character(250) :: msg ! error message
    

    ! copy energy of neutron
    E_in = p % E

    ! determine A
    A = table % awr

    ! determine if scattering is in CM (it should be!)
    if (rxn % TY < 0) then
       ! scattering angle in center-of-mass
       mu_cm = sample_angle(rxn, E_in)
    else
       msg = "Level inelastic scattering should not sample angle &
            &in laboratory system!"
       call error(msg)
    end if

    ! sample outgoing energy in center-of-mass
    call sample_energy(rxn, E_in, E_cm)

    ! determine outgoing energy in lab
    E_lab = E_cm + (E_in + TWO * mu_cm * (A+ONE) * sqrt(E_in * E_cm)) & 
         & / ((A+ONE)*(A+ONE))

    ! determine outgoing angle in lab
    mu_lab = mu_cm * sqrt(E_cm/E_lab) + ONE/(A+ONE) * sqrt(E_in/E_lab)

    ! change direction of particle
    call rotate_angle(p, mu_lab)

    ! change energy of particle
    p % E = E_lab

    ! find energy index, interpolation factor
    call find_energy_index(p)

  end subroutine level_scatter

!=====================================================================
! N_GAMMA
!=====================================================================

  subroutine n_gamma(p)

    type(Particle), pointer :: p

    integer :: cell_num
    character(250) :: msg

    p % alive = .false.
    if (verbosity >= 10) then
       cell_num = cells(p % cell)%uid
       msg = "    Absorbed in cell " // trim(int_to_str(cell_num))
       call message(msg, 10)
    end if

  end subroutine n_gamma

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

    ! determine number of incoming energies
    n = rxn % adist % n_energy

    ! find energy bin and calculate interpolation factor -- if the
    ! energy is outside the range of the tabulated energies, choose
    ! the first or last bins
    if (E < rxn % adist % energy(1)) then
       i = 1
       r = 0.0
    elseif (E > rxn % adist % energy(n)) then
       i = n - 1
       r = 1.0
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

       elseif (interp == LINEARLINEAR) then
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

  subroutine rotate_angle(p, mu)

    type(Particle), pointer :: p
    real(8), intent(in)     :: mu ! cosine of angle in lab

    real(8) :: phi, sinphi, cosphi
    real(8) :: a,b
    real(8) :: u0, v0, w0

    ! Copy original directional cosines
    u0 = p % uvw(1)
    v0 = p % uvw(2)
    w0 = p % uvw(3)

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
       p % uvw(1) = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
       p % uvw(2) = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
       p % uvw(3) = mu*w0 - a*b*cosphi
    else
       b = sqrt(ONE - v0*v0)
       p % uvw(1) = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
       p % uvw(2) = mu*v0 - a*b*cosphi
       p % uvw(3) = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end subroutine rotate_angle
    
!=====================================================================
! SAMPLE_ENERGY
!=====================================================================

  subroutine sample_energy(rxn, E_in, E_out, mu_out)

    type(AceReaction), pointer       :: rxn
    real(8), intent(in)              :: E_in
    real(8), intent(out)             :: E_out
    real(8), intent(inout), optional :: mu_out

    integer :: i           ! index on incoming energy grid
    integer :: k           ! sampled index on outgoing grid
    integer :: l           ! sampled index on incoming grid
    integer :: loc         ! dummy index
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies
    integer :: NET         ! number of outgoing energies
    integer :: INTTp       ! combination of INTT and ND
    integer :: INTT        ! 1 = histogram, 2 = linear-linear
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
    if (rxn % edist % n_interp > 0) then
       msg = "Multiple interpolation regions not supported while &
            &attempting to sampling secondary energy distribution."
       call error(msg)
    end if
       
    ! Determine which secondary energy distribution law to use
    select case (rxn % edist % law)
    case (1)
       ! =============================================================
       ! TABULAR EQUIPROBABLE ENERGY BINS

       ! read number of interpolation regions, incoming energies, and
       ! outgoing energies
       NR  = rxn % edist % data(1)
       NE  = rxn % edist % data(2 + 2*NR)
       NET = rxn % edist % data(3 + 2*NR + NE)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample equiprobable energy bins."
          call error(msg)
       end if

       ! determine index on incoming energy grid and interpolation
       ! factor
       loc = 2 + 2*NR
       i = binary_search(rxn % edist % data(loc+1), NE, E_in)
       r = (E_in - rxn%edist%data(loc+i)) / &
            & (rxn%edist%data(loc+i+1) - rxn%edist%data(loc+i))

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
       E_l_k  = rxn % edist % data(loc+k)
       E_l_k1 = rxn % edist % data(loc+k+1)
       xi2 = rang()
       E_out  = E_l_k + xi2*(E_l_k1 - E_l_k)

       ! TODO: Add scaled interpolation

    case (3)
       ! =============================================================
       ! INELASTIC LEVEL SCATTERING

       E_cm = rxn%edist%data(2) * (E_in - rxn%edist%data(1))
       
       E_out = E_cm

    case (4)
       ! =============================================================
       ! CONTINUOUS TABULAR DISTRIBUTION

       ! read number of interpolation regions and incoming energies 
       NR  = rxn % edist % data(1)
       NE  = rxn % edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample continuous tabular distribution."
          call error(msg)
       end if

       ! find energy bin and calculate interpolation factor -- if the
       ! energy is outside the range of the tabulated energies, choose
       ! the first or last bins
       loc = 2 + 2*NR
       if (E_in < rxn % edist % data(loc+1)) then
          i = 1
          r = 0.0
       elseif (E_in > rxn % edist % energy(loc+NE)) then
          i = NE - 1
          r = 1.0
       else
          i = binary_search(rxn % edist % energy(loc+1), NE, E_in)
          r = (E_in - rxn%edist%energy(loc+i)) / & 
               & (rxn%edist%energy(loc+i+1) - rxn%edist%energy(loc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       xi2 = rang()
       if (r > xi2) then
          l = i + 1
       else
          l = i
       end if

       ! interpolation on grid i for energy E1
       loc = rxn % edist % data(2 + 2*NR + NE + i) + 2 ! start of EOUT(i)
       E_i_1 = rxn%edist%data(loc + 1)
       E_i_K = rxn%edist%data(loc + NE)

       loc = rxn % edist % data(2 + 2*NR + NE + i + 1) + 2 ! start of EOUT(i+1)
       E_i1_1 = rxn%edist%data(loc + 1)
       E_i1_K = rxn%edist%data(loc + NE)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       loc = rxn % edist % data(2 + 2*NR + NE + l)

       ! determine type of interpolation and number of discrete lines
       INTTp = rxn % edist % data(loc + 1)
       NP    = rxn % edist % data(loc + 2)
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
       c_k = rxn % edist % data(loc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = rxn % edist % data(loc + 2*NP + k+1)
          if (xi1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = rxn % edist % data(loc+k)
       p_l_k = rxn % edist % data(loc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (xi1 - c_k)/p_l_k

       elseif (INTT == LINEARLINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          E_l_k1 = rxn % edist % data(loc+k+1)
          p_l_k1 = rxn % edist % data(loc+NP+k+1)

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

       ! Now interpolate between incident enregy bins i and i + 1
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
       NR  = rxn % edist % data(1)
       NE  = rxn % edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Maxwell fission spectrum."
          call error(msg)
       end if

       ! find incident energy bin and calculate interpolation factor
       loc = 2 + 2*NR
       if (E_in < rxn % edist % data(loc+1)) then
          i = 1
          r = 0.0
       elseif (E_in > rxn % edist % energy(loc+NE)) then
          i = NE - 1
          r = 1.0
       else
          i = binary_search(rxn % edist % energy(loc+1), NE, E_in)
          r = (E_in - rxn%edist%energy(loc+i)) / & 
               & (rxn%edist%energy(loc+i+1) - rxn%edist%energy(loc+i))
       end if

       ! determine nuclear temperature from tabulated function
       loc = loc + NE
       T = rxn%edist%data(loc+i) + r * &
            & (rxn%edist%data(loc+i+1) - rxn%edist%data(loc+i))
       
       ! sample maxwell fission spectrum
       E_out = maxwell_spectrum(T)
       
    case (9)
       ! =============================================================
       ! EVAPORATION SPECTRUM

       ! read number of interpolation regions and incoming energies 
       NR  = rxn % edist % data(1)
       NE  = rxn % edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample evaporation spectrum."
          call error(msg)
       end if

       ! find energy bin and calculate interpolation factor -- if the
       ! energy is outside the range of the tabulated energies, choose
       ! the first or last bins
       loc = 2 + 2*NR
       if (E_in < rxn % edist % data(loc+1)) then
          i = 1
          r = 0.0
       elseif (E_in > rxn % edist % energy(loc+NE)) then
          i = NE - 1
          r = 1.0
       else
          i = binary_search(rxn % edist % energy(loc+1), NE, E_in)
          r = (E_in - rxn%edist%energy(loc+i)) / & 
               & (rxn%edist%energy(loc+i+1) - rxn%edist%energy(loc+i))
       end if

       ! determine nuclear temperature from tabulated function
       loc = loc + NE
       T = rxn%edist%data(loc+i) + r * &
            & (rxn%edist%data(loc+i+1) - rxn%edist%data(loc+i))

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
       NR  = rxn % edist % data(1)
       NE  = rxn % edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Watt fission spectrum."
          call error(msg)
       end if

       ! find incident energy bin and calculate interpolation factor
       loc = 2 + 2*NR
       if (E_in < rxn % edist % data(loc+1)) then
          i = 1
          r = 0.0
       elseif (E_in > rxn % edist % energy(loc+NE)) then
          i = NE - 1
          r = 1.0
       else
          i = binary_search(rxn % edist % energy(loc+1), NE, E_in)
          r = (E_in - rxn%edist%energy(loc+i)) / & 
               & (rxn%edist%energy(loc+i+1) - rxn%edist%energy(loc+i))
       end if

       ! determine Watt parameter 'a' from tabulated function
       loc = loc + NE
       Watt_a = rxn%edist%data(loc+i) + r * &
            & (rxn%edist%data(loc+i+1) - rxn%edist%data(loc+i))

       ! read number of interpolation regions and incoming energies
       ! for parameter 'b'
       loc = loc + NE
       NR  = rxn % edist % data(loc + 1)
       NE  = rxn % edist % data(loc + 2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Watt fission spectrum."
          call error(msg)
       end if

       ! find incident energy bin and calculate interpolation factor
       loc = loc + 2 + 2*NR
       if (E_in < rxn % edist % data(loc+1)) then
          i = 1
          r = 0.0
       elseif (E_in > rxn % edist % energy(loc+NE)) then
          i = NE - 1
          r = 1.0
       else
          i = binary_search(rxn % edist % energy(loc+1), NE, E_in)
          r = (E_in - rxn%edist%energy(loc+i)) / & 
               & (rxn%edist%energy(loc+i+1) - rxn%edist%energy(loc+i))
       end if

       ! determine Watt parameter 'b' from tabulated function
       loc = loc + NE
       Watt_b = rxn%edist%data(loc+i) + r * &
            & (rxn%edist%data(loc+i+1) - rxn%edist%data(loc+i))

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
       NR  = rxn % edist % data(1)
       NE  = rxn % edist % data(2 + 2*NR)
       if (NR > 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to sample Kalbach-Mann distribution."
          call error(msg)
       end if

       ! find energy bin and calculate interpolation factor -- if the
       ! energy is outside the range of the tabulated energies, choose
       ! the first or last bins
       loc = 2 + 2*NR
       if (E_in < rxn % edist % data(loc+1)) then
          i = 1
          r = 0.0
       elseif (E_in > rxn % edist % energy(loc+NE)) then
          i = NE - 1
          r = 1.0
       else
          i = binary_search(rxn % edist % energy(loc+1), NE, E_in)
          r = (E_in - rxn%edist%energy(loc+i)) / & 
               & (rxn%edist%energy(loc+i+1) - rxn%edist%energy(loc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       xi2 = rang()
       if (r > xi2) then
          l = i + 1
       else
          l = i
       end if

       ! interpolation on grid i for energy E1
       loc = rxn % edist % data(2 + 2*NR + NE + i) + 2 ! start of EOUT(i)
       E_i_1 = rxn%edist%data(loc + 1)
       E_i_K = rxn%edist%data(loc + NE)

       loc = rxn % edist % data(2 + 2*NR + NE + i + 1) + 2 ! start of EOUT(i+1)
       E_i1_1 = rxn%edist%data(loc + 1)
       E_i1_K = rxn%edist%data(loc + NE)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       loc = rxn % edist % data(2 + 2*NR + NE + l)

       ! determine type of interpolation and number of discrete lines
       INTTp = rxn % edist % data(loc + 1)
       NP    = rxn % edist % data(loc + 2)
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
       c_k = rxn % edist % data(loc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = rxn % edist % data(loc + 2*NP + k+1)
          if (xi1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = rxn % edist % data(loc+k)
       p_l_k = rxn % edist % data(loc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (xi1 - c_k)/p_l_k

          ! Determine Kalbach-Mann parameters
          KM_R = rxn % edist % data(loc + 3*NP + k)
          KM_A = rxn % edist % data(loc + 4*NP + k)

       elseif (INTT == LINEARLINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          E_l_k1 = rxn % edist % data(loc+k+1)
          p_l_k1 = rxn % edist % data(loc+NP+k+1)

          ! Find E prime
          frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
          if (frac == ZERO) then
             E_out = E_l_k + (xi1 - c_k)/p_l_k
          else
             E_out = E_l_k + (sqrt(p_l_k*p_l_k + 2*frac*(xi1 - c_k)) - & 
                  & p_l_k)/frac
          end if

          ! Determine Kalbach-Mann parameters
          R_k  = rxn % edist % data(loc + 3*NP + k)
          R_k1 = rxn % edist % data(loc + 3*NP + k+1)
          A_k  = rxn % edist % data(loc + 4*NP + k)
          A_k1 = rxn % edist % data(loc + 4*NP + k+1)
          
          KM_R = R_k + (R_k1 - R_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
          KM_A = A_k + (A_k1 - A_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
       else
          msg = "Unknown interpolation type: " // trim(int_to_str(INTT))
          call error(msg)
       end if

       ! Now interpolate between incident enregy bins i and i + 1
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

    r1 = rang()
    do
       r2 = rang()
       d = r1*r1 + r2*r2
       if (d < ONE) exit
       r1 = r2
    end do

    r3 = rang()
    r4 = rang()
    E_out = -T*(r1**2 * log(r3) / d + log(r4))

  end function maxwell_spectrum

!=====================================================================
! WATT_SPECTRUM samples the outgoing energy from a Watt
! energy-dependent fission spectrum. Although fitted parameters exist
! for many nuclides, generally the continuous tabular distributions
! (LAW 4) should be used in lieu of the Watt spectrum
!=====================================================================

  function watt_spectrum(a, b) result(E_out)

    real(8), intent(in) :: a
    real(8), intent(in) :: b
    real(8)             :: E_out

    real(8) :: g
    real(8) :: r1, r2

    g = sqrt((1 + a*b/8)**2 - 1) + (1 + a*b/8)
    do
       r1 = log(rang())
       r2 = log(rang())
       E_out = -a*g*r1
       if (((1 - g)*(1 - r1) - r2)**2 < b*E_out) exit
    end do

  end function watt_spectrum

end module physics
