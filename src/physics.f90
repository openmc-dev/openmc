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
    prob = 0.0_8
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
    prob = 0.0_8
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


    vcx = vx/(awr + 1.0_8)
    vcy = vy/(awr + 1.0_8)
    vcz = vz/(awr + 1.0_8)

    ! Transform to CM frame
    vx = vx - vcx
    vy = vy - vcy
    vz = vz - vcz

    vel = sqrt(vx*vx + vy*vy + vz*vz)

    ! Select isotropic direcion -- this is only valid for s-wave
    ! scattering
    phi = 2.0_8*PI*rang()
    mu = 2.0_8*rang() - 1.0_8
    u = mu
    v = sqrt(1.0_8 - mu**2) * cos(phi)
    w = sqrt(1.0_8 - mu**2) * sin(phi)

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
! LEVEL_INELASTIC
!=====================================================================

  subroutine level_inelastic

  end subroutine level_inelastic

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

  subroutine sample_angle(rxn, E)

    type(AceReaction), pointer    :: rxn ! reaction
    real(8),           intent(in) :: E   ! incoming energy

    real(8)        :: xi      ! random number on [0,1)
    integer        :: interp  ! type of interpolation
    integer        :: type    ! angular distribution type
    integer        :: i       ! incoming energy bin
    integer        :: n       ! number of incoming energy bins
    integer        :: loc     ! location in data array
    integer        :: np      ! number of points in cos distribution
    integer        :: bin     ! cosine bin
    real(8)        :: f       ! interpolation factor
    real(8)        :: mu0     ! cosine in bin b
    real(8)        :: mu1     ! cosine in bin b+1
    real(8)        :: mu      ! final cosine sampled
    real(8)        :: c       ! cumulative distribution frequency
    real(8)        :: p0,p1   ! probability distribution
    character(250) :: msg     ! error message

    ! determine number of incoming energies
    n = rxn % adist_n_energy

    ! find energy bin and calculate interpolation factor -- if the
    ! energy is outside the range of the tabulated energies, choose
    ! the first or last bins
    if (E < rxn%adist_energy(1)) then
       i = 1
       f = 0.0
    elseif (E > rxn%adist_energy(n)) then
       i = n - 1
       f = 1.0
    else
       i = binary_search(rxn%adist_energy, n, E)
       f = (E - rxn % adist_energy(i)) / & 
            & (rxn % adist_energy(i+1) - rxn % adist_energy(i))
    end if

    ! Sample between the ith and (i+1)th bin
    if (f > rang()) i = i + 1

    ! check whether this is a 32-equiprobable bin or a tabular
    ! distribution
    loc  = rxn % adist_location(i)
    type = rxn % adist_type(i)
    if (type == ANGLE_ISOTROPIC) then
       mu = 2.0_8 * rang() - 1
    elseif (type == ANGLE_32_EQUI) then
       ! sample cosine bin
       xi = rang()
       bin = 1 + int(32.0_8*xi)

       ! calculate cosine
       mu0 = rxn % adist_data(loc + bin - 1)
       mu1 = rxn % adist_data(loc + bin)
       mu = mu0 + (32.0_8 * xi - bin) * (mu1 - mu0)

    elseif (type == ANGLE_TABULAR) then
       interp = rxn % adist_data(loc)
       np     = rxn % adist_data(loc+1)

       ! determine outgoing cosine bin
       xi = rang()
       do bin = loc+2, loc+1+NP
          c = rxn % adist_data(bin+2*np)
          if (xi > c) exit
       end do

       p0  = rxn % adist_data(bin+np)
       mu0 = rxn % adist_data(bin)
       if (interp == HISTOGRAM) then
          ! Histogram interpolation
          mu = mu0 + (xi - c)/p0

       elseif (interp == LINEARLINEAR) then
          ! Linear-linear interpolation -- not sure how you come about
          ! the formula given in the MCNP manual
          p1  = rxn % adist_data(bin+np+1)
          mu1 = rxn % adist_data(bin+1)
          
          f = (p1 - p0)/(mu1 - mu0)
          if (f == ZERO) then
             mu = mu0 + (xi - c)/p0
          else
             mu = mu0 + (sqrt(p0*p0 + 2*f*(xi - c))-p0)/f
          end if
       else
          msg = "Unknown interpolation type: " // trim(int_to_str(interp))
          call error(msg)
       end if
          
    else
       msg = "Unknown angular distribution type: " // trim(int_to_str(type))
       call error(msg)
    end if
    
  end subroutine sample_angle

!=====================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose
! cosine is mu and through an azimuthal angle sampled uniformly. Note
! that this is done with direct sampling rather than rejection as is
! done in MCNP and SERPENT.
!=====================================================================

  subroutine rotate_angle(mu, u, v, w)

    real(8), intent(in)    :: mu ! cosine of angle
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v
    real(8), intent(inout) :: w

    real(8) :: phi, sinphi, cosphi
    real(8) :: a,b
    real(8) :: u0, v0, w0

    ! Copy original directional cosines
    u0 = u
    v0 = v
    w0 = w

    ! Sample azimuthal angle in [0,2pi)
    phi = 2.0_8 * PI * rang()

    ! Precompute factors to save flops
    sinphi = sin(phi)
    cosphi = cos(phi)
    a = sqrt(1.0_8 - mu*mu)
    b = sqrt(1.0_8 - w*w)

    ! Need to treat special case where sqrt(1 - w**2) is close to zero
    ! by expanding about the v component rather than the w component
    if (b > 1e-10) then
       u = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
       v = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
       w = mu*w0 - a*b*cosphi
    else
       b = sqrt(1.0_8 - v*v)
       u = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
       v = mu*v0 - a*b*cosphi
       w = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end subroutine rotate_angle
    
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
       if (d < 1) exit
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
