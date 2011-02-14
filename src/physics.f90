module physics

  use global
  use geometry,    only: find_cell, dist_to_boundary, cross_boundary
  use types,       only: Neutron
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

  subroutine transport(neut)

    type(Neutron), pointer, intent(inout) :: neut

    character(250) :: msg
    integer :: i, surf
    logical :: alive
 
    real(8) :: d_to_boundary
    real(8) :: d_to_collision
    real(8) :: distance
    real(8) :: Sigma          ! total cross-section
    real(8) :: f              ! interpolation factor
    real(8) :: tmp(3)
    integer :: IE             ! index on energy grid

    ! determine what cell the particle is in
    if (neut%cell == 0) then
       call find_cell(neut)
    end if
    if (verbosity >= 10) then
       msg = "=== Particle " // trim(int_to_str(neut%uid)) // " ==="
       call message(msg, 10)

       i = cells(neut%cell)%uid
       msg = "    Born in cell " // trim(int_to_str(i))
       call message(msg, 10)
    end if

    ! sample energy from Watt fission energy spectrum for U-235
    neut%E = watt_spectrum(0.988_8, 2.249_8)

    ! find energy index, interpolation factor
    IE = binary_search(e_grid, n_grid, neut%E)
    f = (neut%E - e_grid(IE))/(e_grid(IE+1) - e_grid(IE))

    ! Determine material total cross-section
    Sigma = f*cMaterial%total_xs(IE) + (1-f)*cMaterial%total_xs(IE+1)
    neut%IE = IE
    neut%interp = f

    do while (neut%alive)

       ! Find the distance to the nearest boundary
       call dist_to_boundary(neut, d_to_boundary, surf)

       ! Sample a distance to collision
       d_to_collision = -log(rang()) / 1.0 ! Sigma
       distance = min(d_to_boundary, d_to_collision)

       ! Advance neutron
       neut%xyz = neut%xyz + distance*neut%uvw

       ! Add pathlength tallies

       if (d_to_collision > d_to_boundary) then
          neut%surface = surf
          neut%cell = 0
          call cross_boundary(neut)
       else
          ! collision
          call collision(neut)
       end if
       
    end do

  end subroutine transport

!=====================================================================
! COLLISION
!=====================================================================

  subroutine collision(neut)

    type(Neutron), pointer :: neut

    type(AceContinuous), pointer :: table
    type(AceReaction),   pointer :: rxn
    real(8) :: r1
    character(250) :: msg
    integer :: i,j
    integer :: n_isotopes
    integer :: IE
    real(8) :: f, Sigma, total
    real(8) :: density, density_i
    real(8) :: p
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
       IE = table%grid_index(neut%IE)
       f = (neut%E - table%energy(IE))/(table%energy(IE+1) - table%energy(IE))

       Sigma = density_i*(f*table%sigma_t(IE) + (1-f)*(table%sigma_t(IE+1)))
       Sigma_t(i) = Sigma
    end do

    ! sample nuclide
    r1 = rang()
    p = 0.0_8
    total = sum(Sigma_t)
    do i = 1, n_isotopes
       p = p + Sigma_t(i) / total
       if (r1 < p) exit
    end do

    ! Get table, total xs, interpolation factor
    table => xs_continuous(cMaterial%table(i))
    Sigma = Sigma_t(i)
    IE = table%grid_index(neut%IE)
    f = (neut%E - table%energy(IE))/(table%energy(IE+1) - table%energy(IE))
    density = cMaterial%atom_percent(i)*density

    ! free memory
    deallocate(Sigma_t)

    ! sample reaction channel
    r1 = rang()*Sigma
    p = 0.0_8
    do i = 1, table%n_reaction
       rxn => table%reactions(i)
       if (rxn%MT >= 200) cycle
       if (IE < rxn%IE) cycle
       p = p + density * (f*rxn%sigma(IE-rxn%IE+1) + (1-f)*(rxn%sigma(IE-rxn%IE+2)))
       if (r1 < p) exit
    end do
    if (verbosity >= 10) then
       msg = "    " // trim(reaction_name(rxn%MT)) // " with nuclide " // &
            & trim(table%name)
       call message(msg, 10)
    end if

    ! call appropriate subroutine
    select case (rxn%MT)
    case (2)
       call elastic_scatter(neut, table%awr)
    case (102)
       call n_gamma(neut)
    case default
       call elastic_scatter(neut, table%awr)
    end select
    
  end subroutine collision

!=====================================================================
! ELASTIC_SCATTER
!=====================================================================

  subroutine elastic_scatter(neut, awr)

    type(Neutron), pointer :: neut
    real(8), intent(in) :: awr

    real(8) :: phi ! azimuthal angle
    real(8) :: mu  ! cosine of polar angle
    real(8) :: vx, vy, vz
    real(8) :: vcx, vcy ,vcz
    real(8) :: vel
    real(8) :: u, v, w
    real(8) :: E
    integer :: IE

    vel = sqrt(neut%E)

    vx = vel*neut%uvw(1)
    vy = vel*neut%uvw(2)
    vz = vel*neut%uvw(3)


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

    neut%E = E
    neut%uvw(1) = vx/vel
    neut%uvw(2) = vy/vel
    neut%uvw(3) = vz/vel

    ! find energy index, interpolation factor
    IE = binary_search(e_grid, n_grid, E)
    neut%IE = IE
    neut%interp = (E - e_grid(IE))/(e_grid(IE+1) - e_grid(IE))

  end subroutine elastic_scatter

!=====================================================================
! LEVEL_INELASTIC
!=====================================================================

  subroutine level_inelastic

  end subroutine level_inelastic

!=====================================================================
! N_GAMMA
!=====================================================================

  subroutine n_gamma(neut)

    type(Neutron), pointer :: neut

    integer :: cell_num
    character(250) :: msg

    neut%alive = .false.
    if (verbosity >= 10) then
       cell_num = cells(neut%cell)%uid
       msg = "    Absorbed in cell " // trim(int_to_str(cell_num))
       call message(msg, 10)
    end if

  end subroutine n_gamma

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
