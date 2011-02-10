module physics

  use global
  use geometry,    only: find_cell, dist_to_boundary, cross_boundary
  use types,       only: Neutron
  use mcnp_random, only: rang
  use output,      only: error, message
  use search,      only: binary_search

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
    real(8) :: r1
    real(8) :: phi ! azimuthal angle
    real(8) :: mu  ! cosine of polar angle
    character(250) :: msg
    integer :: cell_num
    integer :: i,j
    integer :: n_isotopes
    integer :: IE

    real(8) :: f, Sigma
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

    ! normalize to create a discrete pdf
    Sigma_t = Sigma_t / sum(Sigma_t)

    ! sample nuclide
    r1 = rang()
    p = 0.0_8
    do i = 1, n_isotopes
       p = p + Sigma_t(i)
       if (r1 < p) exit
    end do
    table => xs_continuous(cMaterial%table(i))
    ! print *, 'sampled nuclide ', i

    ! select collision type
    r1 = rang()
    if (r1 <= 0.5) then
       ! scatter
       phi = 2.*pi*rang()
       mu = 2.*rang() - 1
       neut%uvw(1) = mu
       neut%uvw(2) = sqrt(1. - mu**2) * cos(phi)
       neut%uvw(3) = sqrt(1. - mu**2) * sin(phi)
    else
       neut%alive = .false.
       if (verbosity >= 10) then
          cell_num = cells(neut%cell)%uid
          msg = "    Absorbed in cell " // trim(int_to_str(cell_num))
          call message(msg, 10)
       end if
       return
    end if
    
  end subroutine collision

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
