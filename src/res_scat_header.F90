module res_scat_header

  implicit none
  private
  
  type, public :: ResScatterer

    character(10) :: name    ! name of nuclide, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    integer       :: listing ! index in xs_listings
    real(8)       :: awr     ! weight of nucleus in neutron masses
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Energy grid information
    integer :: n_grid
    integer, allocatable :: grid_index(:) ! pointers to union grid
    real(8), allocatable :: energy(:)     ! energy values corresponding to xs

    ! CDF of neutron velocity x cross section
    real(8), allocatable :: xs_cdf(:)

    ! Microscopic elastic cross section
    real(8), allocatable :: elastic(:)

    ! lower cutoff energy for resonance scattering
    real(8) :: E_min

    ! upper cutoff energy for resonance scattering
    real(8) :: E_max

    ! target velocity sampling scheme
    character(16) :: scheme

    ! Type-Bound procedures
    contains
      procedure :: clear => res_scatterer_clear ! Deallocates resonant scatterer
    
  end type ResScatterer

contains

!===============================================================================
! RES_SCATTERER_CLEAR resets and deallocates data in ResScatterer.
!===============================================================================

  subroutine res_scatterer_clear(this)

    class(ResScatterer), intent(inout) :: this ! the ResScatterer object to clear

    integer :: i ! Loop counter

    if (allocated(this % grid_index)) &
      deallocate(this % grid_index)

    if (allocated(this % energy)) &
      deallocate(this % elastic, this % xs_cdf, this % energy)

  end subroutine res_scatterer_clear

end module res_scat_header
