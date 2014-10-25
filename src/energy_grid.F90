module energy_grid

  use constants, only: N_LOG_BINS
  use global

  implicit none

contains

!===============================================================================
! LOGARITHMIC_GRID determines a logarithmic mapping for energies to bounding
! indices on a nuclide energy grid
!===============================================================================

  subroutine logarithmic_grid()

    integer :: i, j, k               ! Loop indices
    integer :: M                     ! Number of equally log-spaced bins
    real(8) :: E_max                 ! Maximum energy in MeV
    real(8) :: E_min                 ! Minimum energy in MeV
    real(8), allocatable :: umesh(:) ! Equally log-spaced energy grid
    type(Nuclide), pointer :: nuc => null()

    ! Set minimum/maximum energies
    E_max = 20.0_8
    E_min = 1.0e-11_8

    ! Determine equal-logarithmic energy spacing
    M = N_LOG_BINS
    log_spacing = log(E_max/E_min)/M

    ! Create equally log-spaced energy grid
    allocate(umesh(0:M))
    umesh(:) = [(i*log_spacing, i=0, M)]

    do i = 1, n_nuclides_total
      ! Allocate logarithmic mapping for nuclide
      nuc => nuclides(i)
      allocate(nuc % grid_index(0:M))

      ! Determine corresponding indices in nuclide grid to energies on
      ! equal-logarithmic grid
      j = 1
      do k = 0, M - 1
        do while (log(nuc%energy(j + 1)/E_min) <= umesh(k))
          j = j + 1
        end do
        nuc % grid_index(k) = j
      end do

      ! Set the last point explicitly so that we don't have out-of-bounds issues
      nuc % grid_index(M) = size(nuc % energy) - 1
    end do

    deallocate(umesh)

  end subroutine logarithmic_grid

end module energy_grid
