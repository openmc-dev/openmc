module energy_grid

  use global
  use list_header, only: ListReal
  use output,      only: write_message


  implicit none

  integer :: grid_method ! how to treat the energy grid
  integer :: n_log_bins  ! number of bins for logarithmic grid
  real(8) :: log_spacing ! spacing on logarithmic grid

contains

!===============================================================================
! LOGARITHMIC_GRID determines a logarithmic mapping for energies to bounding
! indices on a nuclide energy grid
!===============================================================================

  subroutine logarithmic_grid()
    integer :: i, j, k               ! Loop indices
    integer :: t                     ! temperature index
    integer :: M                     ! Number of equally log-spaced bins
    real(8) :: E_max                 ! Maximum energy in MeV
    real(8) :: E_min                 ! Minimum energy in MeV
    real(8), allocatable :: umesh(:) ! Equally log-spaced energy grid

    ! Set minimum/maximum energies
    E_max = energy_max_neutron
    E_min = energy_min_neutron

    ! Determine equal-logarithmic energy spacing
    M = n_log_bins
    log_spacing = log(E_max/E_min)/M

    ! Create equally log-spaced energy grid
    allocate(umesh(0:M))
    umesh(:) = [(i*log_spacing, i=0, M)]

    do i = 1, n_nuclides_total
      associate (nuc => nuclides(i))
        do t = 1, size(nuc % grid)
          ! Allocate logarithmic mapping for nuclide
          allocate(nuc % grid(t) % grid_index(0:M))

          ! Determine corresponding indices in nuclide grid to energies on
          ! equal-logarithmic grid
          j = 1
          do k = 0, M
            do while (log(nuc % grid(t) % energy(j + 1)/E_min) <= umesh(k))
              ! Ensure that for isotopes where maxval(nuc % energy) << E_max
              ! that there are no out-of-bounds issues.
              if (j + 1 == size(nuc % grid(t) % energy)) exit
              j = j + 1
            end do
            nuc % grid(t) % grid_index(k) = j
          end do
        end do
      end associate
    end do

    deallocate(umesh)

  end subroutine logarithmic_grid

end module energy_grid
