module energy_grid

  use global

  implicit none

  integer :: grid_method ! how to treat the energy grid
  integer :: n_log_bins  ! number of bins for logarithmic grid
  real(8) :: log_spacing ! spacing on logarithmic grid

contains

!===============================================================================
! UNIONIZED_GRID creates a unionized energy grid, for the entire problem or for
! each material, composed of the grids from each nuclide in the entire problem,
! or each material, respectively.  Right now, the grid for each nuclide is added
! into a linked list one at a time with an effective insertion sort. Could be
! done with a hash for all energy points and then a quicksort at the end (what
! hash function to use?)
!===============================================================================

  subroutine unionized_grid()

    integer :: i ! index in nuclides array
    integer :: j ! index in materials array
    type(ListReal), pointer, save :: list => null()
    type(Nuclide),  pointer, save :: nuc  => null()
    type(Material), pointer, save :: mat  => null()
    !$omp threadprivate(list, nuc, mat)

    call write_message("Creating unionized energy grid...", 5)

    ! add grid points for each nuclide in the material
    do j = 1, n_materials
      mat => materials(j)
      do i = 1, mat % n_nuclides
        nuc => nuclides(mat % nuclide(i))
        call add_grid_points(list, nuc % energy)
      end do

      ! set size of unionized material energy grid
      mat % n_grid = list % size()

      ! create allocated array from linked list
      allocate(mat % e_grid(mat % n_grid))
      do i = 1, mat % n_grid
        mat % e_grid(i) = list % get_item(i)
      end do

      ! delete linked list and dictionary
      call list % clear()
      deallocate(list)
    end do

    ! Set pointers to unionized energy grid for each nuclide
    call grid_pointers()

  end subroutine unionized_grid
>>>>>>> mit-crpg/openmc/develop

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
    M = n_log_bins
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

!===============================================================================
! GRID_POINTERS creates an array of pointers (ints) for each nuclide to link
! each point on the nuclide energy grid to one on a unionized energy grid
!===============================================================================

  subroutine grid_pointers()

    integer :: i            ! loop index for nuclides
    integer :: j            ! loop index for nuclide energy grid
    integer :: k            ! loop index for materials
    integer :: index_e      ! index on union energy grid
    real(8) :: union_energy ! energy on union grid
    real(8) :: energy       ! energy on nuclide grid
    type(Nuclide),  pointer :: nuc => null()
    type(Material), pointer :: mat => null()

    do k = 1, n_materials
      mat => materials(k)
      allocate(mat % nuclide_grid_index(mat % n_nuclides, mat % n_grid))
      do i = 1, mat % n_nuclides
        nuc => nuclides(mat % nuclide(i))

        index_e = 1
        energy = nuc % energy(index_e)

        do j = 1, mat % n_grid
          union_energy = mat % e_grid(j)
          if (union_energy >= energy .and. index_e < nuc % n_grid) then
            index_e = index_e + 1
            energy = nuc % energy(index_e)
          end if
          mat % nuclide_grid_index(i,j) = index_e - 1
        end do
      end do
    end do

  end subroutine grid_pointers

end module energy_grid
