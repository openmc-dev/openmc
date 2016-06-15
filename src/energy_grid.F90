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
    type(ListReal) :: list
    type(Nuclide),  pointer :: nuc
    type(Material), pointer :: mat

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
    end do

    ! Set pointers to unionized energy grid for each nuclide
    call grid_pointers()

  end subroutine unionized_grid

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
    type(Nuclide), pointer :: nuc

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
      ! Allocate logarithmic mapping for nuclide
      nuc => nuclides(i)
      allocate(nuc % grid_index(0:M))

      ! Determine corresponding indices in nuclide grid to energies on
      ! equal-logarithmic grid
      j = 1
      do k = 0, M
        do while (log(nuc%energy(j + 1)/E_min) <= umesh(k))
          ! Ensure that for isotopes where maxval(nuc % energy) << E_max
          ! that there are no out-of-bounds issues.
          if (j + 1 == nuc % n_grid) then
            exit
          end if
          j = j + 1
        end do
        nuc % grid_index(k) = j
      end do
    end do

    deallocate(umesh)

  end subroutine logarithmic_grid

!===============================================================================
! ADD_GRID_POINTS adds energy points from the 'energy' array into a linked list
! of points already stored from previous arrays.
!===============================================================================

  subroutine add_grid_points(list, energy)

    type(ListReal) :: list
    real(8), intent(in) :: energy(:)

    integer :: i       ! index in energy array
    integer :: n       ! size of energy array
    integer :: current ! current index
    real(8) :: E       ! actual energy value

    i = 1
    n = size(energy)

    ! Set current index to beginning of the list
    current = 1

    do while (i <= n)
      E = energy(i)

      ! If we've reached the end of the grid energy list, add the remaining
      ! energy points to the end
      if (current > list % size()) then
        ! Finish remaining energies
        do while (i <= n)
          call list % append(energy(i))
          i = i + 1
        end do
        exit
      end if

      if (E < list % get_item(current)) then

        ! Insert new energy in this position
        call list % insert(current, E)

        ! Advance index in linked list and in new energy grid
        i = i + 1
        current = current + 1

      elseif (E == list % get_item(current)) then
        ! Found the exact same energy, no need to store duplicates so just
        ! skip and move to next index
        i = i + 1
        current = current + 1
      else
        current = current + 1
      end if

    end do

  end subroutine add_grid_points

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
    type(Nuclide),  pointer :: nuc
    type(Material), pointer :: mat

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
