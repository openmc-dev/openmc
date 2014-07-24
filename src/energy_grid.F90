module energy_grid

  use constants,        only: MAX_LINE_LEN, N_LOG_BINS
  use global
  use list_header,      only: ListReal
  use output,           only: write_message

contains

!===============================================================================
! UNIONIZED_GRID creates a single unionized energy grid combined from each
! nuclide of each material. Right now, the grid for each nuclide is added into a
! linked list one at a time with an effective insertion sort. Could be done with
! a hash for all energy points and then a quicksort at the end (what hash
! function to use?)
!===============================================================================

  subroutine unionized_grid()

    integer :: i ! index in nuclides array
    type(ListReal), pointer :: list => null()
    type(Nuclide),  pointer :: nuc => null()

    message = "Creating unionized energy grid..."
    call write_message(5)

    ! Add grid points for each nuclide in the problem
    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      call add_grid_points(list, nuc % energy)
    end do

    ! Set size of unionized energy grid 
    n_grid = list % size() 

    ! create allocated array from linked list
    allocate(e_grid(n_grid))
    do i = 1, n_grid
      e_grid(i) = list % get_item(i)
    end do

    ! delete linked list and dictionary
    call list % clear()
    deallocate(list)

    ! Set pointers to unionized energy grid for each nuclide
    call grid_pointers()

  end subroutine unionized_grid

!===============================================================================
! ADD_GRID_POINTS adds energy points from the 'energy' array into a linked list
! of points already stored from previous arrays.
!===============================================================================

  subroutine add_grid_points(list, energy)

    type(ListReal), pointer :: list
    real(8), intent(in) :: energy(:)

    integer :: i       ! index in energy array
    integer :: n       ! size of energy array
    integer :: current ! current index 
    real(8) :: E       ! actual energy value

    i = 1
    n = size(energy)

    ! If the original list is empty, we need to allocate the first element and
    ! store first energy point
    if (.not. associated(list)) then
      allocate(list)
      do i = 1, n
        call list % append(energy(i))
      end do
      return
    end if

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
! each point on the nuclide energy grid to one on the unionized energy grid
!===============================================================================

  subroutine grid_pointers()

    integer :: i            ! loop index for nuclides
    integer :: j            ! loop index for nuclide energy grid
    integer :: index_e      ! index on union energy grid
    real(8) :: union_energy ! energy on union grid
    real(8) :: energy       ! energy on nuclide grid
    type(Nuclide), pointer :: nuc => null()

    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      allocate(nuc % grid_index(n_grid))

      index_e = 1
      energy = nuc % energy(index_e)

      do j = 1, n_grid
        union_energy = e_grid(j)
        if (union_energy >= energy .and. index_e < nuc % n_grid) then
          index_e = index_e + 1
          energy = nuc % energy(index_e)
        end if
        nuc % grid_index(j) = index_e - 1
      end do
    end do

  end subroutine grid_pointers

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
