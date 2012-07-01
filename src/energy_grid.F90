module energy_grid

  use constants,        only: MAX_LINE_LEN
  use datatypes,        only: list_insert, list_size, list_delete, &
                              dict_create, dict_get_key, dict_has_key, &
                              dict_add_key, dict_delete
  use datatypes_header, only: ListReal, DictionaryCI
  use global
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

    integer                 :: i ! index in nuclides array
    type(ListReal), pointer :: list => null()
    type(ListReal), pointer :: current => null()
    type(Nuclide),  pointer :: nuc => null()

    message = "Creating unionized energy grid..."
    call write_message(5)

    ! Add grid points for each nuclide in the problem
    do i = 1, n_nuclides_total
       nuc => nuclides(i)
       call add_grid_points(list, nuc % energy)
    end do

    ! create allocated array from linked list
    n_grid = list_size(list)
    allocate(e_grid(n_grid))
    current => list
    do i = 1, n_grid
       e_grid(i) = current % data
       current => current % next
    end do

    ! delete linked list and dictionary
    call list_delete(list)

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

    integer :: i  ! index in energy array
    integer :: n  ! size of energy array
    real(8) :: E  ! actual energy value
    type(ListReal), pointer :: current => null()
    type(ListReal), pointer :: previous => null()
    type(ListReal), pointer :: head => null()
    type(ListReal), pointer :: tmp => null()

    i = 1
    n = size(energy)

    ! if the original list is empty, we need to allocate the first element and
    ! store first energy point
    if (.not. associated(list)) then
       allocate(list)
       current => list
       do i = 1, n
          current % data = energy(i)
          if (i == n) then
             current % next => null()
             return
          end if
          allocate(current % next)
          current => current % next
       end do
    end if

    current => list
    head => list

    do while (i <= n)
       E = energy(i)

       ! If we've reached the end of the grid energy list, add the remaining
       ! energy points to the end
       if (.not. associated(current)) then
          ! finish remaining energies
          do while (i <= n)
             allocate(previous % next)
             current => previous % next
             current % data = energy(i)
             previous => current
             i = i + 1
          end do
          current%next => null()
          exit
       end if
       
       if (E < current % data) then
          ! create new element and insert it in energy grid list
          allocate(tmp)
          tmp % data = E
          tmp % next => current
          if (associated(previous)) then
             previous % next => tmp
             previous => tmp
          else
             previous => tmp
             head => previous
          end if
          nullify(tmp)

          ! advance index
          i = i + 1

       elseif (E == current % data) then
          ! found the exact same energy, no need to store duplicates so just
          ! skip and move to next index
          i = i + 1
       else
          previous => current
          current => current % next
       end if
       
    end do

    ! It's possible that an element was inserted at the front of the list, so we
    ! need to move the list pointer back to the start of the list
    list => head

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

end module energy_grid
