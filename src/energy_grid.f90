module energy_grid

  use global
  use output, only: message

contains

!=====================================================================
! UNIONIZED_GRID creates a single unionized energy grid combined from
! each nuclide of each material. Right now, the grid for each nuclide
! is added into a linked list one at a time with an effective
! insertion sort. Could be done with a hash for all energy points and
! then a quicksort at the end (what hash function to use?)
!=====================================================================

  subroutine unionized_grid()

    type(LinkedListGrid), pointer :: list => null()
    type(LinkedListGrid), pointer :: current => null()

    type(Material),       pointer :: mat => null()
    type(AceContinuous),  pointer :: table => null()
    type(AceReaction),    pointer :: rxn => null()
    integer :: i, j
    integer :: n
    character(100) :: msg

    msg = "Creating unionized energy grid"
    call message(msg, 5)

    ! loop over all materials
    do i = 1, n_materials
       mat => materials(i)
       
       ! loop over all isotopes
       do j = 1, mat%n_isotopes
          table => xs_continuous(mat%table(j))

          ! loop over energy points
          n = size(table%energy)
          call add_grid_points(list, table%energy, n)
       end do
    end do

  end subroutine unionized_grid

!=====================================================================
! ADD_GRID_POINTS adds energy points from the 'energy' array into a
! linked list of points already stored from previous arrays. 
!=====================================================================

  subroutine add_grid_points(list, energy, n_energy)

    type(LinkedListGrid), pointer :: list
    real(8), intent(in) :: energy(n_energy)
    integer, intent(in) :: n_energy

    type(LinkedListGrid), pointer :: current => null()
    type(LinkedListGrid), pointer :: previous => null()
    type(LinkedListGrid), pointer :: head => null()
    type(LinkedListGrid), pointer :: tmp => null()
    integer :: index
    real(8) :: E

    index = 1

    ! if the original list is empty, we need to allocate the first
    ! element and store first energy point
    if (grid_count(list) == 0) then
       allocate(list)
       current => list
       do index = 1, n_energy
          current%energy = energy(index)
          if (index == n_energy) then
             current%next => null()
             return
          end if
          allocate(current%next)
          current => current%next
       end do
    end if

    current => list
    head => list

    E = energy(index)
    do while (index <= n_energy)

       ! If we've reached the end of the grid energy list, add the
       ! remaining energy points to the end
       if (.not. associated(current)) then
          ! finish remaining energies
          do while (index <= n_energy)
             allocate(previous%next)
             current => previous%next
             current%energy = energy(index)
             previous => current
             index = index + 1
          end do
          current%next => null()
          exit
       end if
       
       if (E < current%energy) then
          ! create new element and insert it in energy grid list
          allocate(tmp)
          tmp%energy = E
          tmp%next => current
          if (associated(previous)) then
             previous%next => tmp
             previous => tmp
          else
             previous => tmp
             head => previous
          end if
          nullify(tmp)

          ! advance index
          index = index + 1
          E = energy(index)

       elseif (E == current%energy) then
          ! found the exact same energy, no need to store duplicates
          ! so just skip and move to next index
          index = index + 1
          E = energy(index)
       else
          previous => current
          current => current%next
       end if
       
    end do

    ! It's possible that an element was inserted at the front of the
    ! list, so we need to move the list pointer back to the start of
    ! the list
    list => head

  end subroutine add_grid_points

!=====================================================================
! GRID_COUNT gives the number of energy points in list. This should
! eventually be replaced with a generic list function.
!=====================================================================

  function grid_count(list) result(count)

    type(LinkedListGrid), pointer :: list
    integer :: count

    type(LinkedListGrid), pointer :: current

    ! determine size of list
    if (associated(list)) then
       count = 1
       current => list
       do while (associated(current%next))
          current => current%next
          count = count + 1
       enddo
    else
       count = 0
    end if

  end function grid_count

end module energy_grid
