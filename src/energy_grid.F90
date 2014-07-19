module energy_grid

  use constants,        only: MAX_LINE_LEN
  use global
  use list_header,      only: ListReal
  use output,           only: write_message

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
    type(ListReal), pointer :: list => null()
    type(Nuclide),  pointer :: nuc  => null()
    type(Material), pointer :: mat  => null()

    message = "Creating unionized energy grid(s)..."
    call write_message(5)

    select case(grid_method)
    case(GRID_GLOB_UNION)
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

    case(GRID_MAT_UNION)
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
    end select

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

    select case(grid_method)
    case(GRID_GLOB_UNION)
      do i = 1, n_nuclides_total
        nuc => nuclides(i)
        allocate(nuc % glob_grid_index(n_grid))

        index_e = 1
        energy = nuc % energy(index_e)

        do j = 1, n_grid
          union_energy = e_grid(j)
          if (union_energy >= energy .and. index_e < nuc % n_grid) then
            index_e = index_e + 1
            energy = nuc % energy(index_e)
          end if
          nuc % glob_grid_index(j) = index_e - 1
        end do
      end do

    case(GRID_MAT_UNION)
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
    end select

  end subroutine grid_pointers

end module energy_grid
