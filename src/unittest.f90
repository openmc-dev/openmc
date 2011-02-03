program unittest

  use global
  use energy_grid, only: add_grid_points

  call test_energy_grid()

contains

!=====================================================================
! TEST_ENERGY_GRID is a test of the energy grid unionization
! subroutine on a simplified set of data
!=====================================================================

  subroutine test_energy_grid()

    type(LinkedListGrid), pointer :: list => null()
    type(LinkedListGrid), pointer :: current => null()
    real(8), allocatable :: energy(:)
    integer :: n, i
    integer :: n_list

    n = 10
    allocate(energy(n))
    do i = 1,n
       energy(i) = 1.1248_8*i
    end do

    ! create list
    allocate(list)
    current => list
    current%energy = 4.23_8
    allocate(current%next)
    current => current%next
    current%energy = 7.23_8

    call add_grid_points(list, energy, n)

    current => list
    do while(associated(current))
       print *,current%energy
       current => current%next
    end do

    deallocate(energy)

  end subroutine test_energy_grid

end program unittest
