program test_mesh

  use mesh, only: mesh_indices_to_bin, bin_to_mesh_indices
  use mesh_header, only: StructuredMesh

  implicit none

  integer              :: i, j, k
  integer              :: x, y, z
  integer              :: ijk(3), ijk_new(3)
  integer              :: bin, bin_new
  integer, parameter   :: n = 3
  real(8)              :: xyz(3)
  logical              :: passed
  type(StructuredMesh), pointer :: m => null()

  x = 21
  y = 16
  z = 52

  allocate(m)
  m % n_dimension = n
  allocate(m % dimension(n))
  m % dimension = (/ x, y, z /)

  passed = .true.
  write(*,'(A)', ADVANCE='no') "Testing ijk --> bin:"
  do i = 1, x
     do j = 1, y
        do k = 1, z
           ijk = (/ i, j, k /)
           bin = mesh_indices_to_bin(m, ijk)
           call bin_to_mesh_indices(m, bin, ijk_new)

           ! Check to make sure new indices match old
           if (any(ijk /= ijk_new)) passed = .false.
        end do
     end do
  end do
  call write_result(passed)

  passed = .true.
  write(*,'(A)', ADVANCE='no') "Testing bin --> ijk:"
  do bin = 1, product(m % dimension)
     call bin_to_mesh_indices(m, bin, ijk)
     bin_new = mesh_indices_to_bin(m, ijk)

     ! Check to make sure new bin matches old
     if (bin_new /= bin) passed = .false.
  end do
  call write_result(passed)

contains

  subroutine write_result(passed)

    logical, intent(in) :: passed

    if (passed) then
       write (*,'(5X,A)') char(27) // '[32mPASSED' // char(27) // '[30m'
    else
       write (*,'(5X,A)') char(27) // '[31mFAILED' // char(27) // '[30m'
    end if

  end subroutine write_result

end program test_mesh
