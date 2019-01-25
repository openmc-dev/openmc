module material_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header, only: DictIntInt
  use error
  use nuclide_header
  use particle_header, only: Particle
  use photon_header
  use sab_header
  use stl_vector, only: VectorReal, VectorInt
  use string, only: to_str, to_f_string, to_c_string

  implicit none

  private
  public :: free_memory_material
  public :: material_nuclide
  public :: material_nuclide_size
  public :: material_nuclide_index
  public :: material_atom_density
  public :: material_density_gpcc
  public :: openmc_extend_materials
  public :: openmc_material_get_volume
  public :: material_pointer

  interface
    function material_pointer(mat_ind) bind(C) result(ptr)
      import C_PTR, C_INT32_T
      integer(C_INT32_T), intent(in), value :: mat_ind
      type(C_PTR)                           :: ptr
    end function material_pointer

    function material_id_c(mat_ptr) bind(C, name='material_id') result(id)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: mat_ptr
      integer(C_INT32_T)             :: id
    end function material_id_c

    function material_nuclide(i_mat, idx) bind(C) result(nuc)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value :: i_mat
      integer(C_INT), value :: idx
      integer(C_INT) :: nuc
    end function

    function material_nuclide_size(i_mat) bind(C) result(n)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value :: i_mat
      integer(C_INT) :: n
    end function

    function material_nuclide_index(i_mat, i_nuc) bind(C) result(idx)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value :: i_mat
      integer(C_INT), value :: i_nuc
      integer(C_INT) :: idx
    end function

    function material_atom_density(i_mat, idx) bind(C) result(density)
      import C_INT32_T, C_INT, C_DOUBLE
      integer(C_INT32_T), value :: i_mat
      integer(C_INT), value :: idx
      real(C_DOUBLE) :: density
    end function

    function material_density_gpcc(i_mat) bind(C) result(density)
      import C_INT32_T, C_DOUBLE
      integer(C_INT32_T), value :: i_mat
      real(C_DOUBLE) :: density
    end function

    subroutine extend_materials_c(n) bind(C)
      import C_INT32_T
      integer(C_INT32_T), intent(in), value :: n
    end subroutine extend_materials_c

    function openmc_material_get_volume(index, volume) result(err) bind(C)
      import C_INT32_T, C_DOUBLE, C_INT
      integer(C_INT32_T), value :: index
      real(C_DOUBLE), intent(out) :: volume
      integer(C_INT) :: err
    end function openmc_material_get_volume
  end interface

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type, public :: Material
    type(C_PTR) :: ptr
  contains
    procedure :: id => material_id
    procedure :: calculate_xs => material_calculate_xs
  end type Material

  integer(C_INT32_T), public, bind(C) :: n_materials ! # of materials

  type(Material), public, allocatable, target :: materials(:)

contains

  function material_id(this) result(id)
    class(Material), intent(in) :: this
    integer(C_INT32_T)          :: id
    id = material_id_c(this % ptr)
  end function material_id

!===============================================================================
! MATERIAL_CALCULATE_XS determines the macroscopic cross sections for the
! material the particle is currently traveling through.
!===============================================================================

  subroutine material_calculate_xs(this, p)
    class(Material), intent(in) :: this
    type(Particle),  intent(in) :: p

    interface
      subroutine material_calculate_xs_c(ptr, p) bind(C)
        import C_PTR, Particle
        type(C_PTR), value :: ptr
        type(Particle), intent(in) :: p
      end subroutine
    end interface

    call material_calculate_xs_c(this % ptr, p)
  end subroutine material_calculate_xs

!===============================================================================
! FREE_MEMORY_MATERIAL deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_material()
    interface
      subroutine free_memory_material_c() bind(C)
      end subroutine free_memory_material_c
    end interface
    call free_memory_material_c()
    n_materials = 0
    if (allocated(materials)) deallocate(materials)
  end subroutine free_memory_material

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_materials(n, index_start, index_end) result(err) bind(C)
    ! Extend the materials array by n elements
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), optional, intent(out) :: index_start
    integer(C_INT32_T), optional, intent(out) :: index_end
    integer(C_INT) :: err

    integer :: i
    type(Material), allocatable :: temp(:) ! temporary materials array

    if (n_materials == 0) then
      ! Allocate materials array
      allocate(materials(n))
    else
      ! Allocate materials array with increased size
      allocate(temp(n_materials + n))

      ! Move original materials to temporary array
      temp(1:n_materials) = materials(:)

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=materials)
    end if

    ! Return indices in materials array
    if (present(index_start)) index_start = n_materials + 1
    if (present(index_end)) index_end = n_materials + n
    n_materials = n_materials + n

    ! Extend the C++ materials array and get pointers to the C++ objects
    call extend_materials_c(n)
    do i = n_materials - n, n_materials
      materials(i) % ptr = material_pointer(i - 1)
    end do

    err = 0
  end function openmc_extend_materials

end module material_header
