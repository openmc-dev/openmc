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

    subroutine material_set_id_c(mat_ptr, id, index) &
         bind(C, name='material_set_id')
      import C_PTR, C_INT32_T
      type(C_PTR),        intent(in), value :: mat_ptr
      integer(C_INT32_T), intent(in), value :: id
      integer(C_INT32_T), intent(in), value :: index
    end subroutine material_set_id_c

    function material_fissionable_c(mat_ptr) &
         bind(C, name='material_fissionable') result(fissionable)
      import C_PTR, C_BOOL
      type(C_PTR), intent(in), value :: mat_ptr
      logical(C_BOOL)                :: fissionable
    end function material_fissionable_c

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
    character(len=104)   :: name = ""       ! User-defined name
    integer              :: n_nuclides = 0  ! number of nuclides
    integer, allocatable :: nuclide(:)      ! index in nuclides array
    integer, allocatable :: element(:)      ! index in elements array
    real(8)              :: density         ! total atom density in atom/b-cm
    real(C_DOUBLE), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm
    real(8)              :: density_gpcc    ! total density in g/cm^3

    ! To improve performance of tallying, we store an array (direct address
    ! table) that indicates for each nuclide in the global nuclides(:) array the
    ! index of the corresponding nuclide in the Material % nuclide(:) array. If
    ! it is not present in the material, the entry is set to zero.
    integer, allocatable :: mat_nuclide_index(:)

    ! S(a,b) data
    integer              :: n_sab = 0         ! number of S(a,b) tables
    integer, allocatable :: i_sab_nuclides(:) ! index of corresponding nuclide
    integer, allocatable :: i_sab_tables(:)   ! index in sab_tables
    real(8), allocatable :: sab_fracs(:)      ! how often to use S(a,b)

    ! Temporary names read during initialization
    character(20), allocatable :: names(:)     ! isotope names
    character(20), allocatable :: sab_names(:) ! name of S(a,b) table

    ! Does this material contain fissionable nuclides? Is it depletable?
    logical :: depletable = .false.

    ! enforce isotropic scattering in lab for specific nuclides
    logical :: has_isotropic_nuclides = .false.
    logical, allocatable :: p0(:)

  contains
    procedure :: id => material_id
    procedure :: set_id => material_set_id
    procedure :: fissionable => material_fissionable
    procedure :: init_nuclide_index => material_init_nuclide_index
    procedure :: calculate_xs => material_calculate_xs
  end type Material

  integer(C_INT32_T), public, bind(C) :: n_materials ! # of materials

  type(Material), public, allocatable, target :: materials(:)

  ! Dictionary that maps user IDs to indices in 'materials'
  type(DictIntInt), public :: material_dict

contains

  function material_id(this) result(id)
    class(Material), intent(in) :: this
    integer(C_INT32_T)          :: id
    id = material_id_c(this % ptr)
  end function material_id

  subroutine material_set_id(this, id, index)
    class(Material),    intent(in) :: this
    integer(C_INT32_T), intent(in) :: id
    integer(C_INT32_T), intent(in) :: index
    call material_set_id_c(this % ptr, id, index)
  end subroutine material_set_id

  function material_fissionable(this) result(fissionable)
    class(Material), intent(in) :: this
    logical(C_BOOL)             :: fissionable
    fissionable = material_fissionable_c(this % ptr)
  end function material_fissionable

!===============================================================================
! INIT_NUCLIDE_INDEX creates a mapping from indices in the global nuclides(:)
! array to the Material % nuclides array
!===============================================================================

  subroutine material_init_nuclide_index(this)
    class(Material), intent(inout) :: this

    integer :: i

    ! Allocate nuclide index array and set to zeros
    if (allocated(this % mat_nuclide_index)) &
         deallocate(this % mat_nuclide_index)
    allocate(this % mat_nuclide_index(n_nuclides))
    this % mat_nuclide_index(:) = 0

    ! Assign entries in the index array
    do i = 1, this % n_nuclides
      this % mat_nuclide_index(this % nuclide(i)) = i
    end do
  end subroutine material_init_nuclide_index

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
    call material_dict % clear()
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

!===============================================================================
! Fortran compatibility
!===============================================================================

  function material_isotropic(i_material, i_nuc_mat) result(iso) bind(C)
    integer(C_INT), value :: i_material
    integer(C_INT), value :: i_nuc_mat
    logical(C_BOOL) :: iso

    iso = .false.
    associate (mat => materials(i_material))
      if (mat % has_isotropic_nuclides) then
        iso = mat % p0(i_nuc_mat)
      end if
    end associate
  end function

  function material_element(i_material) result(ptr) bind(C)
    integer(C_INT), value :: i_material
    type(C_PTR) :: ptr
    ptr = C_LOC(materials(i_material) % element(1))
  end function

end module material_header
