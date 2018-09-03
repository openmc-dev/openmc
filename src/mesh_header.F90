module mesh_header

  use, intrinsic :: ISO_C_BINDING

  implicit none

!===============================================================================
! STRUCTUREDMESH represents a tessellation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type :: RegularMesh
    type(C_PTR) :: ptr
  contains
    procedure :: id => regular_id
    procedure :: volume_frac => regular_volume_frac
    procedure :: n_dimension => regular_n_dimension
    procedure :: dimension => regular_dimension
    procedure :: lower_left => regular_lower_left
    procedure :: upper_right => regular_upper_right
    procedure :: width => regular_width

    procedure :: get_bin => regular_get_bin
    procedure :: get_indices => regular_get_indices
    procedure :: get_bin_from_indices => regular_get_bin_from_indices
    procedure :: get_indices_from_bin => regular_get_indices_from_bin
  end type RegularMesh

  interface
    function openmc_extend_meshes(n, index_start, index_end) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value, intent(in) :: n
      integer(C_INT32_T), optional, intent(out) :: index_start
      integer(C_INT32_T), optional, intent(out) :: index_end
      integer(C_INT) :: err
    end function openmc_extend_meshes

    function openmc_get_mesh_index(id, index) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value :: id
      integer(C_INT32_T), intent(out) :: index
      integer(C_INT) :: err
    end function openmc_get_mesh_index

    function openmc_mesh_get_id(index, id) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value       :: index
      integer(C_INT32_T), intent(out) :: id
      integer(C_INT) :: err
    end function openmc_mesh_get_id

    function openmc_mesh_set_id(index, id) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      integer(C_INT32_T), value, intent(in) :: id
      integer(C_INT) :: err
    end function openmc_mesh_set_id

    function openmc_mesh_get_dimension(index, dims, n) result(err) bind(C)
      import C_INT32_T, C_PTR, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      type(C_PTR),               intent(out) :: dims
      integer(C_INT),            intent(out) :: n
      integer(C_INT) :: err
    end function openmc_mesh_get_dimension

    function openmc_mesh_set_dimension(index, n, dims) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      integer(C_INT),     value, intent(in) :: n
      integer(C_INT),            intent(in) :: dims(n)
      integer(C_INT) :: err
    end function openmc_mesh_set_dimension

    function openmc_mesh_get_params(index, ll, ur, width, n) result(err) bind(C)
      import C_INT32_T, C_PTR, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      type(C_PTR), intent(out) :: ll
      type(C_PTR), intent(out) :: ur
      type(C_PTR), intent(out) :: width
      integer(C_INT), intent(out) :: n
      integer(C_INT) :: err
    end function openmc_mesh_get_params

    function openmc_mesh_set_params(index, n, ll, ur, width) result(err) bind(C)
      import C_INT32_T, C_INT, C_DOUBLE
      integer(C_INT32_T), value, intent(in) :: index
      integer(C_INT),     value, intent(in) :: n
      real(C_DOUBLE), intent(in), optional :: ll(n)
      real(C_DOUBLE), intent(in), optional :: ur(n)
      real(C_DOUBLE), intent(in), optional :: width(n)
      integer(C_INT) :: err
    end function openmc_mesh_set_params

    function mesh_id(ptr) result(id) bind(C)
      import C_PTR, C_INT32_T
      type(C_PTR), value :: ptr
      integer(C_INT32_T) :: id
    end function

    function mesh_volume_frac(ptr) result(volume_frac) bind(C)
      import C_PTR, C_DOUBLE
      type(C_PTR), value :: ptr
      real(C_DOUBLE) :: volume_frac
    end function

    function mesh_n_dimension(ptr) result(n) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: ptr
      integer(C_INT) :: n
    end function

    function mesh_dimension(ptr, i) result(d) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: ptr
      integer(C_INT), value :: i
      integer(C_INT) :: d
    end function

    function mesh_lower_left(ptr, i) result(ll) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: ptr
      integer(C_INT), value :: i
      real(C_DOUBLE) :: ll
    end function

    function mesh_upper_right(ptr, i) result(ur) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: ptr
      integer(C_INT), value :: i
      real(C_DOUBLE) :: ur
    end function

    function mesh_width(ptr, i) result(w) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: ptr
      integer(C_INT), value :: i
      real(C_DOUBLE) :: w
    end function

    pure function mesh_get_bin(ptr, xyz) result(bin) bind(C)
      import C_PTR, C_DOUBLE, C_INT
      type(C_PTR), value :: ptr
      real(C_DOUBLE), intent(in) :: xyz(*)
      integer(C_INT) :: bin
    end function

    pure function mesh_get_bin_from_indices(ptr, ijk) result(bin) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: ptr
      integer(C_INT), intent(in) :: ijk(*)
      integer(C_INT) :: bin
    end function

    pure subroutine mesh_get_indices(ptr, xyz, ijk, in_mesh) bind(C)
      import C_PTR, C_DOUBLE, C_INT, C_BOOL
      type(C_PTR), value :: ptr
      real(C_DOUBLE), intent(in) :: xyz(*)
      integer(C_INT), intent(out) :: ijk(*)
      logical(C_BOOL), intent(out) :: in_mesh
    end subroutine

    pure subroutine mesh_get_indices_from_bin(ptr, bin, ijk) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: ptr
      integer(C_INT), value :: bin
      integer(C_INT), intent(out) :: ijk(*)
    end subroutine

    function mesh_ptr(i) result(ptr) bind(C)
      import C_INT, C_PTR
      integer(C_INT), value :: i
      type(C_PTR) :: ptr
    end function

    subroutine read_meshes(node_ptr) bind(C)
      import C_PTR
      type(C_PTR) :: node_ptr
    end subroutine

    function n_meshes() result(n) bind(C)
      import C_INT
      integer(C_INT) :: n
    end function
  end interface

contains

  function meshes(i) result(m)
    integer, intent(in) :: i
    type(RegularMesh) :: m

    m % ptr = mesh_ptr(i)
  end function

  function regular_id(this) result(id)
    class(RegularMesh), intent(in) :: this
    integer(C_INT32_T) :: id
    id = mesh_id(this % ptr)
  end function

  function regular_volume_frac(this) result(volume_frac)
    class(RegularMesh), intent(in) :: this
    real(C_DOUBLE) :: volume_frac
    volume_frac = mesh_volume_frac(this % ptr)
  end function

  function regular_n_dimension(this) result(n)
    class(RegularMesh), intent(in) :: this
    integer(C_INT) :: n
    n = mesh_n_dimension(this % ptr)
  end function

  function regular_dimension(this, i) result(d)
    class(RegularMesh), intent(in) :: this
    integer(C_INT), intent(in) :: i
    integer(C_INT) :: d
    d = mesh_dimension(this % ptr, i)
  end function

  function regular_lower_left(this, i) result(ll)
    class(RegularMesh), intent(in) :: this
    integer(C_INT), intent(in) :: i
    real(C_DOUBLE) :: ll
    ll = mesh_lower_left(this % ptr, i)
  end function

  function regular_upper_right(this, i) result(ur)
    class(RegularMesh), intent(in) :: this
    integer(C_INT), intent(in) :: i
    real(C_DOUBLE) :: ur
    ur = mesh_upper_right(this % ptr, i)
  end function

  function regular_width(this, i) result(w)
    class(RegularMesh), intent(in) :: this
    integer(C_INT), intent(in) :: i
    real(C_DOUBLE) :: w
    w = mesh_width(this % ptr, i)
  end function

!===============================================================================
! GET_MESH_BIN determines the tally bin for a particle in a structured mesh
!===============================================================================

  pure subroutine regular_get_bin(this, xyz, bin)
    class(RegularMesh), intent(in) :: this
    real(8), intent(in)           :: xyz(:) ! coordinates
    integer, intent(out)          :: bin    ! tally bin

    bin = mesh_get_bin(this % ptr, xyz)
  end subroutine

!===============================================================================
! GET_MESH_INDICES determines the indices of a particle in a structured mesh
!===============================================================================

  pure subroutine regular_get_indices(this, xyz, ijk, in_mesh)
    class(RegularMesh), intent(in) :: this
    real(8), intent(in)           :: xyz(:)  ! coordinates to check
    integer, intent(out)          :: ijk(:)  ! indices in mesh
    logical, intent(out)          :: in_mesh ! were given coords in mesh?

    logical(C_BOOL) :: in_mesh_

    call mesh_get_indices(this % ptr, xyz, ijk, in_mesh_)
    in_mesh = in_mesh_
  end subroutine regular_get_indices

!===============================================================================
! MESH_INDICES_TO_BIN maps (i), (i,j), or (i,j,k) indices to a single bin number
! for use in a TallyObject results array
!===============================================================================

  pure function regular_get_bin_from_indices(this, ijk) result(bin)
    class(RegularMesh), intent(in) :: this
    integer, intent(in)           :: ijk(:)
    integer                       :: bin

    bin = mesh_get_bin_from_indices(this % ptr, ijk)
  end function regular_get_bin_from_indices

!===============================================================================
! BIN_TO_MESH_INDICES maps a single mesh bin from a TallyObject results array to
! (i), (i,j), or (i,j,k) indices
!===============================================================================

  pure subroutine regular_get_indices_from_bin(this, bin, ijk)
    class(RegularMesh), intent(in) :: this
    integer, intent(in)           :: bin
    integer, intent(out)          :: ijk(:)

    call mesh_get_indices_from_bin(this % ptr, bin, ijk)
  end subroutine regular_get_indices_from_bin

end module mesh_header
