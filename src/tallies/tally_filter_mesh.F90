module tally_filter_mesh

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use mesh_header
  use tally_filter_header

  implicit none
  private
  public :: openmc_mesh_filter_get_mesh
  public :: openmc_mesh_filter_set_mesh

!===============================================================================
! MESHFILTER indexes the location of particle events to a regular mesh.  For
! tracklength tallies, it will produce multiple valid bins and the bin weight
! will correspond to the fraction of the track length that lies in that bin.
!===============================================================================

  type, public, extends(CppTallyFilter) :: MeshFilter
  contains
    procedure :: mesh => get_mesh
  end type MeshFilter

contains

  function get_mesh(this) result(mesh)
    class(MeshFilter) :: this
    integer :: mesh
    interface
      function mesh_filter_get_mesh(filt) result(index_mesh) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT)     :: index_mesh
      end function mesh_filter_get_mesh
    end interface
    mesh = mesh_filter_get_mesh(this % ptr)
  end function get_mesh

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_mesh_filter_get_mesh(index, index_mesh) result(err) bind(C)
    ! Get the mesh for a mesh filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), intent(out)       :: index_mesh
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshFilter)
        index_mesh = f % mesh()
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set mesh on a non-mesh filter.")
      end select
    end if
  end function openmc_mesh_filter_get_mesh


  function openmc_mesh_filter_set_mesh(index, index_mesh) result(err) bind(C)
    ! Set the mesh for a mesh filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: index_mesh
    integer(C_INT) :: err

    interface
      subroutine mesh_filter_set_mesh(filt, mesh) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: mesh
      end subroutine mesh_filter_set_mesh
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshFilter)
        if (index_mesh >= 0 .and. index_mesh < n_meshes()) then
          call mesh_filter_set_mesh(f % ptr, index_mesh)
          f % n_bins = f % n_bins_cpp()
        else
          err = E_OUT_OF_BOUNDS
          call set_errmsg("Index in 'meshes' array is out of bounds.")
        end if
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set mesh on a non-mesh filter.")
      end select
    end if
  end function openmc_mesh_filter_set_mesh

end module tally_filter_mesh
