module tally_filter_meshsurface

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use mesh_header
  use tally_filter_header

  implicit none
  private
  public :: openmc_meshsurface_filter_get_mesh
  public :: openmc_meshsurface_filter_set_mesh

!===============================================================================
! MESHFILTER indexes the location of particle events to a regular mesh.  For
! tracklength tallies, it will produce multiple valid bins and the bin weight
! will correspond to the fraction of the track length that lies in that bin.
!===============================================================================

  type, public, extends(CppTallyFilter) :: MeshSurfaceFilter
  end type MeshSurfaceFilter

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_meshsurface_filter_get_mesh(index, index_mesh) result(err) bind(C)
    ! Get the mesh for a mesh surface filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), intent(out)       :: index_mesh
    integer(C_INT) :: err

    interface
      function meshsurface_filter_get_mesh(filt) result(index_mesh) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT)     :: index_mesh
      end function meshsurface_filter_get_mesh
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshSurfaceFilter)
        index_mesh = meshsurface_filter_get_mesh(f % ptr)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set mesh on a non-mesh filter.")
      end select
    end if
  end function openmc_meshsurface_filter_get_mesh


  function openmc_meshsurface_filter_set_mesh(index, index_mesh) result(err) bind(C)
    ! Set the mesh for a mesh surface filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: index_mesh
    integer(C_INT) :: err

    interface
      subroutine meshsurface_filter_set_mesh(filt, mesh) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: mesh
      end subroutine meshsurface_filter_set_mesh
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshSurfaceFilter)
        if (index_mesh >= 0 .and. index_mesh < n_meshes()) then
          call meshsurface_filter_set_mesh(f % ptr, index_mesh)
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
  end function openmc_meshsurface_filter_set_mesh

end module tally_filter_meshsurface
