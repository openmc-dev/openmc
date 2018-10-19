module tally_filter_meshsurface

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use mesh_header
  use tally_filter_header

  implicit none

  interface
    function openmc_meshsurface_filter_set_mesh(index, index_mesh) result(err) &
         bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      integer(C_INT32_T), value, intent(in) :: index_mesh
      integer(C_INT) :: err
    end function
  end interface

!===============================================================================
! MESHFILTER indexes the location of particle events to a regular mesh.  For
! tracklength tallies, it will produce multiple valid bins and the bin weight
! will correspond to the fraction of the track length that lies in that bin.
!===============================================================================

  type, extends(CppTallyFilter) :: MeshSurfaceFilter
  end type MeshSurfaceFilter

end module tally_filter_meshsurface
