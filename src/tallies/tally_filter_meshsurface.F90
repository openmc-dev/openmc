module tally_filter_meshsurface

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_cpp

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

  type, extends(CppTallyFilter) :: MeshSurfaceFilter
  end type MeshSurfaceFilter

end module tally_filter_meshsurface
