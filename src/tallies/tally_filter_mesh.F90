module tally_filter_mesh

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_cpp

  implicit none

  interface
    function openmc_mesh_filter_set_mesh(index, index_mesh) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      integer(C_INT32_T), value, intent(in) :: index_mesh
      integer(C_INT) :: err
    end function
  end interface

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

end module tally_filter_mesh
