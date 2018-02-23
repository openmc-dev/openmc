module cad_header

  use, intrinsic :: ISO_C_BINDING
  use hdf5

  implicit none

  interface
     subroutine load_cad_geometry_c() bind(C)
     end subroutine load_cad_geometry_c
  end interface
     
contains

  subroutine load_cad_geometry()
    call load_cad_geometry_c()
  end subroutine load_cad_geometry
  
end module cad_header
