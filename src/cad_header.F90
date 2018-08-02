#ifdef CAD

module cad_header

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    subroutine load_cad_geometry_c() bind(C)
    end subroutine load_cad_geometry_c

    subroutine free_memory_cad_c() bind(C)
    end subroutine free_memory_cad_c

  end interface

contains

  subroutine load_cad_geometry()
    call load_cad_geometry_c()
  end subroutine load_cad_geometry

  subroutine free_memory_cad()
    call free_memory_cad_c()
  end subroutine free_memory_cad

end module cad_header

#endif
