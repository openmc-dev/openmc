#ifdef DAGMC

module dagmc_header

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    subroutine load_dagmc_geometry_c() bind(C)
    end subroutine load_dagmc_geometry_c

    subroutine free_memory_dagmc_c() bind(C)
    end subroutine free_memory_dagmc_c

  end interface

contains

  subroutine load_dagmc_geometry()
    call load_dagmc_geometry_c()
  end subroutine load_dagmc_geometry

  subroutine free_memory_dagmc()
    call free_memory_dagmc_c()
  end subroutine free_memory_dagmc

end module dagmc_header

#endif
