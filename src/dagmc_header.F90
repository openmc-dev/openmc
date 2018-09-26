#ifdef DAGMC

module dagmc_header

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    subroutine load_dagmc_geometry() bind(C)
    end subroutine load_dagmc_geometry

    subroutine free_memory_dagmc() bind(C)
    end subroutine free_memory_dagmc

  end interface

end module dagmc_header

#endif
