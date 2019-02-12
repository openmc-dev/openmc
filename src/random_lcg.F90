module random_lcg

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    function prn() result(pseudo_rn) bind(C)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE) :: pseudo_rn
    end function prn

    subroutine set_particle_seed(id) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: id
    end subroutine set_particle_seed

    subroutine prn_set_stream(n) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value :: n
    end subroutine prn_set_stream
  end interface
end module random_lcg
