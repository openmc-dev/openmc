module random_lcg

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    function prn() result(pseudo_rn) bind(C)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE) :: pseudo_rn
    end function prn

    function future_prn(n) result(pseudo_rn) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: n
      real(C_DOUBLE)            :: pseudo_rn
    end function future_prn

    subroutine set_particle_seed(id) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: id
    end subroutine set_particle_seed

    subroutine advance_prn_seed(n) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: n
    end subroutine advance_prn_seed

    subroutine prn_set_stream(n) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value :: n
    end subroutine prn_set_stream

    function openmc_get_seed() result(seed) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T) :: seed
    end function openmc_get_seed

    subroutine openmc_set_seed(new_seed) bind(C)
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: new_seed
    end subroutine openmc_set_seed
  end interface
end module random_lcg
