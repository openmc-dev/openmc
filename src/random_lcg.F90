module random_lcg

  use, intrinsic :: ISO_C_BINDING

  implicit none

  integer(C_INT64_T), bind(C) :: seed

  interface
    function prn() result(pseudo_rn) bind(C, name='prn')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE) :: pseudo_rn
    end function prn

    function future_prn(n) result(pseudo_rn) bind(C, name='future_prn')
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: n
      real(C_DOUBLE)            :: pseudo_rn
    end function future_prn

    subroutine set_particle_seed(id) bind(C, name='set_particle_seed')
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: id
    end subroutine set_particle_seed

    subroutine advance_prn_seed(n) bind(C, name='advance_prn_seed')
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: n
    end subroutine advance_prn_seed

    subroutine prn_set_stream(n) bind(C, name='prn_set_stream')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value :: n
    end subroutine prn_set_stream

    function openmc_set_seed(new_seed) result(err) &
         bind(C, name='openmc_set_seed')
      use ISO_C_BINDING
      implicit none
      integer(C_INT64_T), value :: new_seed
      integer(C_INT)            :: err
    end function openmc_set_seed
  end interface
end module random_lcg
