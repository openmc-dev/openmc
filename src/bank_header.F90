module bank_header

  use, intrinsic :: ISO_C_BINDING

  implicit none

!===============================================================================
! BANK is used for storing fission sites in eigenvalue calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type, bind(C) :: Bank
    real(C_DOUBLE) :: wgt           ! weight of bank site
    real(C_DOUBLE) :: xyz(3)        ! location of bank particle
    real(C_DOUBLE) :: uvw(3)        ! diretional cosines
    real(C_DOUBLE) :: E             ! energy / energy group if in MG mode.
    integer(C_INT) :: delayed_group ! delayed group
    integer(C_INT) :: particle      ! particle type (neutron, photon, etc.)
  end type Bank

  interface
    function openmc_fission_bank(ptr, n) result(err) bind(C)
      import C_PTR, C_INT64_T, C_INT
      type(C_PTR), intent(out) :: ptr
      integer(C_INT64_T), intent(out) :: n
      integer(C_INT) :: err
    end function

    function fission_bank_delayed_group(i) result(g) bind(C)
      import C_INT64_T, C_INT
      integer(C_INT64_T), value :: i
      integer(C_INT) :: g
    end function

    function fission_bank_E(i) result(E) bind(C, name='fission_bank_E')
      import C_INT64_T, C_DOUBLE
      integer(C_INT64_T), value :: i
      real(C_DOUBLE) :: E
    end function

    function fission_bank_wgt(i) result(wgt) bind(C)
      import C_INT64_T, C_DOUBLE
      integer(C_INT64_T), value :: i
      real(C_DOUBLE) :: wgt
    end function

    subroutine source_bank_xyz(i, xyz) bind(C)
      import C_INT64_T, C_DOUBLE
      integer(C_INT64_T), value :: i
      real(C_DOUBLE), intent(in) :: xyz(*)
    end subroutine

    function source_bank_E(i) result(E) bind(C, name='source_bank_E')
      import C_INT64_T, C_DOUBLE
      integer(C_INT64_T), value :: i
      real(C_DOUBLE) :: E
    end function

    function source_bank_wgt(i) result(wgt) bind(C)
      import C_INT64_T, C_DOUBLE
      integer(C_INT64_T), value :: i
      real(C_DOUBLE) :: wgt
    end function

    subroutine source_bank_set_wgt(i, wgt) bind(C)
      import C_INT64_T, C_DOUBLE
      integer(C_INT64_T), value :: i
      real(C_DOUBLE), value :: wgt
    end subroutine
  end interface

end module bank_header
