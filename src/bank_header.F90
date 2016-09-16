module bank_header

  use, intrinsic :: ISO_C_BINDING
  use constants, only: N_STREAMS

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
    integer(C_LONG_LONG) :: prn_seed(N_STREAMS) ! prn seeds to continue site
    integer(C_INT) :: delayed_group ! delayed group
  end type Bank

end module bank_header
