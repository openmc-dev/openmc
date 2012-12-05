module bank_header

  implicit none

!===============================================================================
! BANK is used for storing fission sites in eigenvalue calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type Bank
    ! The 'sequence' attribute is used here to ensure that the data listed
    ! appears in the given order. This is important for MPI purposes when bank
    ! sites are sent from one processor to another.
    sequence

    real(8)    :: wgt    ! weight of bank site
    real(8)    :: xyz(3) ! location of bank particle
    real(8)    :: uvw(3) ! diretional cosines
    real(8)    :: E      ! energy
  end type Bank

end module bank_header
