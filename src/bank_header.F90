module bank_header

  implicit none

!===============================================================================
! BANK is used for storing fission sites in criticality calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type Bank
     sequence
     integer(8) :: id     ! Unique ID
     real(8)    :: xyz(3) ! location of bank particle
     real(8)    :: uvw(3) ! diretional cosines
     real(8)    :: E      ! energy
  end type Bank

end module bank_header
