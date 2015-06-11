module edits

  use output,     only: write_coords
  use unresolved, only: wigner_dist

  implicit none

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SAMPLEWIGNERDIST samples the Wigner distribution for level spacings using
! a constant mean spacing in order to create a curve for comparison with the
! analytical result
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine sampleWignerDist(n)

    integer :: n ! number of samples
    integer :: i ! iteration index
    real(8) :: D_mean    ! average level spacing in eV
    real(8) :: D_vals(n) ! sampled level spacings
    character(80) :: temp_str

    D_mean = 20.0_8

    do i = 1, n
      D_vals(i) = wigner_dist(D_mean)
    end do

    temp_str = 'wigner-dist-samples.dat'
    call write_coords(99, temp_str, n, n, dble([(i,i=1,n)]),&
      D_vals(:)/D_mean)

  end subroutine sampleWignerDist

end module edits
