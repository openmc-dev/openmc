module URR_edits

  use output, only: write_coords

  implicit none
  private

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

    D_mean = 20.0_8

    do i = 1, n
      D_vals(i) = wigner_surmise(D_mean)
    end do

    call write_coords(99, 'wigner-dist-samples.dat',&
         dble([(i, i = 1, n)]), D_vals(:) / D_mean)

  end subroutine sampleWignerDist

end module URR_edits
