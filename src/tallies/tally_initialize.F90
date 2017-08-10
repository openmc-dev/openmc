module tally_initialize

  use constants
  use global
  use tally_header, only: TallyObject

  implicit none
  private
  public :: configure_tallies

contains

!===============================================================================
! CONFIGURE_TALLIES initializes several data structures related to tallies. This
! is called after the basic tally data has already been read from the
! tallies.xml file.
!===============================================================================

  subroutine configure_tallies()

    ! Allocate global tallies
    allocate(global_tallies(3, N_GLOBAL_TALLIES))
    global_tallies(:,:) = ZERO

    call setup_tally_arrays()

  end subroutine configure_tallies

!===============================================================================
! SETUP_TALLY_ARRAYS allocates and populates several member arrays of the
! TallyObject derived type, including stride, filter_matches, and results.
!===============================================================================

  subroutine setup_tally_arrays()

    integer :: i                 ! loop index for tallies
    integer :: j                 ! loop index for filters
    integer :: n                 ! temporary stride
    integer :: i_filt            ! filter index
    type(TallyObject), pointer :: t

    TALLY_LOOP: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Allocate stride
      allocate(t % stride(size(t % filter)))

      ! The filters are traversed in opposite order so that the last filter has
      ! the shortest stride in memory and the first filter has the largest
      ! stride
      n = 1
      STRIDE: do j = size(t % filter), 1, -1
        i_filt = t % filter(j)
        t % stride(j) = n
        n = n * filters(i_filt) % obj % n_bins
      end do STRIDE

      ! Set total number of filter and scoring bins
      t % total_filter_bins = n
      t % total_score_bins = t % n_score_bins * t % n_nuclide_bins

      ! Allocate results array
      allocate(t % results(3, t % total_score_bins, t % total_filter_bins))
      t % results(:,:,:) = ZERO

    end do TALLY_LOOP

  end subroutine setup_tally_arrays

end module tally_initialize
