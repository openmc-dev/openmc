module tally_initialize

  use constants
  use global
  use tally_header, only: TallyObject

  implicit none
  private
  public :: configure_tallies
  public :: add_tallies

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
! TallyObject derived type, including stride, matching_bins, and results.
!===============================================================================

  subroutine setup_tally_arrays()

    integer :: i                 ! loop index for tallies
    integer :: j                 ! loop index for filters
    integer :: n                 ! temporary stride
    integer :: max_n_filters = 0 ! maximum number of filters
    type(TallyObject), pointer :: t

    TALLY_LOOP: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Allocate stride and matching_bins arrays
      allocate(t % stride(size(t % filters)))
      max_n_filters = max(max_n_filters, size(t % filters))

      ! The filters are traversed in opposite order so that the last filter has
      ! the shortest stride in memory and the first filter has the largest
      ! stride

      n = 1
      STRIDE: do j = size(t % filters), 1, -1
        t % stride(j) = n
        n = n * t % filters(j) % obj % n_bins
      end do STRIDE

      ! Set total number of filter and scoring bins
      t % total_filter_bins = n
      t % total_score_bins = t % n_score_bins * t % n_nuclide_bins

      ! Allocate results array
      allocate(t % results(3, t % total_score_bins, t % total_filter_bins))
      t % results(:,:,:) = ZERO

    end do TALLY_LOOP

    ! Allocate array for matching filter bins
!$omp parallel
    allocate(matching_bins(max_n_filters))
    allocate(filter_weights(max_n_filters))
!$omp end parallel

  end subroutine setup_tally_arrays

!===============================================================================
! ADD_TALLIES extends the tallies array with a new group of tallies and assigns
! pointers to each group. This is called once for user tallies, once for CMFD
! tallies, etc.
!===============================================================================

  subroutine add_tallies(tally_group, n)

    character(*), intent(in) :: tally_group ! name of tally group
    integer,      intent(in) :: n           ! number of tallies to add

    type(TallyObject), allocatable :: temp(:) ! temporary tallies array

    if (n_tallies == 0) then
      ! Allocate tallies array
      allocate(tallies(n))
    else
      ! Allocate tallies array with increased size
      allocate(temp(n_tallies + n))

      ! Copy original tallies to temporary array
      temp(1:n_tallies) = tallies

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=tallies)

    end if

    ! Set index for ths tally group
    select case(tally_group)
    case ("user")
      i_user_tallies = n_tallies
    case ("cmfd")
      i_cmfd_tallies = n_tallies
    end select

    ! Set n_tallies
    n_tallies = size(tallies)

    ! Reassign pointers for each group -- after the call to move_alloc, any
    ! pointers that were associated with tallies before become unassociated
    if (i_user_tallies >= 0) then
      user_tallies => tallies(i_user_tallies+1 : i_user_tallies+n_user_tallies)
    end if
    if (i_cmfd_tallies >= 0) then
      cmfd_tallies => tallies(i_cmfd_tallies+1 : i_cmfd_tallies+n_cmfd_tallies)
    end if

  end subroutine add_tallies

end module tally_initialize
