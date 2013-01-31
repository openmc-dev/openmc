module tally_initialize

  use constants
  use global
  use tally_header, only: TallyObject, TallyMapElement, TallyMapItem

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

    call setup_tally_arrays()
    call setup_tally_maps()

  end subroutine configure_tallies

!===============================================================================
! SETUP_TALLY_ARRAYS allocates and populates several member arrays of the
! TallyObject derived type, including stride, matching_bins, and results.
!===============================================================================

  subroutine setup_tally_arrays()

    integer :: i ! loop index for tallies
    integer :: j ! loop index for filters
    integer :: n ! temporary stride
    type(TallyObject), pointer :: t => null()

    TALLY_LOOP: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Allocate stride and matching_bins arrays
      allocate(t % stride(t % n_filters))
      allocate(t % matching_bins(t % n_filters))

      ! The filters are traversed in opposite order so that the last filter has
      ! the shortest stride in memory and the first filter has the largest
      ! stride

      n = 1
      STRIDE: do j = t % n_filters, 1, -1
        t % stride(j) = n
        n = n * t % filters(j) % n_bins
      end do STRIDE

      ! Set total number of filter and scoring bins
      t % total_filter_bins = n
      t % total_score_bins = t % n_score_bins * t % n_nuclide_bins

      ! Allocate results array
      allocate(t % results(t % total_score_bins, t % total_filter_bins))

    end do TALLY_LOOP

  end subroutine setup_tally_arrays

!===============================================================================
! SETUP_TALLY_MAPS creates a map that allows a quick determination of which
! tallies and bins need to be scored to when a particle makes a collision. This
! subroutine also sets the stride attribute for each tally as well as allocating
! storage for the results array.
!===============================================================================

  subroutine setup_tally_maps()

    integer :: i    ! loop index for tallies
    integer :: j    ! loop index for filters
    integer :: k    ! loop index for bins
    integer :: bin  ! filter bin entries
    integer :: type ! type of tally filter
    type(TallyObject), pointer :: t => null()

    ! allocate tally map array -- note that we don't need a tally map for the
    ! energy_in and energy_out filters
    allocate(tally_maps(N_FILTER_TYPES - 3))

    ! allocate list of items for each different filter type
    allocate(tally_maps(FILTER_UNIVERSE) % items(n_universes))
    allocate(tally_maps(FILTER_MATERIAL) % items(n_materials))
    allocate(tally_maps(FILTER_CELL)     % items(n_cells))
    allocate(tally_maps(FILTER_CELLBORN) % items(n_cells))
    allocate(tally_maps(FILTER_SURFACE)  % items(n_surfaces))

    TALLY_LOOP: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! No need to set up tally maps for surface current tallies
      if (t % type == TALLY_SURFACE_CURRENT) cycle

      FILTER_LOOP: do j = 1, t % n_filters
        ! Determine type of filter
        type = t % filters(j) % type

        if (type == FILTER_CELL .or. type == FILTER_SURFACE .or. &
             type == FILTER_MATERIAL .or. type == FILTER_UNIVERSE .or. &
             type == FILTER_CELLBORN) then

          ! Add map elements
          BIN_LOOP: do k = 1, t % filters(j) % n_bins
            bin = t % filters(j) % int_bins(k)
            call add_map_element(tally_maps(type) % items(bin), i, k)
          end do BIN_LOOP
        end if

      end do FILTER_LOOP

    end do TALLY_LOOP

  end subroutine setup_tally_maps

!===============================================================================
! ADD_MAP_ELEMENT adds a pair of tally and bin indices to the list for a given
! cell/surface/etc.
!===============================================================================

  subroutine add_map_element(item, index_tally, index_bin)

    type(TallyMapItem), intent(inout) :: item
    integer, intent(in) :: index_tally ! index in tallies array
    integer, intent(in) :: index_bin   ! index in bins array

    integer :: n                       ! size of elements array
    type(TallyMapElement), allocatable :: temp(:)

    if (.not. allocated(item % elements)) then
      allocate(item % elements(1))
      item % elements(1) % index_tally = index_tally
      item % elements(1) % index_bin   = index_bin
    else
      ! determine size of elements array
      n = size(item % elements)

      ! allocate temporary storage and copy elements
      allocate(temp(n+1))
      temp(1:n) = item % elements

      ! move allocation back to main array
      call move_alloc(FROM=temp, TO=item%elements)

      ! set new element
      item % elements(n+1) % index_tally = index_tally
      item % elements(n+1) % index_bin   = index_bin
    end if

  end subroutine add_map_element

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
