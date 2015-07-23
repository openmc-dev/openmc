module ndpp_initialize

  use ace_header,   only: Nuclide, XsListing
  use constants
  use dict_header,  only: DictCharInt
  use error,        only: fatal_error, warning
  use global
  use ndpp_header,  only: Ndpp
  use ndpp_ops,     only: ndpp_read
  use output,       only: write_message
  use search
  use string,       only: ends_with, to_lower, starts_with, to_str
  use xml_interface

  implicit none

  ! ndpp_lib.xml preprocessed data listings and associated data.
  type(XsListing), allocatable, target :: ndpp_listings(:)
  type(DictCharInt)                    :: ndpp_listing_dict

contains

!===============================================================================
! READ_NDPP_LIB reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_ndpp_data()
    type(Nuclide), pointer    :: nuc ! Current working nuclide
    type(SAlphaBeta), pointer :: sab ! Current working SAB table
    type(XsListing), pointer  :: ndpp_listing ! The NDPP data listings
    integer :: i_listing     ! index in ndpp_listings array
    integer :: i_nuclide     ! index in nuclides
    integer :: i_sab         ! index in sab_tables
    integer :: scatt_type    ! Whether or not legendre or tabular data
    logical :: get_scatt     ! Flag for whether or not to get scatt data
    logical :: get_nuscatt   ! Flag for whether or not to get nuscatt data
    logical :: get_chi_t     ! Flag for whether or not to get total chi data
    logical :: get_chi_p     ! Flag for whether or not to get prompt chi data
    logical :: get_chi_d     ! Flag for whether or not to get delayed chi data
    integer :: scatt_order   ! Number of scatter moments requested in tallies
    integer :: ndpp_groups   ! Number of groups in the NDPP library

    ! allocate arrays for NDPP data storage
    allocate(ndpp_nuc_data(n_nuclides_total))
    allocate(ndpp_sab_data(n_sab_tables))

    ! Read the ndpp_lib.xml file
    call read_ndpp_xml(scatt_type, ndpp_groups)

    ! Determine which data is required to be stored as well as the maximum
    ! orders for scatt and nuscatt
    call which_data(scatt_type, get_scatt, get_nuscatt, get_chi_t, get_chi_p, &
                    get_chi_d, scatt_order)

    ! Parse through each nuclide in the model and read in the corresponding
    ! NDPP data.
    ! During this we will make sure the temperatures of ACE and NDPP match,
    ! as well as making sure the NDPP data exists at all.

    do i_nuclide = 1, n_nuclides_total
      nuc => nuclides(i_nuclide)
      i_listing = ndpp_listing_dict % get_key(adjustl(trim(nuc % name)))
      if (.not. ndpp_listing_dict % has_key(adjustl(trim(nuc % name)))) then
        ! Could not find ndpp_lib.xml file
        call fatal_error(trim(nuc % name) // " does not exist in " // &
             "NDPP XML file: '" // trim(ndpp_lib) // "'!")
      end if
      ! Read the NDPP data and also check that the temperatures match
      ndpp_listing => ndpp_listings(i_listing)

      ! display message
      call write_message("Loading NDPP data library: " // &
           ndpp_listing % name, 6)

      call ndpp_read(ndpp_nuc_data(i_nuclide), ndpp_listing, get_scatt, &
                     get_nuscatt, get_chi_t, get_chi_p, get_chi_d,  &
                     scatt_order, .True.)
    end do

    do i_sab = 1, n_sab_tables
      sab => sab_tables(i_sab)
      i_listing = ndpp_listing_dict % get_key(adjustl(trim(sab % name)))
      if (.not. ndpp_listing_dict % has_key(adjustl(trim(sab % name)))) then
        ! Could not find ndpp_lib.xml file
        call fatal_error(trim(sab % name) // " does not exist in " // &
             "NDPP XML file: '" // trim(ndpp_lib) // "'!")
      end if
      ! Read the NDPP data and also check that the temperatures match
      ndpp_listing => ndpp_listings(i_listing)

      ! display message
      call write_message("Loading NDPP data library: " // &
           ndpp_listing % name, 6)

      call ndpp_read(ndpp_sab_data(i_sab), ndpp_listing, get_scatt, &
                     get_nuscatt, get_chi_t, get_chi_p, get_chi_d, &
                     scatt_order, .False.)
    end do

    call remove_temperature_independent_data()

    if (allocated(ndpp_listings)) then
      deallocate(ndpp_listings)
    end if
    call ndpp_listing_dict % clear()

    ! Now allocate ndpp_outgoing, our `scratch` variable to store the combined
    ! and interpolated elastic + inelastic data. This is program global and has
    ! been declared threadprivate.
    allocate(ndpp_outgoing(0: n_threads, scatt_order, ndpp_groups))

  end subroutine read_ndpp_data

!===============================================================================
! REMOVE_TEMPERATURE_INDEPENDENT_DATA sweeps through the NDPP data looking
! for nuclides with multiple temperature data sets. If multiplies are found,
! then only one temperature is kept and the rest are replaced with pointers to
! the one remaining temperature. In this way, the tempeature-independent data
! is not unnecessarily duplicated.
!===============================================================================

  subroutine remove_temperature_independent_data()
    integer :: i_nuclide  ! index in ndpp_nuc_data
    integer :: i_other    ! index in ndpp_nuc_data to use while looking for
                          ! other temperatures of same ZAID
    logical, allocatable :: skip(:) ! Data safe to skip (already done)
    type(Ndpp), pointer  :: data    ! Current data to change
    type(Ndpp), pointer  :: ref     ! Current data to point to

    allocate(skip(n_nuclides_total))
    skip(:) = .False.

    ! Parse through each nuclide to find temperature duplicates
    ! During this we will make sure the temperatures of ACE and NDPP match,
    ! as well as making sure the NDPP data exists at all.

    do i_nuclide = 1, n_nuclides_total
      skip(i_nuclide) = .True.
      ref => ndpp_nuc_data(i_nuclide)
      do i_other = 1, n_nuclides_total
        ! Dont do this one if we already have taken care of it, or it is our
        ! reference
        if (skip(i_other)) then
          cycle
        end if
        data => ndpp_nuc_data(i_other)
        ! Now check to see if zaid is the same as the reference
        if (ref % zaid == data % zaid) then
          ! We have a match so just deallocate this guy's temp-indep. data
          ! and point it to the data of i_nuclide
          if (associated(data % inel)) then
            deallocate(data % inel)
            data % inel => ref % inel
          end if
          if (associated(data % nuinel)) then
            deallocate(data % nuinel)
            data % nuinel => ref % nuinel
          end if
          if (associated(data % inel_Ein)) then
            deallocate(data % inel_Ein)
            data % inel_Ein => ref % inel_Ein
          end if
          if (associated(data % inel_Ein_srch)) then
            deallocate(data % inel_Ein_srch)
            data % inel_Ein_srch => ref % inel_Ein_srch
          end if

          ! Ok, good work team, now mark skip as true so we dont do it again
          skip(i_other) = .True.

          ! Clean up the pointers
          nullify(data)
          nullify(ref)
        end if
      end do
    end do

    deallocate(skip)

  end subroutine remove_temperature_independent_data

!===============================================================================
! WHICH_DATA looks at the requested tallies and determines which data is needed
! from the NDPP libraries so that memory utilization is kept to a minimum.
!===============================================================================

  subroutine which_data(scatt_type, get_scatt, get_nuscatt, get_chi_t, &
                        get_chi_p, get_chi_d, scatt_order)
    integer, intent(in)  :: scatt_type  ! Flag for legendre or tabular data
    logical, intent(out) :: get_scatt   ! Flag for if scatt data is needed
    logical, intent(out) :: get_nuscatt ! Flag for if nuscatt data is needed
    logical, intent(out) :: get_chi_t   ! Flag for if total chi data is needed
    logical, intent(out) :: get_chi_p   ! Flag for if prompt chi data is needed
    logical, intent(out) :: get_chi_d   ! Flag for if delay chi data is needed
    integer, intent(out) :: scatt_order ! Number of scatter moments requested

    type(TallyObject), pointer :: t
    integer :: i ! Tally index
    integer :: j ! Score bin index
    integer :: k ! User score bin index

    ! Initialize the flags and orders
    get_scatt = .false.
    get_nuscatt = .false.
    get_chi_t = .false.
    get_chi_p = .false.
    get_chi_d = .false.
    scatt_order = 0

    ! Step through each tally and score and determine which types are present
    ! and the orders for scatt and nuscatt
    TALLY_LOOP: do i = 1, n_tallies
      t => tallies(i)
      j = 0
      SCORE_LOOP: do k = 1, t % n_user_score_bins
        j = j + 1
        select case (t % score_bins(j))
        case (SCORE_NDPP_SCATT_N, SCORE_NDPP_SCATT_PN, SCORE_NDPP_SCATT_YN)
          get_scatt = .true.
          if (t % moment_order(j) > scatt_order) then
            scatt_order = t % moment_order(j)
          end if

        case (SCORE_NDPP_NU_SCATT_N, SCORE_NDPP_NU_SCATT_PN, SCORE_NDPP_NU_SCATT_YN)
          get_nuscatt = .true.
          if (t % moment_order(j) > scatt_order) then
            scatt_order = t % moment_order(j)
          end if

        case (SCORE_NDPP_CHI)
          get_chi_t = .true.

        case (SCORE_NDPP_CHI_P)
          get_chi_p = .true.

        case (SCORE_NDPP_CHI_D)
          get_chi_d = .true.

        end select
      end do SCORE_LOOP
    end do TALLY_LOOP

    ! Adjust the order so array sizing is correct
    if (scatt_type == SCATT_TYPE_LEGENDRE) then
      scatt_order = scatt_order + 1
    end if

  end subroutine which_data

!===============================================================================
! READ_NDPP_XML reads information from a ndpp_lib.xml file. This
! file contains a listing of the NDPP-processed data which may be used.
!===============================================================================

  subroutine read_ndpp_xml(scatt_type, ndpp_groups)

    integer, intent(out)  :: scatt_type  ! Legendre or tabular data
    integer, intent(out)  :: ndpp_groups ! number of groups in NDPP data

    integer :: i, j, k     ! loop indices
    logical :: file_exists ! does ndpp_lib.xml exist?
    integer :: filetype    ! default file type
    logical :: nuscatter   ! is nuscatter data present?
    logical :: chi_present ! is chi data present?
    character(MAX_WORD_LEN)  :: directory   ! directory with cross sections
    character(MAX_LINE_LEN)  :: temp_str
    type(TallyObject), pointer :: t
    integer :: i_filter    ! Index of filter containing energyin or energyout
    ! We can use the same XSListing type for our ndpp data since the NDPP
    ! data is a subset of whats in cross_sections.xml
    type(XsListing), pointer :: listing
    type(Node), pointer      :: doc
    type(Node), pointer      :: node_ndpp
    type(NodeList), pointer  :: node_ndpp_list
    real(8), allocatable     :: ndpp_energy_bins(:)
    integer                  :: order       ! ndpp_lib.xml's scattering order

    ! Check if ndpp_lib.xml exists
    inquire(FILE=ndpp_lib, EXIST=file_exists)
    if (.not. file_exists) then
      ! Could not find ndpp_lib.xml file
      call fatal_error("NDPP Library XML file '" // trim(ndpp_lib) // &
           "' does not exist!")
    end if

    call write_message("Reading NDPP Library XML file...", 5)

    ! Parse ndpp_lib.xml file
    call open_xmldoc(doc, ndpp_lib)

    if (check_for_node(doc, "directory")) then
      ! Copy directory information if present
      call get_node_value(doc, "directory", directory)
    else
      ! If no directory is listed in ndpp_lib.xml, by default select the
      ! directory in which the ndpp_lib.xml file resides
      i = index(ndpp_lib, "/", BACK=.true.)
      directory = ndpp_lib(1:i)
    end if

    ! determine whether binary, ascii, or hdf5
    temp_str = ''
    if (check_for_node(doc, "filetype")) &
         call get_node_value(doc, "filetype", temp_str)
    if (trim(temp_str) == 'ascii') then
      filetype = ASCII
    elseif (trim(temp_str) == 'binary') then
      filetype = BINARY
    elseif (trim(temp_str) == 'hdf5') then
      filetype = H5
    else
      call fatal_error("Unknown filetype in " // trim(ndpp_lib) // &
                       ": " // trim(temp_str))
    end if

    ! Test metadata to ensure this library matches the problem definition
    ! First test to see if we have a legendre scattering type and not tabular
    ! (This will be removed at an undefined date when OpenMC supports
    ! tabular scattering distributions as opposed to just Legendres.
    ! Note that the authors have no immediate plans for this).
    if (check_for_node(doc, "scatt_type")) &
         call get_node_value(doc, "scatt_type", scatt_type)
    if (scatt_type /= SCATT_TYPE_LEGENDRE) then
      call fatal_error("Invalid Scattering Type represented in NDPP &
           &data. Rerun NDPP with the Legendre scattering type set.")
    end if

    ! Test to ensure the scattering order is less than the maximum
    order = 0
    if (check_for_node(doc, "scatt_order")) &
         call get_node_value(doc, "scatt_order", order)
    if (order > MAX_ANG_ORDER) then
      call fatal_error("Invalid scattering order of " // &
           trim(to_str(order)) // " requested.")
    end if

    ! Check the energy bin structure
    if (check_for_node(doc, "energy_bins")) then
      ndpp_groups = get_arraysize_double(doc, "energy_bins")
      allocate(ndpp_energy_bins(ndpp_groups))
      ndpp_groups = ndpp_groups - 1
      call get_node_array(doc, "energy_bins", ndpp_energy_bins)
    end if

    ! Check if nuscatter is present, obtain if it is
    temp_str = ""
    if (check_for_node(doc, "nuscatter")) then
      call get_node_value(doc, "nuscatter", temp_str)
      if (trim(temp_str) == 'true') then
        nuscatter = .true.
      else
        nuscatter = .false.
      end if
    else
      nuscatter = .false.
    end if

    ! Check if chi is present, obtain if it is
    temp_str = ""
    if (check_for_node(doc, "chi_present")) then
      call get_node_value(doc, "chi_present", temp_str)
      if (trim(temp_str) == 'true') then
        chi_present = .true.
      else
        chi_present = .false.
      end if
    else
      chi_present = .false.
    end if

    ! Check the tallies to ensure that the energy group structure, scattering
    ! order requested are valid (i.e., groups match, orders are less than in
    ! the library)
    TALLY_LOOP: do i = 1, n_tallies
      t => tallies(i)
      j = 0
      SCORE_LOOP: do k = 1, t % n_user_score_bins
        j = j + 1

        ! Check to see if the correct filters and orders are requested
        select case (t % score_bins(j))
        case (SCORE_NDPP_SCATT_N, SCORE_NDPP_NU_SCATT_N)
          ! We found the correct score, get comparing!
          ! First check the scattering order
          if (order < t % moment_order(j)) then
            call fatal_error("Invalid scattering order of " // &
                 trim(to_str(t % moment_order(j))) // &
                 " requested. Order requested is larger than &
                 &provided in the library (" // trim(to_str(order)) // ")!")
          end if

          ! Compare the energyin and energyout filters of this tally to the
          ! energy_bins_ metadata of the NDPP library.
          ! Check the energyin filter first.
          i_filter = t % find_filter(FILTER_ENERGYIN)
          ! We have already checked to ensure some energyin filter exists,
          ! reducing the cases to check. First we check the size, so that we
          ! can use the Fortran intrinsic ALL to check the actual values
          if (t % filters(i_filter) % n_bins /= ndpp_groups) then
            call fatal_error("Number of groups in NDPP Library do " // &
                 "not match that requested in tally!")
          end if
          ! Now we can check the actual group boundaries.
          if (any(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
            call fatal_error("NDPP Library group structure does not " // &
                 "match that requested in tally!")
          end if
          ! Repeat the same steps as above,
          ! but this time for the energyout filter
          i_filter = t % find_filter(FILTER_ENERGYOUT)
          if (t % filters(i_filter) % n_bins /= ndpp_groups) then
            call fatal_error("Number of groups in NDPP Library do not " // &
                 "match that requested in tally!")
          end if
          if (any(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
            call fatal_error("NDPP Library group structure does not " // &
                 "match that requested in tally!")
          end if
        case (SCORE_NDPP_SCATT_PN, SCORE_NDPP_NU_SCATT_PN)
          ! We found the correct score, get comparing!
          ! First check the scattering order
          if (order < t % moment_order(j)) then
            call fatal_error("Invalid scattering order of " // &
                 trim(to_str(t % moment_order(j))) // &
                 " requested. Order requested is larger than " // &
                 "provided in the library (" // &
                 trim(to_str(order)) // ")!")
          end if

          ! Compare the energyin and energyout filters of this tally to the
          ! energy_bins_ metadata of the NDPP library.
          ! Check the energyin filter first.
          i_filter = t % find_filter(FILTER_ENERGYIN)
          ! We have already checked to ensure some energyin filter exists,
          ! reducing the cases to check. First we check the size, so that we
          ! can use the Fortran intrinsic ALL to check the actual values
          if (t % filters(i_filter) % n_bins /= ndpp_groups) then
            call fatal_error("Number of groups in NDPP Library do " // &
                 "not match that requested in tally!")
          end if
          ! Now we can check the actual group boundaries.
          if (any(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
            call fatal_error("NDPP Library group structure does not " // &
                 "match that requested in tally!")
          end if
          ! Repeat the same steps as above,
          ! but this time for the energyout filter
          i_filter = t % find_filter(FILTER_ENERGYOUT)
          if (t % filters(i_filter) % n_bins /= ndpp_groups) then
            call fatal_error("Number of groups in NDPP Library do " // &
                 "not match that requested in tally!")
          end if
          if (any(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
            call fatal_error("NDPP Library group structure does not " // &
                 "match that requested in tally!")
          end if

        case (SCORE_NDPP_SCATT_YN, SCORE_NDPP_NU_SCATT_YN)
          ! We found the correct score, get comparing!
          ! First check the scattering order
          if (order < t % moment_order(j)) then
            call fatal_error("Invalid scattering order of " // &
                 trim(to_str(t % moment_order(j))) // &
                 " requested. Order requested is larger than " // &
                 "provided in the library (" // &
                 trim(to_str(order)) // ")!")
          end if

          ! Compare the energyin and energyout filters of this tally to the
          ! energy_bins_ metadata of the NDPP library.
          ! Check the energyin filter first.
          i_filter = t % find_filter(FILTER_ENERGYIN)
          ! We have already checked to ensure some energyin filter exists,
          ! reducing the cases to check. First we check the size, so that we
          ! can use the Fortran intrinsic ALL to check the actual values
          if (t % filters(i_filter) % n_bins /= ndpp_groups) then
            call fatal_error("Number of groups in NDPP Library do " // &
                 "not match that requested in tally " // to_str(t % id) // ".")
          end if
          ! Now we can check the actual group boundaries.
          if (any(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
            call fatal_error("NDPP Library group structure does " // &
                 "not match that requested in tally " // to_str(t % id) // ".")
          end if
          ! Repeat the same steps as above,
          ! but this time for the energyout filter
          i_filter = t % find_filter(FILTER_ENERGYOUT)
          if (t % filters(i_filter) % n_bins /= ndpp_groups) then
            call fatal_error("Number of groups in NDPP Library do " // &
                 "not match that requested in tally " // to_str(t % id) // ".")
          end if
          if (any(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
            call fatal_error("NDPP Library group structure does " // &
                 "not match that requested in tally " // to_str(t % id) // ".")
          end if

        case (SCORE_NDPP_CHI, SCORE_NDPP_CHI_P, SCORE_NDPP_CHI_D)
          ! Check that the group structure matches
          i_filter = t % find_filter(FILTER_ENERGYOUT)
          if (t % filters(i_filter) % n_bins /= ndpp_groups) then
            call fatal_error("Number of groups in NDPP Library do " // &
                 "not match that requested in tally " // to_str(t % id) // ".")
          end if
          if (any(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
            call fatal_error("NDPP Library group structure does not " // &
                 "match that requested in tally " // to_str(t % id) // ".")
          end if
        end select

        ! Adjust indices within the score_bins array (j)
        select case (t % score_bins(j))
        case (SCORE_NDPP_SCATT_PN, SCORE_NDPP_NU_SCATT_PN)
          j = j + t % moment_order(j)
          cycle SCORE_LOOP ! Skip the others to save cycles

        case (SCORE_NDPP_SCATT_YN, SCORE_NDPP_NU_SCATT_YN)
          j = j + (t % moment_order(j) + 1)**2
          cycle SCORE_LOOP ! Skip the others to save cycles

        end select
      end do SCORE_LOOP
    end do TALLY_LOOP

    ! Get node list of all <ndpp_table> entries
    call get_node_list(doc, "ndpp_table", node_ndpp_list)
    n_listings = get_list_size(node_ndpp_list)

    ! Allocate ndpp_listings array
    if (n_listings == 0) then
      call fatal_error("No NDPP table listings present in ndpp_lib.xml file!")
    else
      allocate(ndpp_listings(n_listings))
    end if

    do i = 1, n_listings
      listing => ndpp_listings(i)
      ! Get pointer to ace table XML node
      call get_list_item(node_ndpp_list, i, node_ndpp)

      ! copy a number of attributes
      call get_node_value(node_ndpp, "name", listing % name)
      if (check_for_node(node_ndpp, "alias")) &
           call get_node_value(node_ndpp, "alias", listing % alias)
      call get_node_value(node_ndpp, "zaid", listing % zaid)
      call get_node_value(node_ndpp, "awr", listing % awr)
      if (check_for_node(node_ndpp, "temperature")) &
           call get_node_value(node_ndpp, "temperature", listing % kT)
      call get_node_value(node_ndpp, "location", listing % location)

      ! determine type of cross section
      if (ends_with(listing % name, 'c')) then
        listing % type = ACE_NEUTRON
      elseif (ends_with(listing % name, 't')) then
        listing % type = ACE_THERMAL
      end if

      ! set filetype
      if (check_for_node(node_ndpp, "filetype")) then
        temp_str = ''
        call get_node_value(node_ndpp, "filetype", temp_str)
        if (temp_str == 'ascii') then
          listing % filetype = ASCII
        else if (temp_str == 'binary') then
          listing % filetype = BINARY
        end if
      else
        listing % filetype = filetype
      end if

      ! determine metastable state
      if (.not.check_for_node(node_ndpp, "metastable")) then
        listing % metastable = .false.
      else
        listing % metastable = .true.
      end if

      ! determine path of ndpp data table
      if (check_for_node(node_ndpp, "path")) then
        call get_node_value(node_ndpp, "path", temp_str)
      else
        call fatal_error("Path missing for isotope " // listing % name)
      end if

      if (starts_with(temp_str, '/')) then
        listing % path = trim(temp_str)
      else
        if (ends_with(directory,'/')) then
          listing % path = trim(directory) // trim(temp_str)
        else
          listing % path = trim(directory) // '/' // trim(temp_str)
        end if
      end if

      ! create dictionary entry for both name and alias
      call ndpp_listing_dict % add_key(listing % name, i)
      if (check_for_node(node_ndpp, "alias")) then
        call ndpp_listing_dict % add_key(listing % alias, i)
      end if
    end do

    ! Close ndpp_lib XML file
    call close_xmldoc(doc)

  end subroutine read_ndpp_xml

end module ndpp_initialize
