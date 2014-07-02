module ndpp

  use ace_header,       only: Nuclide, GrpTransfer, XsListing
  use constants
  use dict_header,      only: DictCharInt, DICT_NULL
  use error,            only: fatal_error, warning
  use global
  use output,           only: write_message
  use search
  use string,           only: ends_with, lower_case, starts_with, to_str
  use xml_interface

  implicit none

  ! Module global data:

  ! ndpp_lib.xml preprocessed data listings and associated data.
  type(XsListing),  allocatable, target :: ndpp_listings(:)
  type(DictCharInt)                     :: ndpp_listing_dict

contains

!===============================================================================
! READ_NDPP_LIB reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_ndpp_data()
    type(Nuclide), pointer    :: nuc => null() ! Current working nuclide
    type(SAlphaBeta), pointer :: sab => null() ! Current working SAB table
    type(XsListing), pointer  :: ndpp_listing => null()
    integer :: i_listing     ! index in ndpp_listings array
    integer :: i_nuclide     ! index in nuclides
    integer :: i_sab         ! index in sab_tables
    integer :: scatt_type    ! Whether or not legendre or tabular data
    logical :: get_scatt     ! Flag for whether or not to get scatt data
    logical :: get_nuscatt   ! Flag for whether or not to get nuscatt data
    logical :: get_chi_t     ! Flag for whether or not to get total chi data
    logical :: get_chi_p     ! Flag for whether or not to get prompt chi data
    logical :: get_chi_d     ! Flag for whether or not to get delayed chi data
    integer :: scatt_order   ! Number of moments requested in tallies for scatter
    integer :: ndpp_groups   ! Number of groups in the NDPP library

    ! First lets go read the ndpp_lib.xml file
    call read_ndpp_xml(scatt_type, ndpp_groups)

    ! Determine which data is required to be stored as well as the maximum orders
    ! for scatt and nuscatt
    call which_data(scatt_type, get_scatt, get_nuscatt, get_chi_t, get_chi_p, &
                    get_chi_d, scatt_order)

    ! Parse through each nuclide in the model and read in the corresponding
    ! NDPP data.
    ! During this we will make sure the temperatures of ACE and NDPP match,
    ! as well as making sure the NDPP data exists at all.

    do i_nuclide = 1, n_nuclides_total
      nuc => nuclides(i_nuclide)
      i_listing = ndpp_listing_dict % get_key(adjustl(trim(nuc % name)))
      if (i_listing == DICT_NULL) then
        ! Could not find ndpp_lib.xml file
        message = trim(nuc % name) // " does not exist in NDPP XML file: '" // &
                  trim(ndpp_lib) // "'!"
        call fatal_error()
      end if
      ! read_ndpp_table will read the NDPP data and also check that
      ! the temperatures match
      ndpp_listing => ndpp_listings(i_listing)
      call read_ndpp_table(ndpp_listing, get_scatt, get_nuscatt, get_chi_t, &
                           get_chi_p, get_chi_d, scatt_order, NUC=nuc)
    end do

    do i_sab = 1, n_sab_tables
      sab => sab_tables(i_sab)
      i_listing = ndpp_listing_dict % get_key(adjustl(trim(sab % name)))
      if (i_listing == DICT_NULL) then
        ! Could not find ndpp_lib.xml file
        message = trim(sab % name) // " does not exist in NDPP XML file: '" // &
                  trim(ndpp_lib) // "'!"
        call fatal_error()
      end if
      ! read_ndpp_table will read the NDPP data and also check that
      ! the temperatures match
      ndpp_listing => ndpp_listings(i_listing)
      call read_ndpp_table(ndpp_listing, get_scatt, get_nuscatt, get_chi_t, &
                           get_chi_p, get_chi_d, scatt_order, SAB=sab)
    end do

    if (allocated(ndpp_listings)) then
      deallocate(ndpp_listings)
    end if
    call ndpp_listing_dict % clear()

    ! Now allocate ndpp_outgoing, our `scratch` variable to store the combined
    ! and interpolated elastic + inelastic data.  This is program global and has
    ! been declared threadprivate.
    allocate(ndpp_outgoing(0: n_threads, scatt_order, ndpp_groups))

  end subroutine read_ndpp_data

!===============================================================================
! WHICH_DATA looks at the requested tallies and determines which data is needed
! from the NDPP libraries so that memory utilization is kept to a minimum.
!===============================================================================

  subroutine which_data(scatt_type, get_scatt, get_nuscatt, get_chi_t, &
                        get_chi_p, get_chi_d, scatt_order)
    integer, intent(in)  :: scatt_type  ! Whether or not legendre or tabular data
    logical, intent(out) :: get_scatt   ! Flag for whether or not to get scatt data
    logical, intent(out) :: get_nuscatt ! Flag for whether or not to get nuscatt data
    logical, intent(out) :: get_chi_t   ! Flag for whether or not to get total chi data
    logical, intent(out) :: get_chi_p   ! Flag for whether or not to get prompt chi data
    logical, intent(out) :: get_chi_d   ! Flag for whether or not to get delayed chi data
    integer, intent(out) :: scatt_order ! Number of moments requested in tallies for scatter

    type(TallyObject), pointer :: t => null()
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
          case (SCORE_NDPP_SCATT_N, SCORE_NDPP_SCATT_PN)
            get_scatt = .true.
            if (t % moment_order(j) > scatt_order) then
              scatt_order = t % moment_order(j)
            end if

            if (t % score_bins(j) == SCORE_NDPP_SCATT_PN) then
              j = j + t % moment_order(j)
              cycle SCORE_LOOP ! Skip the others to save cycles
            end if

          case (SCORE_NDPP_NU_SCATT_N, SCORE_NDPP_NU_SCATT_PN)
            get_nuscatt = .true.
            if (t % moment_order(j) > scatt_order) then
              scatt_order = t % moment_order(j)
            end if

            if (t % score_bins(j) == SCORE_NDPP_NU_SCATT_PN) then
              j = j + t % moment_order(j)
              cycle SCORE_LOOP ! Skip the others to save cycles
            end if

          case (SCORE_NDPP_CHI)
            get_chi_t = .true.
          case (SCORE_NDPP_CHI_P)
            get_chi_p = .true.
          case (SCORE_NDPP_CHI_D)
            get_chi_d = .true.

          ! Cycle through j if necessary
          case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN, SCORE_FLUX_YN, &
                SCORE_TOTAL_YN, SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN)
            j = j + t % moment_order(j)
            cycle SCORE_LOOP ! Skip the others to save cycles
        end select
      end do SCORE_LOOP
    end do TALLY_LOOP

    ! Adjust the order so array sizing is correct
    if (scatt_type == SCATT_TYPE_LEGENDRE) then
      scatt_order = scatt_order + 1
    end if

  end subroutine which_data

!===============================================================================
! READ_NDPP_XML reads information from a cross_sections.xml file. This
! file contains a listing of the ACE cross sections that may be used.
!===============================================================================

  subroutine read_ndpp_xml(scatt_type, ndpp_groups)

    integer, intent(out)  :: scatt_type  ! Whether or not legendre or tabular data
    integer, intent(out)  :: ndpp_groups ! number of groups in NDPP data

    integer :: i, j, k     ! loop indices
    logical :: file_exists ! does ndpp_lib.xml exist?
    integer :: filetype    ! default file type
    logical :: nuscatter   ! is nuscatter data present?
    logical :: chi_present ! is chi data present?
    character(MAX_WORD_LEN)  :: directory   ! directory with cross sections
    character(MAX_LINE_LEN)  :: temp_str
    type(TallyObject), pointer :: t => null()
    integer :: i_filter    ! Index of filter which contains energyin or energyout
    ! We can use the same XSListing type for our ndpp data since the NDPP
    ! data is a subset of whats in cross_sections.xml
    type(XsListing), pointer :: listing => null()
    type(Node), pointer      :: doc => null()
    type(Node), pointer      :: node_ndpp => null()
    type(NodeList), pointer  :: node_ndpp_list => null()
    real(8), allocatable     :: ndpp_energy_bins(:)
    integer                  :: order       ! ndpp_lib.xml's scattering order

    ! Check if ndpp_lib.xml exists
    inquire(FILE=ndpp_lib, EXIST=file_exists)
    if (.not. file_exists) then
       ! Could not find ndpp_lib.xml file
       message = "NDPP Library XML file '" // trim(ndpp_lib) // &
                 "' does not exist!"
       call fatal_error()
    end if

    message = "Reading NDPP Library XML file..."
    call write_message(5)

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
       message = "Unknown filetype in " // trim(ndpp_lib) // &
        ": " // trim(temp_str)
       call fatal_error()
    end if

    ! Test metadata to ensure this library matches the problem definition
    ! First test to see if we have a legendre scattering type and not tabular
    !!! (This will be removed at an undefined data when OpenMC supports
    !!! tabular scattering distributions as opposed to Legendres.
    if (check_for_node(doc, "scatt_type")) &
      call get_node_value(doc, "scatt_type", scatt_type)
    if (scatt_type /= SCATT_TYPE_LEGENDRE) then
      message = "Invalid Scattering Type represented in NDPP data. Rerun " // &
                "NDPP with the Legendre scattering type set."
      call fatal_error()
    end if

    ! Test to ensure the scattering order is less than the maximum
    order = 0
    if (check_for_node(doc, "scatt_order")) &
      call get_node_value(doc, "scatt_order", order)
    if (order > MAX_ANG_ORDER) then
      message = "Invalid scattering order of " // trim(to_str(order)) // &
                " requested."
      call fatal_error()
    end if

    ! Get the energy bin structure
    if (check_for_node(doc, "energy_bins")) then
      ndpp_groups = get_arraysize_double(doc, "energy_bins")
      allocate(ndpp_energy_bins(ndpp_groups))
      ndpp_groups = ndpp_groups - 1
      call get_node_array(doc, "energy_bins", ndpp_energy_bins)
    end if

    ! Get if nuscatter is present
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

    ! Get if chi is present
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
        select case (t % score_bins(j))
          case (SCORE_SCATTER_PN)
            j = j + t % moment_order(j)
            cycle SCORE_LOOP ! Skip the others to save cycles
          case (SCORE_NDPP_SCATT_N, SCORE_NDPP_NU_SCATT_N)
            ! We found the correct score, get comparing!
            ! First check the scattering order
            if (order < t % moment_order(j)) then
              message = "Invalid scattering order of " // &
                        trim(to_str(t % moment_order(j))) // " requested. Order " // &
                        "requested is larger than provided in the library (" // &
                        trim(to_str(order)) // ")!"
              call fatal_error()
            end if

            ! Compare the energyin and energyout filters of this tally to the
            ! energy_bins_ metadata of the NDPP library.
            ! Check the energyin filter first.
            i_filter = t % find_filter(FILTER_ENERGYIN)
            ! We have already checked to ensure some energyin filter exists,
            ! reducing the cases to check. First we check the size, so that we
            ! can use the Fortran intrinsic ALL to check the actual values
            if (t % filters(i_filter) % n_bins /= ndpp_groups) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            ! Now we can check the actual group boundaries.
            if (all(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            ! Repeat the same steps as above, but this time for the energyout filter
            i_filter = t % find_filter(FILTER_ENERGYOUT)
            if (t % filters(i_filter) % n_bins /= ndpp_groups) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            if (all(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
          case (SCORE_NDPP_SCATT_PN, SCORE_NDPP_NU_SCATT_PN)
            ! We found the correct score, get comparing!
            ! First check the scattering order
            if (order < t % moment_order(j)) then
              message = "Invalid scattering order of " // &
                        trim(to_str(t % moment_order(j))) // " requested. Order " // &
                        "requested is larger than provided in the library (" // &
                        trim(to_str(order)) // ")!"
              call fatal_error()
            end if

            ! Compare the energyin and energyout filters of this tally to the
            ! energy_bins_ metadata of the NDPP library.
            ! Check the energyin filter first.
            i_filter = t % find_filter(FILTER_ENERGYIN)
            ! We have already checked to ensure some energyin filter exists,
            ! reducing the cases to check. First we check the size, so that we
            ! can use the Fortran intrinsic ALL to check the actual values
            if (t % filters(i_filter) % n_bins /= ndpp_groups) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            ! Now we can check the actual group boundaries.
            if (all(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            ! Repeat the same steps as above, but this time for the energyout filter
            i_filter = t % find_filter(FILTER_ENERGYOUT)
            if (t % filters(i_filter) % n_bins /= ndpp_groups) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            if (all(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if

            j = j + t % moment_order(j)
            cycle SCORE_LOOP ! Skip the others to save cycles

          case (SCORE_NDPP_CHI, SCORE_NDPP_CHI_P, SCORE_NDPP_CHI_D)
            ! Check that the group structure matches
            i_filter = t % find_filter(FILTER_ENERGYOUT)
            if (t % filters(i_filter) % n_bins /= ndpp_groups) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            if (all(t % filters(i_filter) % real_bins /= ndpp_energy_bins)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
        end select
      end do SCORE_LOOP
    end do TALLY_LOOP

    ! Get node list of all <ndpp_table> entries
    call get_node_list(doc, "ndpp_table", node_ndpp_list)
    n_listings = get_list_size(node_ndpp_list)

    ! Allocate ndpp_listings array
    if (n_listings == 0) then
       message = "No NDPP table listings present in ndpp_lib.xml file!"
       call fatal_error()
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

       ! set filetype, record length, and number of entries
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
         message = "Path missing for isotope " // listing % name
         call fatal_error()
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

!===============================================================================
! READ_NDPP_TABLE reads information from a pre-processed nuclear data file
! as generated by the NDPP program.  The result of this subroutine is that the
! Nuclide % ndpp_el (or SAlphaBeta % ndpp_el) array will be allocated and
! filled with the corresponding pre-processed scattering (and TODO: Chi) data
! from NDPP. Versions of this routine exist for.  These are currently written
! with lots of branching statements throughout to maximize the number of
! lines that can be shared.  This si done with the expectation that in the
! future, nuc and sab will be extended types of a base type.
!===============================================================================

  subroutine read_ndpp_table(listing, get_scatt, get_nuscatt, get_chi_t, &
                             get_chi_p, get_chi_d, scatt_order, nuc, sab)
    type(XsListing),  pointer, intent(in) :: listing    ! Current NDPP data
    logical, intent(in) :: get_scatt   ! Whether or not to get scatt data
    logical, intent(in) :: get_nuscatt ! Whether or not to get nuscatt data
    logical, intent(in) :: get_chi_t   ! Whether or not to get total chi data
    logical, intent(in) :: get_chi_p   ! Whether or not to get prompt chi data
    logical, intent(in) :: get_chi_d   ! Whether or not to get delayed chi data
    integer, intent(in) :: scatt_order ! Number of moments requested in tallies
    type(Nuclide),    pointer, optional, intent(inout) :: nuc ! Current nuclide
    type(SAlphaBeta), pointer, optional, intent(inout) :: sab ! Current Sab data

    integer       :: in = 7        ! file unit
    integer       :: i             ! loop index for data records
    integer       :: location      ! location of NDPP table
    integer       :: filetype      ! Ascii, Binary or HDF5 filetype
    logical       :: file_exists   ! does NDPP library exist?
    character(7)  :: readable      ! is NDPP library readable?
    character(10) :: name          ! name of NDPP table
    real(8)       :: kT            ! temperature of table
    real(8)       :: dkT           ! difference in temperature of table
    integer       :: NG            ! Number of energy groups in library
    character(MAX_FILE_LEN) :: filename ! path to NDPP data library
    real(8), allocatable    :: energy_bins(:) ! Energy group structure
    integer       :: mu_bins       ! NUmber of angular points used
    integer       :: scatt_type    ! Type of scattering data, discarded
    integer       :: lib_order     ! Order of scattering data in library
    integer       :: nuscatter     ! Flag as to if nuscatter data is present
    integer       :: chi_present   ! Flag as to if chi data is present
    integer       :: gmin, gmax    ! Min and max possible group transfers
    integer       :: g             ! Group index
    real(8)       :: thin_tol      ! Thinning tolerance used in lib, discarded
    integer       :: NEin, iE      ! Number of incoming energies and the index
    real(8), allocatable :: temp_outgoing(:,:) ! Temporary storage of scatt data
    logical       :: is_nuc        ! Is our data a nuc or an sab?
    integer       :: NP            ! Number of precursors groups

    ! Set is_nuc based on optional params
    if ((present(nuc)) .and. (present(sab))) then
      message = "Invalid usage of READ_NDPP_TABLE! Cannot accept both nuc &
                &and sab arguments!"
      call fatal_error()
    else if ((.not. present(nuc)) .and. (.not. present(sab))) then
      message = "Invalid usage of READ_NDPP_TABLE! Subroutine must be provided &
                & a nuc or an sab argument!"
      call fatal_error()
    else if (present(nuc)) then
      is_nuc = .true.
    else if (present(sab)) then
      is_nuc = .false.
    end if


    ! determine path, record length, and location of table
    filename      = listing % path
    filetype      = listing % filetype
    location      = listing % location

    ! Check if ACE library exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
      message = "NDPP library '" // trim(filename) // "' does not exist!"
      call fatal_error()
    elseif (readable(1:3) == 'NO') then
      message = "NDPP library '" // trim(filename) // "' is not readable! &
           &Change file permissions with chmod command."
      call fatal_error()
    end if

    ! display message
    message = "Loading NDPP data library: " // listing % name
    call write_message(6)

    if (listing % filetype == ASCII) then
      ! =======================================================================
      ! READ NDPP DATA IN ASCII FORMAT

      ! Find location of table
      open(UNIT=in, FILE=filename, STATUS='old', ACTION='read')
      rewind(UNIT=in)
      do i = 1, location - 1
        read(UNIT=in, FMT=*)
      end do

      ! Read the header information to make sure this is the correct file
      read(UNIT=in, FMT='(A20,1PE20.12,I20,A20)') name, kT, NG
      dkT = abs(kT - listing % kT) / kT
      if ((adjustl(trim(name)) /= listing % name) .or. (dkT < 0.01_8)) then
        message = "NDPP library '" // trim(filename) // "' does not &
                  &contain the correct data set where expected!"
        call fatal_error()
      end if
      ! Get the energy bins (not checking, will assume from here on out we have
      ! the right data set)
      allocate(energy_bins(NG + 1))
      read(UNIT=in, FMT=*) energy_bins
      ! The next line is scatt_type, lib_order, nuscatter, chi_present
      read(UNIT=in, FMT='(I20,I20,I20,I20)') scatt_type, lib_order, &
        nuscatter, chi_present
      ! Finally, mu_bins, thin_tol
      read(UNIT=in, FMT='(I20,1PE20.12)') mu_bins, thin_tol

      ! set lib_order to the right number for allocating the outgoing array
      if (scatt_type == SCATT_TYPE_LEGENDRE) lib_order = lib_order + 1

      ! Start with elastic data
      ! Get Ein information

      ! First find our data length, doesnt depend on nuc or sab
      read(UNIT=in, FMT='(I20)') NEin
      if (is_nuc) then
        ! Now allocate all that will be filled
        allocate(nuc % ndpp_el_Ein(NEin))
        allocate(nuc % ndpp_el_Ein_srch(NG + 1))
        allocate(nuc % ndpp_el(NEin))

        ! Now read in el_Ein and el_Ein_srch
        read(UNIT=in, FMT=*) nuc % ndpp_el_Ein_srch
        read(UNIT=in, FMT=*) nuc % ndpp_el_Ein

      else
        ! Now allocate all that will be filled
        allocate(sab % ndpp_el_Ein(NEin))
        allocate(sab % ndpp_el_Ein_srch(NG + 1))
        allocate(sab % ndpp_el(NEin))

        ! Now read in el_Ein and el_Ein_srch
        read(UNIT=in, FMT=*) sab % ndpp_el_Ein
        read(UNIT=in, FMT=*) sab % ndpp_el_Ein_srch
      end if

      ! Get the elastic moments themselves
      do iE = 1, NEin
        ! get gmin and gmax
        read(UNIT=in, FMT='(I20,I20)') gmin, gmax

        if ((gmin > 0) .and. (gmax > 0)) then
          ! Then we can allocate the space. Do it to ndpp_el_order
          ! since this is the largest order requested in the tallies.
          ! Since we only need to store up to the maximum, we also need to have
          ! an array for reading the file which we can later truncate to fit
          ! in to nuc/sab % ndpp_el(iE) % outgoing.
          allocate(temp_outgoing(lib_order, gmin : gmax))

          ! Now we have a space to store the data, get it.
          read(UNIT=in, FMT=*) temp_outgoing
          ! And copy in to nuc/sab % ndpp_el
          if (is_nuc) then
            allocate(nuc % ndpp_el(iE) % outgoing(scatt_order, gmin : gmax))
            nuc % ndpp_el(iE) % outgoing(:, gmin : gmax) = &
              temp_outgoing(1 : scatt_order, gmin : gmax)
          else
            allocate(sab % ndpp_el(iE) % outgoing(scatt_order, gmin : gmax))
            sab % ndpp_el(iE) % outgoing(:, gmin : gmax) = &
              temp_outgoing(1 : scatt_order, gmin : gmax)
          end if
          deallocate(temp_outgoing)
        end if
      end do

      ! The remainder only apply to nuclides (inelastic, nu-inelastic and chi data)
      if (is_nuc) then
        read(UNIT=in, FMT='(I20)') NEin
        if (NEin > 0) then
          ! Now allocate all that will be filled
          allocate(nuc % ndpp_inel_Ein(NEin))
          allocate(nuc % ndpp_inel_Ein_srch(NG + 1))
          allocate(nuc % ndpp_sigma_inel(NEin))
          allocate(nuc % ndpp_inel(NEin))

          ! Now read in inel_Ein, inel_Ein_srch, and sigma_inel
          read(UNIT=in, FMT=*) nuc % ndpp_inel_Ein
          read(UNIT=in, FMT=*) nuc % ndpp_inel_Ein_srch
          read(UNIT=in, FMT=*) nuc % ndpp_sigma_inel

          ! Get the inelastic moments themselves
          do iE = 1, NEin
            ! get gmin and gmax
            read(UNIT=in, FMT='(I20,I20)') gmin, gmax

            if ((gmin > 0) .and. (gmax > 0)) then
              ! Then we can allocate the space. Do it to ndpp_el_order
              ! since this is the largest order requested in the tallies.
              ! Since we only need to store up to the maximum, we also need to have
              ! an array for reading the file which we can later truncate to fit
              ! in to nuc/sab % ndpp_el(iE) % outgoing.
              allocate(temp_outgoing(lib_order, gmin : gmax))

              ! Now we have a space to store the data, get it.
              read(UNIT=in, FMT=*) temp_outgoing
              ! And copy in to nuc % ndpp_el
              allocate(nuc % ndpp_el(iE) % outgoing(scatt_order, gmin : gmax))
              nuc % ndpp_el(iE) % outgoing(:, gmin : gmax) = &
                temp_outgoing(1 : scatt_order, gmin : gmax)
              deallocate(temp_outgoing)
            end if
          end do

          ! Get nu-scatter, if needed
          if (nuscatter == 1) then
            allocate(nuc % ndpp_nuinel(NEin))
            do iE = 1, NEin
              ! get gmin and gmax
              read(UNIT=in, FMT='(I20,I20)') gmin, gmax

              if ((gmin > 0) .and. (gmax > 0)) then
                ! Then we can allocate the space. Do it to ndpp_el_order
                ! since this is the largest order requested in the tallies.
                ! Since we only need to store up to the maximum, we also need to have
                ! an array for reading the file which we can later truncate to fit
                ! in to nuc/sab % ndpp_el(iE) % outgoing.
                allocate(temp_outgoing(lib_order, gmin : gmax))

                ! Now we have a space to store the data, get it.
                read(UNIT=in, FMT=*) temp_outgoing
                allocate(nuc % ndpp_nuinel(iE) % outgoing(scatt_order, &
                  gmin : gmax))
                nuc % ndpp_nuinel(iE) % outgoing(:, gmin : gmax) = &
                  temp_outgoing(1 : scatt_order, gmin : gmax)
                deallocate(temp_outgoing)
              end if
            end do
          end if
        end if

        ! Get chi(E_{in}) data if provided
        if (chi_present == 1) then
          ! Get Ein grid and number of precursors
          read(UNIT=in, FMT='(I20,I20)') NEin, NP

          ! Allocate data as needed
          allocate(nuc % ndpp_chi_Ein(NEin))
          allocate(nuc % ndpp_chi(size(energy_bins) - 1, NEin))
          allocate(nuc % ndpp_chi_p(size(energy_bins) - 1, NEin))
          allocate(nuc % ndpp_chi_d(NP, size(energy_bins) - 1, NEin))

          ! Get Ein Grid
          read(UNIT=in, FMT=*) nuc % ndpp_chi_Ein

          ! Get Chi-Total
          read(UNIT=in, FMT=*) nuc % ndpp_chi

          ! Get Chi-Prompt
          read(UNIT=in, FMT=*) nuc % ndpp_chi_p

          ! Get Chi-Delayed
          do i = 1, NP
            do iE = 1, NEin
              do g = 1, size(energy_bins) - 1
                read(UNIT=in, FMT=*) nuc % ndpp_chi_d(i, g, iE)
              end do
            end do
          end do
        end if
      end if

      close(UNIT=in)

    else if (listing % filetype == BINARY) then
      ! =======================================================================
      ! READ NDPP DATA IN BINARY FORMAT

      ! Open file
      open(UNIT=in, FILE=filename, STATUS='old', ACTION='read', ACCESS='stream')
      rewind(UNIT=in)
      !!! right now binary files dont have the capability to do location /= 1!

      ! Read the header information to make sure this is the correct file
      read(UNIT=in) name, kT, NG
      dkT = abs(kT - listing % kT) / kT
      if ((adjustl(trim(name)) /= listing % name) .or. (dkT > 0.001_8)) then
        message = "NDPP library '" // trim(filename) // "' does not &
                  &contain the correct data set where expected!"
        call fatal_error()
      end if
      ! Get the energy bins (not checking, will assume from here on out we have
      ! the right data set), and the other meta information
      allocate(energy_bins(NG + 1))
      read(UNIT=in) energy_bins, scatt_type, lib_order, &
        nuscatter, chi_present, mu_bins, thin_tol

      ! set lib_order to the right number for allocating the outgoing array
      if (scatt_type == SCATT_TYPE_LEGENDRE) lib_order = lib_order + 1

      ! Get \sigma_{s,g'->g,l}(E_{in}) data
      ! Get elastic Ein information

      ! First find our data length, doesnt depend on nuc or sab
      read(UNIT=in) NEin
      if (is_nuc) then
        ! Now allocate all that will be filled
        allocate(nuc % ndpp_el_Ein(NEin))
        allocate(nuc % ndpp_el_Ein_srch(NG + 1))
        allocate(nuc % ndpp_el(NEin))

        ! Now read in el_Ein and el_Ein_srch
        read(UNIT=in) nuc % ndpp_el_Ein
        read(UNIT=in) nuc % ndpp_el_Ein_srch

      else
        ! Now allocate all that will be filled
        allocate(sab % ndpp_el_Ein(NEin))
        allocate(sab % ndpp_el_Ein_srch(NG + 1))
        allocate(sab % ndpp_el(NEin))

        ! Now read in el_Ein and el_Ein_srch
        read(UNIT=in) sab % ndpp_el_Ein
        read(UNIT=in) sab % ndpp_el_Ein_srch
      end if

      ! Now we can get the elastic moments themselves one-by-one due to sparse
      ! storage nature.
      do iE = 1, NEin
        ! get gmin and gmax
        read(UNIT=in) gmin, gmax

        if ((gmin > 0) .and. (gmax > 0)) then
          ! Then we can allocate the space. Do it to ndpp_el_order
          ! since this is the largest order requested in the tallies.
          ! Since we only need to store up to the maximum, we also need to have
          ! an array for reading the file which we can later truncate to fit
          ! in to nuc/sab % ndpp_el(iE) % outgoing.
          allocate(temp_outgoing(lib_order, gmin : gmax))

          ! Now we have a space to store the data, get it.
          read(UNIT=in) temp_outgoing
          ! And copy in to nuc/sab % ndpp_el
          if (is_nuc) then
            allocate(nuc % ndpp_el(iE) % outgoing(scatt_order, gmin : gmax))
            nuc % ndpp_el(iE) % outgoing(:, gmin : gmax) = &
              temp_outgoing(1 : scatt_order, gmin : gmax)
          else
            allocate(sab % ndpp_el(iE) % outgoing(scatt_order, gmin : gmax))
            sab % ndpp_el(iE) % outgoing(:, gmin : gmax) = &
              temp_outgoing(1 : scatt_order, gmin : gmax)
          end if
          deallocate(temp_outgoing)
        end if
      end do

      ! The remainder only apply to nuclides (inelastic, nu-inelastic and chi data)
      if (is_nuc) then
        read(UNIT=in) NEin
        if (NEin > 0) then
          ! Now allocate all that will be filled
          allocate(nuc % ndpp_inel_Ein(NEin))
          allocate(nuc % ndpp_inel_Ein_srch(NG + 1))
          allocate(nuc % ndpp_sigma_inel(NEin))
          allocate(nuc % ndpp_inel(NEin))

          ! Now read in inel_Ein, inel_Ein_srch, and sigma_inel
          read(UNIT=in) nuc % ndpp_inel_Ein
          read(UNIT=in) nuc % ndpp_inel_Ein_srch
          read(UNIT=in) nuc % ndpp_sigma_inel

          ! Get the inelastic moments themselves
          do iE = 1, NEin
            ! get gmin and gmax
            read(UNIT=in) gmin, gmax

            if ((gmin > 0) .and. (gmax > 0)) then
              ! Then we can allocate the space. Do it to ndpp_inel_order
              ! since this is the largest order requested in the tallies.
              ! Since we only need to store up to the maximum, we also need to have
              ! an array for reading the file which we can later truncate to fit
              ! in to nuc/sab % ndpp_inel(iE) % outgoing.
              allocate(temp_outgoing(lib_order, gmin : gmax))

              ! Now we have a space to store the data, get it.
              read(UNIT=in) temp_outgoing
              ! And copy in to nuc % ndpp_inel
              allocate(nuc % ndpp_inel(iE) % outgoing(scatt_order, gmin : gmax))
              nuc % ndpp_inel(iE) % outgoing(:, gmin : gmax) = &
                temp_outgoing(1 : scatt_order, gmin : gmax)
              deallocate(temp_outgoing)
            end if
          end do

          ! Get nu-scatter, if needed
          if (nuscatter == 1) then
            allocate(nuc % ndpp_nuinel(NEin))
            do iE = 1, NEin
              ! get gmin and gmax
              read(UNIT=in) gmin, gmax

              if ((gmin > 0) .and. (gmax > 0)) then
                ! Then we can allocate the space. Do it to ndpp_el_order
                ! since this is the largest order requested in the tallies.
                ! Since we only need to store up to the maximum, we also need to have
                ! an array for reading the file which we can later truncate to fit
                ! in to nuc/sab % ndpp_nuinel(iE) % outgoing.
                allocate(temp_outgoing(lib_order, gmin : gmax))

                ! Now we have a space to store the data, get it.
                read(UNIT=in) temp_outgoing
                ! And copy in to nuc % ndpp_nuinel
                allocate(nuc % ndpp_nuinel(iE) % outgoing(scatt_order, &
                  gmin : gmax))
                nuc % ndpp_nuinel(iE) % outgoing(:, gmin : gmax) = &
                  temp_outgoing(1 : scatt_order, gmin : gmax)
                deallocate(temp_outgoing)
              end if
            end do
          end if
        end if

        ! Get chi data, if provided
        if (chi_present == 1) then
          ! Get Ein grid and Number of Precursors
          read(UNIT=in) NEin, NP

          ! Allocate data as needed
          allocate(nuc % ndpp_chi_Ein(NEin))
          allocate(nuc % ndpp_chi(size(energy_bins) - 1, NEin))
          allocate(nuc % ndpp_chi_p(size(energy_bins) - 1, NEin))
          allocate(nuc % ndpp_chi_d(NP, size(energy_bins) - 1, NEin))

          ! Get Ein Grid
          read(UNIT=in) nuc % ndpp_chi_Ein

          ! Get Chi-Total
          read(UNIT=in) nuc % ndpp_chi

          ! Get Chi-Prompt
          read(UNIT=in) nuc % ndpp_chi_p

          ! Get Chi-Delayed
          do i = 1, NP
            do iE = 1, NEin
              do g = 1, size(energy_bins) - 1
                read(UNIT=in) nuc % ndpp_chi_d(i, g, iE)
              end do
            end do
          end do
        end if
      end if

      close(UNIT=in)

    else if (listing % filetype == H5) then
      ! =======================================================================
      ! READ NDPP DATA IN HDF5 FORMAT
      !!! TBI
    end if ! No default needed - read_ndpp_xml() already checked filetype

    ! Finally, the above code read in all data since NDPP libs are sequential
    ! access; now we can deallocate what we do not need
    if (is_nuc) then
      if (.not. get_scatt) then
        if (allocated(nuc % ndpp_inel)) then
          deallocate(nuc % ndpp_inel)
        end if
      end if
      if (.not. get_nuscatt) then
        if (allocated(nuc % ndpp_nuinel)) then
          deallocate(nuc % ndpp_nuinel)
        end if
      end if
      if ((.not. get_scatt) .and. (.not. get_nuscatt)) then
        if (allocated(nuc % ndpp_el_Ein)) then
          deallocate(nuc % ndpp_el_Ein)
        end if
        if (allocated(nuc % ndpp_el_Ein_srch)) then
          deallocate(nuc % ndpp_el_Ein_srch)
        end if
        if (allocated(nuc % ndpp_inel_Ein)) then
          deallocate(nuc % ndpp_inel_Ein)
        end if
        if (allocated(nuc % ndpp_inel_Ein_srch)) then
          deallocate(nuc % ndpp_inel_Ein_srch)
        end if
        if (allocated(nuc % ndpp_sigma_inel)) then
          deallocate(nuc % ndpp_sigma_inel)
        end if
      end if
      if (.not. get_chi_t) then
        if (allocated(nuc % ndpp_chi)) then
          deallocate(nuc % ndpp_chi)
        end if
      end if
      if (.not. get_chi_p) then
        if (allocated(nuc % ndpp_chi_p)) then
          deallocate(nuc % ndpp_chi_p)
        end if
      end if
      if (.not. get_chi_d) then
        if (allocated(nuc % ndpp_chi_d)) then
          deallocate(nuc % ndpp_chi_d)
        end if
      end if
      if ((.not. get_chi_t) .and. (.not. get_chi_p) .and. (.not. get_chi_d)) then
        if (allocated(nuc % ndpp_chi_Ein)) then
          deallocate(nuc % ndpp_chi_Ein)
        end if
      end if
    else ! Do similar for S(a,b) tables
      if ((.not. get_scatt) .and. (.not. get_nuscatt)) then
        if (allocated(sab % ndpp_el)) then
          deallocate(sab % ndpp_el)
        end if
        if (allocated(sab % ndpp_el_Ein)) then
          deallocate(sab % ndpp_el_Ein)
        end if
        if (allocated(sab % ndpp_el_Ein_srch)) then
          deallocate(sab % ndpp_el_Ein_srch)
        end if
      end if
    end if

  end subroutine read_ndpp_table

end module ndpp
