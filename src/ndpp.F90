module ndpp

  use ace_header,       only: Nuclide, GrpTransfer, XsListing
  use constants
  use dict_header,      only: DictCharInt, DICT_NULL
  use error,            only: fatal_error, warning
  use global
  use output,           only: write_message
  use string,           only: ends_with, lower_case, starts_with, to_str

  implicit none
    
contains

!===============================================================================
! READ_NDPP_LIB reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_ndpp_data()
    type(Nuclide), pointer :: nuc => null()
    character(12)  :: name  ! name of isotope, e.g. 92235.03c
    character(12)  :: alias ! alias of nuclide, e.g. U-235.03c
    integer :: i_listing    ! index in ndpp_listings array
    integer :: i_nuclide    ! index in nuclides
    ! Parse through each nuclide in the model and read in the corresponding
    ! NDPP data.
    ! During this we will make sure the temperatures of ACE and NDPP match,
    ! as well as making sure the NDPP data exists at all.
    
    do i_nuclide = 1, n_nuclides_total
      nuc => nuclides(i_nuclide)
      i_listing = ndpp_listing_dict % get_key(nuc % name)
      if (i_listing == DICT_NULL) then
        ! Could not find ndpp_lib.xml file
        message = nuc % name // "does not exist in NDPP XML file: '" //  &
                  trim(integrated_scatt_lib) // "!"
        call fatal_error()
      end if
      ! read_ndpp_table will populate nuclide % int_scatt and also check that
      ! the temperatures match
!~       call read_ndpp_table(i_nuclide, i_listing)
            
    end do
    
  end subroutine read_ndpp_data


!===============================================================================
! READ_CROSS_SECTIONS_XML reads information from a cross_sections.xml file. This
! file contains a listing of the ACE cross sections that may be used.
!===============================================================================

  subroutine read_ndpp_xml()

    use xml_data_ndpp_lib_t

    integer :: i, j, k     ! loop indices
    integer :: filetype    ! default file type
    integer :: recl        ! default record length
    integer :: entries     ! default number of entries
    logical :: file_exists ! does cross_sections.xml exist?
    character(MAX_WORD_LEN)  :: directory   ! directory with cross sections
    type(TallyObject), pointer :: t => null()
    integer :: i_filter    ! Index of filter which contains energyin or energyout
    ! We can use the same XSListing type for our ndpp data since the NDPP
    ! data is a subset of whats in cross_sections.xml
    type(XsListing), pointer :: listing => null()
    

    ! Check if cross_sections.xml exists
    inquire(FILE=integrated_scatt_lib, EXIST=file_exists)
    if (.not. file_exists) then
       ! Could not find ndpp_lib.xml file
       message = "NDPP XML file '" // trim(integrated_scatt_lib) // &
            "' does not exist!"
       call fatal_error()
    end if
    
    message = "Reading NDPP Library XML file..."
    call write_message(5)

    ! Initialize variables that may go unused
    directory_ = ""
    filetype_ = ""
    record_length_ = 0
    entries_ = 0

    ! Parse cross_sections.xml file
    call read_xml_file_ndpp_lib_t(integrated_scatt_lib)

    if (len_trim(directory_) > 0) then
       ! Copy directory information if present
       directory = trim(directory_)
    else
       ! If no directory is listed in the NDPP xml, by default select the
       ! directory in which the NDPP xml file resides
       i = index(path_cross_sections, "/", BACK=.true.)
       directory = integrated_scatt_lib(1:i)
    end if

    ! determine whether binary, ascii, or hdf5
    call lower_case(filetype_)
    if (filetype_ == 'ascii') then
       filetype = ASCII
    elseif (filetype_ == 'binary') then
       filetype = BINARY
    elseif (filetype_ == 'hdf5') then
       filetype = HDF5
    else
       message = "Unknown filetype in ndpp_library.xml: " // trim(filetype_)
       call fatal_error()
    end if

    ! copy default record length and entries for binary files
    recl = record_length_
    entries = entries_
    
    ! Test metadata to ensure this library matches the problem definition
    ! First test to see if we have a legendre scattering type and not tabular
    ! (This will be removed at an undefined data when OpenMC supports
    ! tabular scattering distributions as opposed to Legendres.
    if (scatt_type_ /= SCATT_TYPE_LEGENDRE) then
      message = "Invalid Scattering Type represented in NDPP data. Rerun " // &
                "NDPP with the Legendre scattering type set."
      call fatal_error()
    end if
    
    ! Test to ensure the scattering order is less than the maximum
    if (scatt_order_ > SCATT_ORDER_MAX) then
      message = "Invalid scattering order of " // trim(to_str(scatt_order_)) // &
                " requested. Setting to the maximum permissible value, " // &
                trim(to_str(SCATT_ORDER_MAX))
      call fatal_error()
    end if
    
    ! Test that the energy group structure matches that requested in 
    ! tallies.xml.
    TALLY_LOOP: do i = 1, n_tallies
      t => tallies(i)
      j = 0
      SCORE_LOOP: do k = 1, t % n_user_score_bins
        j = j + 1
        select case (t % score_bins(j))
          case (SCORE_SCATTER_PN)
            j = j + t % scatt_order(j)
            cycle SCORE_LOOP ! Skip the others which will only waste cycles
          case (SCORE_INTSCATT_PN)
            ! We found a tally with the right kind of score, compare the
            ! energyin and energyout filters of this tally to the energy_bins_
            ! metadata of the NDPP library.
            ! Check the energyin filter first.
            i_filter = t % find_filter(FILTER_ENERGYIN)
            ! We have already checked to ensure some energyin filter exists,
            ! reducing the cases to check. First we check the size, so that we
            ! can use the Fortran intrinsic ALL to check the actual values
            if (t % filters(i_filter) % n_bins /= size(energy_bins_)) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!. "
              call fatal_error()
            end if
            ! Now we can check the actual group boundaries.
            if (all(t % filters(i_filter) % real_bins /= energy_bins_)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!."
              call fatal_error()
            end if
            ! Repeat the same steps as above, but this time for the energyout filter
            i_filter = t % find_filter(FILTER_ENERGYOUT)
            if (t % filters(i_filter) % n_bins /= size(energy_bins_)) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!. "
              call fatal_error()
            end if
            if (all(t % filters(i_filter) % real_bins /= energy_bins_)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!."
              call fatal_error()
            end if
            
            exit SCORE_LOOP !we found what we want!
        end select
      end do SCORE_LOOP
    end do TALLY_LOOP

    ! Allocate ndpp_listings array
    if (.not. associated(ndpp_tables_)) then
       message = "No NDPP table listings present in ndpp_lib.xml file!"
       call fatal_error()
    else
       n_listings = size(ndpp_tables_)
       allocate(ndpp_listings(n_listings))
    end if

    do i = 1, n_listings
       listing => ndpp_listings(i)

       ! copy a number of attributes
       listing % name       = trim(ndpp_tables_(i) % name)
       listing % alias      = trim(ndpp_tables_(i) % alias)
       listing % zaid       = ndpp_tables_(i) % zaid
       listing % awr        = ndpp_tables_(i) % awr
       listing % kT         = ndpp_tables_(i) % temperature
       listing % location   = ndpp_tables_(i) % location

       ! determine type of cross section
       if (ends_with(listing % name, 'c')) then
          listing % type = ACE_NEUTRON
       elseif (ends_with(listing % name, 't')) then
          listing % type = ACE_THERMAL
       end if

       ! set filetype, record length, and number of entries
       listing % filetype = filetype
       listing % recl     = recl
       listing % entries  = entries

       ! determine metastable state
       if (ndpp_tables_(i) % metastable == 0) then
          listing % metastable = .false.
       else
          listing % metastable = .true.
Title	T	P	Status	Assignee	Created	Last updated Actions
#
       end if

       ! determine path of cross section table
       if (starts_with(ndpp_tables_(i) % path, '/')) then
          listing % path = ndpp_tables_(i) % path
       else
          if (ends_with(directory,'/')) then
             listing % path = trim(directory) // trim(ndpp_tables_(i) % path)
          else
             listing % path = trim(directory) // '/' // &
              trim(ndpp_tables_(i) % path)
          end if
       end if

       ! create dictionary entry for both name and alias
       call ndpp_listing_dict % add_key(listing % name, i)
       call ndpp_listing_dict % add_key(listing % alias, i)
    end do

  end subroutine read_ndpp_xml

end module ndpp
