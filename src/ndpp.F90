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
    type(Nuclide), pointer :: nuc => null() ! Current working nuclide 
    integer :: i_listing    ! index in ndpp_listings array
    integer :: i_nuclide    ! index in nuclides
    type(XsListing), pointer :: ndpp_listing => null()
    
    ! This is probably the best location in the code to do this:
    ! Check to see if S(a,b) tables are in use; if so, print error and quit.
    if (n_sab_tables > 0) then
      message = "Pre-processed scattering kernels do not yet support S(a,b)" // &
        " scattering tables! S(a,b) collisions will be treated with the" // &
        " analog Pn tally method."
      call warning()
    end if
    
    ! First lets go read the ndpp_lib.xml file
    call read_ndpp_xml()
    
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
                  trim(integrated_scatt_lib) // "'!"
        call fatal_error()
      end if
      ! read_ndpp_table will populate nuclide % int_scatt and also check that
      ! the temperatures match
      ndpp_listing => ndpp_listings(i_listing)
      call read_ndpp_table(nuc, ndpp_listing) 
    end do
    
  end subroutine read_ndpp_data

!===============================================================================
! READ_NDPP_XML reads information from a cross_sections.xml file. This
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
    integer :: max_tally_order = 0
    

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
       message = "Unknown filetype in " // trim(integrated_scatt_lib) // &
        ": " // trim(filetype_)
       call fatal_error()
    end if

    ! copy default record length and entries for binary files
    recl = record_length_
    entries = entries_
    
    ! Test metadata to ensure this library matches the problem definition
    ! First test to see if we have a legendre scattering type and not tabular
    !!! (This will be removed at an undefined data when OpenMC supports
    !!! tabular scattering distributions as opposed to Legendres.
    if (scatt_type_ /= SCATT_TYPE_LEGENDRE) then
      message = "Invalid Scattering Type represented in NDPP data. Rerun " // &
                "NDPP with the Legendre scattering type set."
      call fatal_error()
    end if
    
    ! Test to ensure the scattering order is less than the maximum
    if (scatt_order_ > SCATT_ORDER_MAX) then
      message = "Invalid scattering order of " // trim(to_str(scatt_order_)) // &
                " requested."
      call fatal_error()
    end if
    
    ! Check the tallies to ensure that the energy group structure, scattering 
    ! order requested are valid (i.e., groups match, orders are less than in 
    ! the library), and check set all tallies with tracklength estimators and
    ! int-scatter-pn scores to analog so that S(a,b) will be handled.
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
            ! We found the correct score, get comparing!
            ! First check the scattering order
            if (scatt_order_ < t % scatt_order(j)) then
              message = "Invalid scattering order of " // &
                        trim(to_str(scatt_order_)) // " requested. Order " // &
                        "requested is larger than provided in the library (" // &
                        trim(to_str(scatt_order_)) // ")!"
              call fatal_error()
            end if
            ! Find the maximum scattering order requested
            if (t % scatt_order(j) > max_tally_order) &
              max_tally_order = t  % scatt_order(j)
            
            ! Compare the energyin and energyout filters of this tally to the 
            ! energy_bins_ metadata of the NDPP library.
            ! Check the energyin filter first.
            i_filter = t % find_filter(FILTER_ENERGYIN)
            ! We have already checked to ensure some energyin filter exists,
            ! reducing the cases to check. First we check the size, so that we
            ! can use the Fortran intrinsic ALL to check the actual values
            if (t % filters(i_filter) % n_bins /= size(energy_bins_) - 1) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            ! Now we can check the actual group boundaries.
            if (all(t % filters(i_filter) % real_bins /= energy_bins_)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            ! Repeat the same steps as above, but this time for the energyout filter
            i_filter = t % find_filter(FILTER_ENERGYOUT)
            if (t % filters(i_filter) % n_bins /= size(energy_bins_) - 1) then
              message = "Number of groups in NDPP Library do not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            if (all(t % filters(i_filter) % real_bins /= energy_bins_)) then
              message = "NDPP Library group structure does not match that " // &
                        "requested in tally!"
              call fatal_error()
            end if
            
            ! Check to see if the model has S(a,b) data.
            if (n_sab_tables > 0) then
              ! If it does, we can't use tracklength estimators yet with 
              ! int-scatter-pn.  So set the estimator for this tally to analog.
              if (t % estimator == ESTIMATOR_TRACKLENGTH) then
                t % estimator = ESTIMATOR_ANALOG
                message = "Setting the estimator of Tally " // &
                          trim(adjustl(to_str(i))) // " to analog due to the" // &
                          " presence of S(a,b) tables."
                call warning()
              end if
            end if
            
            exit SCORE_LOOP ! we found what we want!
        end select
      end do SCORE_LOOP
    end do TALLY_LOOP

    ! Store the number of groups
    integrated_scatt_groups = size(energy_bins_) - 1
    
    ! Store the order as the maximum requested in tallies
    if (scatt_type_ == SCATT_TYPE_LEGENDRE) then
      integrated_scatt_order = max_tally_order + 1
    else if (scatt_type_ == SCATT_TYPE_TABULAR) then
      ! This one uses scatt_order since there technically is no maximum here
      integrated_scatt_order = scatt_order_ 
    end if

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
  
!===============================================================================
! READ_NDPP_TABLE reads information from a pre-processed nuclear data file
! as generated by the NDPP program.  The result of this subroutine is that the
! Nuclide % int_scatt array will be allocated and filled with the corresponding
! pre-processed scattering (and chi???!!!) data from NDPP
!===============================================================================
  
  subroutine read_ndpp_table(nuc, listing)
    type(Nuclide),   pointer, intent(inout) :: nuc     ! Current nuclide      
    type(XsListing), pointer, intent(in)    :: listing ! Current nuc's NDPP data
    
    integer       :: in = 7        ! file unit
    integer       :: i             ! loop index for data records
    integer       :: location      ! location of NDPP table
    integer       :: filetype      ! Ascii, Binary or HDF5 filetype
    logical       :: file_exists   ! does NDPP library exist?
    character(7)  :: readable      ! is NDPP library readable?
    character(10) :: name          ! name of NDPP table
    real(8)       :: kT            ! temperature of table
    integer       :: NG            ! Number of energy groups in library
    character(MAX_FILE_LEN) :: filename ! path to NDPP data library
    real(8), allocatable    :: energy_bins(:) ! Energy group structure
    integer       :: mu_bins       ! NUmber of angular points used
    integer       :: scatt_type    ! Type of scattering data, discarded
    integer       :: scatt_order   ! Order of scattering data
    integer       :: chi_present   ! Flag as to if chi data is present
    integer       :: gmin, gmax    ! Min and max possible group transfers
    real(8)       :: thin_tol      ! Thinning tolerance used in lib, discarded
    integer       :: NEin, iE      ! Number of incoming energies and the index
    real(8), allocatable :: Ein(:) ! Incoming energies
    real(8), allocatable :: temp_outgoing(:,:) ! Temporary storage of scatt data
    
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
      if ((adjustl(trim(name)) /= listing % name) .or. (kT /= listing % kT)) then 
        message = "NDPP library '" // trim(filename) // "' does not &
                  &contain the correct data set where expected!"
        call fatal_error()
      end if
      ! Get the energy bins (not checking, will assume from here on out we have
      ! the right data set)
      allocate(energy_bins(NG + 1))
      read(UNIT=in, FMT=*) energy_bins
      ! The next line is scatt_type, scatt_order, chi_present, thin_tol, which
      ! we already have. discard.
      read(UNIT=in, FMT='(I20,I20,I20,1PE20.12)') scatt_type, scatt_order, &
        chi_present, thin_tol
      ! Finally, mu_bins, which we can also discard.
      read(UNIT=in, FMT=*)
      
      ! set scatt_order to the right number for allocating the outgoing array
      if (scatt_type == SCATT_TYPE_LEGENDRE) scatt_order = scatt_order + 1
      
      ! Now we get to the meat of the data. 
      ! Start with \chi(E_{in}) data
      if (chi_present == 1) then
        !!! TBI
      else if (chi_present /= 0) then
        message = "NDPP library '" // trim(filename) // "' contains an invalid &
                  & value of chi_present!"
        call fatal_error()
      end if
      
      ! Move to \sigma_{s,g'->g,l}(E_{in}) data
      ! Get Ein information
      read(UNIT=in, FMT=*) NEin
      allocate(Ein(NEin))
      read(UNIT=in, FMT=*) Ein
      !!! Right now the Ein information is the same as in nuc % energy, so
      !!! as for now we do nothing with Ein.
      
      ! Get the moments themselves
      allocate(nuc % int_scatt(NEin))
      do iE = 1, NEin
        ! get gmin and gmax
        read(UNIT=in, FMT=*) gmin, gmax
        nuc % int_scatt(iE) % gmin = gmin
        nuc % int_scatt(iE) % gmax = gmax
                
        if ((gmin > ZERO) .and. (gmax > ZERO)) then
          ! Then we can allocate the space. Do it to integrated_scatt_order
          ! since this is the largest order requested in the tallies.
          ! Since we only need to store up to the maximum, we also need to have
          ! an array for reading the file which we can later truncate to fit
          ! in to nuc % int_scatt(iE) % outgoing.
          allocate(temp_outgoing(scatt_order, gmin : gmax))
          
          allocate(nuc % int_scatt(iE) % outgoing(integrated_scatt_order, &
            gmin : gmax))
          ! Now we have a space to store the data, get it.
          read(UNIT=in, FMT=*) temp_outgoing
          ! And copy in to nuc % int_scatt
          nuc % int_scatt(iE) % outgoing(:, gmin : gmax) = &
            temp_outgoing(1 : integrated_scatt_order, gmin : gmax)
          deallocate(temp_outgoing)
        end if
      end do
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
      if ((adjustl(trim(name)) /= listing % name) .or. (kT /= listing % kT)) then 
        message = "NDPP library '" // trim(filename) // "' does not &
                  &contain the correct data set where expected!"
        call fatal_error()
      end if
      ! Get the energy bins (not checking, will assume from here on out we have
      ! the right data set), and the other meta information
      allocate(energy_bins(NG + 1))
      read(UNIT=in) energy_bins, scatt_type, scatt_order, &
        chi_present, thin_tol, mu_bins
      
      ! set scatt_order to the right number for allocating the outgoing array
      if (scatt_type == SCATT_TYPE_LEGENDRE) scatt_order = scatt_order + 1
      
      ! Now we get to the meat of the data. 
      ! Start with \chi(E_{in}) data
      if (chi_present == 1) then
        !!! TBI
      else if (chi_present /= 0) then
        message = "NDPP library '" // trim(filename) // "' contains an invalid &
                  & value of chi_present!"
        call fatal_error()
      end if
      
      ! Move to \sigma_{s,g'->g,l}(E_{in}) data
      ! Get Ein information
      read(UNIT=in) NEin
      allocate(Ein(NEin))
      read(UNIT=in) Ein
      !!! Right now the Ein information is the same as in nuc % energy, so
      !!! as for now we do nothing with Ein.
      
      ! Get the moments themselves
      allocate(nuc % int_scatt(NEin))
      do iE = 1, NEin
        ! get gmin and gmax
        read(UNIT=in) gmin, gmax
        nuc % int_scatt(iE) % gmin = gmin
        nuc % int_scatt(iE) % gmax = gmax
                
        if ((gmin > ZERO) .and. (gmax > ZERO)) then
          ! Then we can allocate the space. Do it to integrated_scatt_order
          ! since this is the largest order requested in the tallies.
          ! Since we only need to store up to the maximum, we also need to have
          ! an array for reading the file which we can later truncate to fit
          ! in to nuc % int_scatt(iE) % outgoing.
          allocate(temp_outgoing(scatt_order, gmin : gmax))
          
          allocate(nuc % int_scatt(iE) % outgoing(integrated_scatt_order, &
            gmin : gmax))
          ! Now we have a space to store the data, get it.
          read(UNIT=in) temp_outgoing
          ! And copy in to nuc % int_scatt
          nuc % int_scatt(iE) % outgoing(:, gmin : gmax) = &
            temp_outgoing(1 : integrated_scatt_order, gmin : gmax)
          deallocate(temp_outgoing)
        end if
      end do
      close(UNIT=in)
      
    else if (listing % filetype == HDF5) then
      ! =======================================================================
      ! READ NDPP DATA IN HDF5 FORMAT
      !!! TBI
    end if ! No error handling needed; read_ndpp_xml() already checked filetype
    
  end subroutine read_ndpp_table

end module ndpp
