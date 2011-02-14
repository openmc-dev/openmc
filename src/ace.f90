module ace

  use global
  use output, only: error, message
  use string, only: lower_case
  use fileio, only: read_line, read_data, skip_lines
  use string, only: split_string, str_to_real
  use data_structures, only: dict_create, dict_add_key, dict_has_key, &
       &                     dict_get_key, dict_delete

  integer :: NXS(16)
  integer :: JXS(32)
  real(8), allocatable :: XSS(:)
  integer :: XSS_index

  private :: NXS
  private :: JXS
  private :: XSS

contains

!=====================================================================
! READ_XS reads all the cross sections for the problem and stores them
! in xs_continuous and xs_thermal arrays
!=====================================================================

  subroutine read_xs()

    type(Material),      pointer :: mat => null()
    type(xsData),        pointer :: iso => null()
    type(AceContinuous), pointer :: ace_cont => null()
    type(AceThermal),    pointer :: ace_thermal => null()
    integer :: i, j
    integer :: index
    character(10) :: key
    character(250) :: msg
    integer :: n
    integer :: index_continuous
    integer :: index_thermal
    type(DictionaryCI),   pointer :: temp_dict

    call dict_create(ace_dict)

    ! determine how many continuous-energy tables and how many S(a,b)
    ! thermal scattering tables there are
    index_continuous = 0
    index_thermal = 0
    do i = 1, n_materials
       mat => materials(i)
       do j = 1, mat%n_isotopes
          index = mat%isotopes(j)
          key = xsdatas(index)%id
          n = len_trim(key)
          call lower_case(key)
          select case (key(n:n))
          case ('c')
             if (.not. dict_has_key(ace_dict, key)) then
                index_continuous = index_continuous + 1
                call dict_add_key(ace_dict, key, index_continuous)
                mat%table(j) = index_continuous
             else
                mat%table(j) = dict_get_key(ace_dict, key)
             end if
          case ('t')
             n_thermal = n_thermal + 1
          case default
             msg = "Unknown cross section table type: " // key
             call error(msg)
          end select
       end do
    end do

    n_continuous = index_continuous
    n_thermal = index_thermal

    ! allocate arrays for ACE table storage
    allocate(xs_continuous(n_continuous))
    allocate(xs_thermal(n_thermal))

    ! loop over all nuclides in xsdata
    call dict_create(temp_dict)

    index_continuous = 0
    index_thermal = 0
    do i = 1, n_materials
       mat => materials(i)
       do j = 1, mat%n_isotopes
          index = mat%isotopes(j)
          key = xsdatas(index)%id
          n = len_trim(key)
          call lower_case(key)
          select case (key(n:n))
          case ('c')
             if (.not. dict_has_key(temp_dict, key)) then
                index_continuous = index_continuous + 1
                call read_ACE_continuous(index_continuous, index)
             end if
          case ('t')
             n_thermal = n_thermal + 1
          end select
       end do
    end do

    ! delete dictionary
    call dict_delete(temp_dict)
       
  end subroutine read_xs

!=====================================================================
! READ_ACE_CONTINUOUS reads in a single ACE continuous-energy cross
! section table
!=====================================================================

  subroutine read_ACE_continuous(index_table, index)

    integer, intent(in) :: index_table
    integer, intent(in) :: index

    type(AceContinuous), pointer :: table => null()
    integer :: i
    integer :: in = 7
    integer :: ioError
    integer :: words_per_line
    integer :: lines
    integer :: n
    logical :: file_exists
    logical :: found_xs
    character(7) :: readable
    character(250) :: msg, line
    character(32) :: words(max_words)
    character(100) :: filename
    character(10) :: tablename
    real(8) :: kT

    filename = xsdatas(index)%path
    tablename = xsdatas(index)%id

    table => xs_continuous(index_table)

    ! Check if input file exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
       msg = "ACE library '" // trim(filename) // "' does not exist!"
       call error(msg)
    elseif (readable(1:3) == 'NO') then
       msg = "ACE library '" // trim(filename) // "' is not readable! &
            &Change file permissions with chmod command."
       call error(msg)
    end if

    ! display message
    msg = "Loading ACE cross section table: " // tablename
    call message(msg, 6)

    ! open file
    open(file=filename, unit=in, status='old', & 
         & action='read', iostat=ioError)
    if (ioError /= 0) then
       msg = "Error while opening file: " // filename
       call error(msg)
    end if

    found_xs = .false.
    do while (.not. found_xs)
       call read_line(in, line, ioError)
       if (ioError < 0) then
          msg = "Could not find ACE table " // tablename // "."
          call error(msg)
       end if
       call split_string(line, words, n)
       if (trim(words(1)) == trim(tablename)) then
          found_xs = .true.
          table%name = words(1)
          table%awr = str_to_real(words(2))
          kT = str_to_real(words(3))
          table%temp = kT / K_BOLTZMANN
       end if
       
       ! Skip 5 lines
       call skip_lines(in, 5, ioError)

       ! Read NXS data
       lines = 2
       words_per_line = 8
       call read_data(in, NXS, 16, lines, words_per_line)

       ! Read JXS data
       lines = 4
       call read_data(in, JXS, 32, lines, words_per_line)

       ! Calculate how many data points and lines in the XSS array
       n = NXS(1)
       lines = (n + 3)/4

       if (found_xs) then
          ! allocate storage for XSS array
          allocate(XSS(n))
          
          ! Read XSS
          words_per_line = 4
          call read_data(in, XSS, n, lines, words_per_line)
       else
          call skip_lines(in, lines, ioError)
       end if

    end do

    call parse_ESZ(table)
    call parse_MTR(table)

    ! Free memory from XSS array
    if(allocated(XSS)) deallocate(XSS)
    if(associated(table)) nullify(table)

    close(unit=in)

  end subroutine read_ACE_continuous

!=====================================================================
! PARSE_ESZ - reads through the ESZ block. This block contains the
! energy grid, total xs, absorption xs, elastic scattering xs, and
! heating numbers.
!=====================================================================

  subroutine parse_ESZ(table)

    type(AceContinuous), pointer, intent(inout) :: table

    integer :: NE

    ! determine number of energy points
    NE = NXS(3)
    table%n_grid = NE

    ! allocate storage for arrays
    allocate(table%energy(NE))
    allocate(table%sigma_t(NE))
    allocate(table%sigma_a(NE))
    allocate(table%sigma_el(NE))
    allocate(table%heating(NE))

    ! read data from XSS -- right now the total, absorption and
    ! elastic scattering are read in to these special arrays, but in
    ! reality, it should be necessary to only store elastic scattering
    ! and possible total cross-section for total material xs
    ! generation.
    XSS_index = 1
    table%energy = get_real(NE)
    table%sigma_t = get_real(NE)
    table%sigma_a = get_real(NE)
    table%sigma_el = get_real(NE)
    table%heating = get_real(NE)
    
  end subroutine parse_ESZ

!=====================================================================
! PARSE_MTR - Get the list of reaction MTs for this cross-section
! table. The MT values are somewhat arbitrary. Also read in Q-values,
! neutron multiplicities, and cross-sections.
! =====================================================================

  subroutine parse_MTR(table)

    type(AceContinuous), pointer :: table

    type(AceReaction), pointer :: rxn
    integer :: LMT   ! index of MT list in XSS
    integer :: NMT   ! Number of reactions
    integer :: JXS4  ! index of Q values in XSS
    integer :: JXS5  ! index of neutron multiplicities in XSS
    integer :: JXS7  ! index of reactions cross-sections in XSS
    integer :: LXS   ! 
    integer :: LOCA  ! 
    integer :: NE    ! number of energies for reaction
    
    LMT = JXS(3)
    JXS4 = JXS(4)
    JXS5 = JXS(5)
    LXS  = JXS(6)
    JXS7 = JXS(7)
    NMT = NXS(4)

    ! allocate array of reactions. Add one since we need to include an
    ! elastic scattering channel
    table%n_reaction = NMT + 1
    allocate(table%reactions(NMT+1))

    ! Store elastic scattering cross-section on reaction one
    rxn => table%reactions(1)
    rxn%MT      = 2
    rxn%Q_value = 0.0_8
    rxn%TY      = 1
    rxn%IE      = 1
    allocate(rxn%sigma(table%n_grid))
    rxn%sigma = table%sigma_el

    do i = 1, NMT
       rxn => table%reactions(i+1)

       ! read MT number, Q-value, and neutrons produced
       rxn%MT      = XSS(LMT+i-1)
       rxn%Q_value = XSS(JXS4+i-1)
       rxn%TY      = XSS(JXS5+i-1)

       ! read cross section values
       LOCA = XSS(LXS+i-1)
       rxn%IE = XSS(JXS7 + LOCA - 1)
       NE = XSS(JXS7 + LOCA)
       allocate(rxn%sigma(NE))
       XSS_index = JXS7 + LOCA + 1
       rxn%sigma = get_real(NE)
    end do

  end subroutine parse_MTR

!=====================================================================
! GET_INT returns an array of integers read from the current position
! in the XSS array
!=====================================================================

    function get_int(n_values) result(array)

      integer, intent(in) :: n_values
      integer :: array(n_values)

      array = int(XSS(XSS_index:XSS_index + n_values - 1))
      XSS_index = XSS_index + n_values
      
    end function get_int

!=====================================================================
! GET_REAL returns an array of real(8)s read from the current position
! in the XSS array
!=====================================================================

    function get_real(n_values) result(array)

      integer, intent(in) :: n_values
      real(8) :: array(n_values)

      array = XSS(XSS_index:XSS_index + n_values - 1)
      XSS_index = XSS_index + n_values

    end function get_real

!=====================================================================
! READ_ACE_THERMAL reads in a single ACE S(a,b) thermal scattering
! cross section table
!=====================================================================

  subroutine read_ACE_thermal(filename)

    character(*), intent(in) :: filename

  end subroutine read_ACE_thermal

end module ace
