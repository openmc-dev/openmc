module cross_section

  use string, only: split_string, str_to_real
  use global, only: str_to_int, max_words, xsdata_dict, xsdatas
  use data_structures, only: Dictionary, dict_create, dict_add_key, dict_get_key
  use output, only: error
  use types, only: xsData

  implicit none

contains

!=====================================================================
! READ_XSDATA reads the data in a SERPENT xsdata file and builds a
! dictionary to find cross-section information later on.
!=====================================================================

  subroutine read_xsdata( path )

    character(*), intent(in) :: path

    type(xsData), pointer :: iso => null()
    character(250) :: line, msg
    character(100) :: words(max_words)
    character(100) :: filename
    integer :: i, n
    integer :: in = 7
    logical :: file_exists
    character(7) :: readable
    integer :: count
    integer :: index
    integer :: ioError

    ! Construct filename
    n = len(path)
    if (path(n:n) == '/') then
       filename = trim(path) // 'xsdata'
    else
       filename = trim(path) // '/xsdata'
    end if

    ! Check if xsdata exists and is readable
    inquire( FILE=filename, EXIST=file_exists, READ=readable )
    if ( .not. file_exists ) then
       msg = "Cross section summary '" // trim(filename) // "' does not exist!"
       call error( msg )
    elseif ( readable(1:3) /= 'YES' ) then
       msg = "Cross section summary '" // trim(filename) // "' is not readable!" &
            & // "Change file permissions with chmod command."
       call error( msg )
    end if

    ! open file
    open(file=filename, unit=in, status='old', action='read')

    ! determine how many lines
    count = 0
    do
       read(unit=in, fmt='(A250)', iostat=ioError) line
       if ( ioError < 0 ) then
          ! reached end of file
          exit
       elseif ( ioError > 0 ) then
          msg = "Unknown error while reading file: " // filename
          close(unit=in)
          call error( msg )
       end if
       count = count + 1
    end do
    allocate(xsdatas(count))

    ! create dictionary
    call dict_create(xsdata_dict)

    ! read actual lines
    index = 0
    rewind(in)
    do
       read(unit=in, fmt='(A250)', iostat=ioError) line
       if ( ioError < 0 ) exit
       index = index + 1
       call split_string(line, words, n)
       if ( n == 0 ) cycle ! skip blank line

       ! Check to make sure there are enough arguments
       if ( n < 9 ) then
          msg = "Not enough arguments on xsdata line: " // line
          close(unit=in)
          call error( msg )
       end if

       iso => xsdatas(index)

       ! store data
       iso%alias = words(1)
       iso%id = words(2)
       iso%type = str_to_int(words(3))
       iso%zaid = str_to_int(words(4))
       iso%isomeric = str_to_int(words(5))
       iso%awr = str_to_real(words(6))
       iso%temp = str_to_real(words(7))
       iso%binary = str_to_int(words(8))
       iso%path = words(9)

       ! create dictionary entry
       call dict_add_key(xsdata_dict, iso%alias, index)
    end do

    close(unit=in)

  end subroutine read_xsdata

end module cross_section
