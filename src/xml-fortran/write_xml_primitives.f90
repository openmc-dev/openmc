! write_xml_prims.f90 - Write routines for primitive data
!
! $Id: write_xml_prims.f90,v 1.2 2007/12/27 05:13:59 arjenmarkus Exp $
!
! Arjen Markus
!
! General information:
! This module is part of the XML-Fortran library. Its
! purpose is to write individual items to an XML
! file using the right tag. It is used by the code generated
! by the make_xml_reader program.
!
module write_xml_primitives
    use xmlparse
    implicit none

!   interface write_to_xml
!       module procedure write_to_xml_integers
!       module procedure write_to_xml_reals
!       module procedure write_to_xml_doubles
!       module procedure write_to_xml_logicals
!       module procedure write_to_xml_words
!   end interface
    interface write_to_xml_word
       module procedure write_to_xml_string
    end interface
    interface write_to_xml_line
       module procedure write_to_xml_string
    end interface

contains

! write_to_xml_integer --
!    Routine to write a single integer to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    value       Value to be written
!
subroutine write_to_xml_integer( info, tag, indent, value )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    integer, intent(in)                          :: value

    character(len=100)                           :: indentation

    indentation = ' '
    write( info%lun, '(4a,i0,3a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>', value, '</', trim(tag), '>'

end subroutine write_to_xml_integer

! write_to_xml_integer_1dim --
!    Routine to write an array of integers to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    values      Values to be written
!
subroutine write_to_xml_integer_1dim( info, tag, indent, values )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    integer, dimension(:), intent(in)            :: values

    integer                                      :: i

    do i = 1,size(values)
        call write_to_xml_integer( info, tag, indent, values(i) )
    enddo

end subroutine write_to_xml_integer_1dim

! write_to_xml_real --
!    Routine to write a single real value (single precision) to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    value       Value to be written
!
subroutine write_to_xml_real( info, tag, indent, value )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    real, intent(in)                             :: value

    character(len=100)                           :: indentation
    character(len=12)                            :: buffer

    indentation = ' '
    write( buffer, '(1pg12.4)' ) value
    write( info%lun, '(8a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>', trim(adjustl(buffer)), '</', trim(tag), '>'

end subroutine write_to_xml_real

! write_to_xml_real_1dim --
!    Routine to write an array of reals to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    values      Values to be written
!
subroutine write_to_xml_real_1dim( info, tag, indent, values )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    real, dimension(:), intent(in)               :: values

    integer                                      :: i

    do i = 1,size(values)
        call write_to_xml_real( info, tag, indent, values(i) )
    enddo

end subroutine write_to_xml_real_1dim

! write_to_xml_double --
!    Routine to write one real value (double precision) to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    value       Value to be written
!
subroutine write_to_xml_double( info, tag, indent, value )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    real(kind=kind(1.0d0)), intent(in)           :: value

    character(len=100)                           :: indentation
    character(len=16)                            :: buffer

    indentation = ' '
    write( buffer, '(1pg16.7)' ) value
    write( info%lun, '(8a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>', trim(adjustl(buffer)), '</', trim(tag), '>'

end subroutine write_to_xml_double

! write_to_xml_double_1dim --
!    Routine to write an array of double precision reals to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    values      Values to be written
!
subroutine write_to_xml_double_1dim( info, tag, indent, values )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    real(kind=kind(1.0d00)), dimension(:), intent(in) :: values

    integer                                      :: i

    do i = 1,size(values)
        call write_to_xml_double( info, tag, indent, values(i) )
    enddo

end subroutine write_to_xml_double_1dim

! write_to_xml_string --
!    Routine to write one string to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    value       Value to be written
!
subroutine write_to_xml_string( info, tag, indent, value )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    character(len=*), intent(in)                 :: value

    character(len=100)                           :: indentation

    !
    ! NOTE: No guards against <, >, & and " yet!
    ! NOTE: difference needed between words and lines?
    !
    indentation = ' '
    write( info%lun, '(8a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>', trim(value), '</', trim(tag), '>'

end subroutine write_to_xml_string

! write_to_xml_word_1dim --
!    Routine to write an array of single words to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    value       Value to be written
!
subroutine write_to_xml_word_1dim( info, tag, indent, values )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    character(len=*), dimension(:), intent(in)   :: values

    integer                                      :: i

    do i = 1,size(values)
        call write_to_xml_string( info, tag, indent, values(i) )
    enddo
end subroutine write_to_xml_word_1dim

! write_to_xml_string_1dim --
!    Routine to write an array of strings to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    values      Values to be written
!
subroutine write_to_xml_string_1dim( info, tag, indent, values )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    character(len=*), dimension(:), intent(in)   :: values

    integer                                      :: i

    do i = 1,size(values)
        call write_to_xml_string( info, tag, indent, values(i) )
    enddo

end subroutine write_to_xml_string_1dim

! write_to_xml_logical --
!    Routine to write one logical to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    value       Value to be written
!
subroutine write_to_xml_logical( info, tag, indent, value )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    logical, intent(in)                          :: value

    character(len=100)                           :: indentation

    indentation = ' '
    if ( value ) then
        write( info%lun, '(8a)' ) indentation(1:min(indent,100)), &
            '<', trim(tag), '>true</', trim(tag), '>'
    else
        write( info%lun, '(8a)' ) indentation(1:min(indent,100)), &
            '<', trim(tag), '>false</', trim(tag), '>'
    endif

end subroutine write_to_xml_logical

! write_to_xml_logical_1dim --
!    Routine to write an array of logicals to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    values      Values to be written
!
subroutine write_to_xml_logical_1dim( info, tag, indent, values )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    logical, dimension(:), intent(in)            :: values

    integer                                      :: i

    do i = 1,size(values)
        call write_to_xml_logical( info, tag, indent, values(i) )
    enddo

end subroutine write_to_xml_logical_1dim

! write_to_xml_integer_array --
!    Routine to write an array of integers to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    array       Values to be written
!
subroutine write_to_xml_integer_array( info, tag, indent, array )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    integer, dimension(:), intent(in)            :: array

    character(len=100)                           :: indentation
    integer                                      :: i, i2, j

    indentation = ' '

    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>'
    do i = 1,size(array),10
        i2 = min( i + 9, size(array) )
        write( info%lun, '(a,10i12)' ) indentation(1:min(indent+4,100)), &
            ( array(j) ,j = i,i2 )
    enddo
    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '</', trim(tag), '>'

end subroutine write_to_xml_integer_array

! write_to_xml_real_array --
!    Routine to write an array of single precision reals to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    array       Values to be written
!
subroutine write_to_xml_real_array( info, tag, indent, array )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    real, dimension(:), intent(in)               :: array

    character(len=100)                           :: indentation
    integer                                      :: i, i2, j

    indentation = ' '

    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>'
    do i = 1,size(array),10
        i2 = min( i + 9, size(array) )
        write( info%lun, '(a,10g12.4)' ) indentation(1:min(indent+4,100)), &
            ( array(j) ,j = i,i2 )
    enddo
    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '</', trim(tag), '>'

end subroutine write_to_xml_real_array

! write_to_xml_double_array --
!    Routine to write an array of double precision reals to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    array       Values to be written
!
subroutine write_to_xml_double_array( info, tag, indent, array )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    real(kind=kind(1.0d0)), dimension(:), intent(in) :: array

    character(len=100)                           :: indentation
    integer                                      :: i, i2, j

    indentation = ' '

    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>'
    do i = 1,size(array),5
        i2 = min( i + 4, size(array) )
        write( info%lun, '(a,5g20.7)' ) indentation(1:min(indent+4,100)), &
            ( array(j) ,j = i,i2 )
    enddo
    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '</', trim(tag), '>'

end subroutine write_to_xml_double_array

! write_to_xml_logical_array --
!    Routine to write an array of logicals to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    array       Values to be written
!
subroutine write_to_xml_logical_array( info, tag, indent, array )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    logical, dimension(:), intent(in)            :: array

    character(len=100)                           :: indentation
    integer                                      :: i, i2, j

    indentation = ' '

    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>'
    do i = 1,size(array),10
        i2 = min( i + 9, size(array) )
        write( info%lun, '(a,10a)' ) indentation(1:min(indent+4,100)), &
            ( merge('true  ', 'false ', array(j)) ,j = i,i2 )
    enddo
    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '</', trim(tag), '>'

end subroutine write_to_xml_logical_array

! write_to_xml_word_array --
!    Routine to write an array of words to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    array       Values to be written
!
subroutine write_to_xml_word_array( info, tag, indent, array )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    character(len=*), dimension(:), intent(in)   :: array

    character(len=100)                           :: indentation
    integer                                      :: i, i2, j

    indentation = ' '

    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>'
    do i = 1,size(array),10
        i2 = min( i + 9, size(array) )
        write( info%lun, '(a,20a)' ) indentation(1:min(indent+4,100)), &
            ( trim(array(j)) , ' ' ,j = i,i2 )
    enddo
    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '</', trim(tag), '>'

end subroutine write_to_xml_word_array

! write_to_xml_line_array --
!    Routine to write an array of lines to the XML file
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question
!    indent      Number of spaces for indentation
!    array       Values to be written
!
subroutine write_to_xml_line_array( info, tag, indent, array )
    type(XML_PARSE), intent(in)                  :: info
    character(len=*), intent(in)                 :: tag
    integer, intent(in)                          :: indent
    logical, dimension(:), intent(in)            :: array

    character(len=100)                           :: indentation
    integer                                      :: i, i2, j

    indentation = ' '

    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '<', trim(tag), '>'
    do i = 1,size(array)
        write( info%lun, '(a)' ) indentation(1:min(indent+4,100)), &
            array(i)
    enddo
    write( info%lun, '(4a)' ) indentation(1:min(indent,100)), &
        '</', trim(tag), '>'

end subroutine write_to_xml_line_array

end module write_xml_primitives
