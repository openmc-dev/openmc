! read_xml_prims.f90 - Read routines for primitive data
!
! $Id: read_xml_prims.f90,v 1.7 2007/12/07 10:38:41 arjenmarkus Exp $
!
! Arjen Markus
!
! General information:
! This module is part of the XML-Fortran library. Its
! purpose is to help read individual items from an XML
! file into the variables that have been connected to
! the various tags. It is used by the code generated
! by the make_xml_reader program.
!
! Because the routines differ mostly by the type of the
! output variable, the body is included, to prevent
! too much repeated blocks of code with all the maintenance
! issues that causes.
!
module read_xml_primitives
   use xmlparse
   implicit none

   private :: read_from_buffer
   private :: read_from_buffer_integers
   private :: read_from_buffer_reals
   private :: read_from_buffer_doubles
   private :: read_from_buffer_logicals
   private :: read_from_buffer_words

   interface read_from_buffer
      module procedure read_from_buffer_integers
      module procedure read_from_buffer_reals
      module procedure read_from_buffer_doubles
      module procedure read_from_buffer_logicals
      module procedure read_from_buffer_words
   end interface

contains

! skip_until_endtag --
!    Routine to read the XML file until the end tag is encountered
!
! Arguments:
!    info        The XML file data structure
!    tag         The tag in question
!    attribs     Array of attributes and their values
!    data        Array of strings, representing the data
!    error       Has an error occurred?
!
subroutine skip_until_endtag( info, tag, attribs, data, error )
   type(XML_PARSE), intent(inout)                  :: info
   character(len=*), intent(in)                    :: tag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   character(len=*), dimension(:), intent(inout)   :: data
   logical, intent(out)                            :: error

   integer                                         :: noattribs
   integer                                         :: nodata
   integer                                         :: ierr
   logical                                         :: endtag
   character(len=len(tag))                         :: newtag

   error = .true.
   do
      call xml_get( info, newtag, endtag, attribs, noattribs, &
                    data, nodata )
      if ( xml_error(info) ) then
         error = .true.
         exit
      endif
      if ( endtag .and. newtag == tag ) then
         exit
      endif
   enddo
end subroutine skip_until_endtag

! read_xml_integer --
!    Routine to read a single integer from the parsed data
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question (error message only)
!    endtag      End tag found? (Dummy argument, actually)
!    attribs     Array of attributes and their values
!    noattribs   Number of attributes found
!    data        Array of strings, representing the data
!    nodata      Number of data strings
!    var         Variable to be filled
!    has_var     Has the variable been set?
!
subroutine read_xml_integer( info, tag, endtag, attribs, noattribs, data, nodata, &
                             var, has_var )
   integer, intent(inout)                       :: var

   include 'read_xml_scalar.inc'

end subroutine read_xml_integer

! read_xml_line --
!    Routine to read a single line of text from the parsed data
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question (error message only)
!    endtag      End tag found? (Dummy argument, actually)
!    attribs     Array of attributes and their values
!    noattribs   Number of attributes found
!    data        Array of strings, representing the data
!    nodata      Number of data strings
!    var         Variable to be filled
!    has_var     Has the variable been set?
!
subroutine read_xml_line( info, tag, endtag, attribs, noattribs, data, nodata, &
                          var, has_var )
   type(XML_PARSE), intent(inout)               :: info
   character(len=*), intent(in)                 :: tag
   logical, intent(inout)                       :: endtag
   character(len=*), dimension(:,:), intent(in) :: attribs
   integer, intent(in)                          :: noattribs
   character(len=*), dimension(:), intent(in)   :: data
   integer, intent(in)                          :: nodata
   character(len=*), intent(inout)              :: var
   logical, intent(inout)                       :: has_var

   character(len=len(attribs(1,1)))             :: buffer
   integer                                      :: idx
   integer                                      :: ierr

   !
   ! The value can be stored in an attribute value="..." or in
   ! the data
   !
   has_var = .false.
   idx = xml_find_attrib( attribs, noattribs, 'value', buffer )
   if ( idx > 0 ) then
      var     = buffer
      has_var = .true.
   else
      do idx = 1,nodata
         if ( data(idx) /= ' ' ) then
            var = data(idx)
            has_var = .true.
            exit
         endif
      enddo
   endif
end subroutine read_xml_line

! read_xml_real, ... --
!    See read_xml_integer for an explanation
!
subroutine read_xml_real( info, tag, endtag, attribs, noattribs, data, nodata, &
                          var, has_var )
   real, intent(inout)                          :: var

   include 'read_xml_scalar.inc'

end subroutine read_xml_real

subroutine read_xml_double( info, tag, endtag, attribs, noattribs, data, nodata, &
                            var, has_var )
   real(kind=kind(1.0d00)), intent(inout)       :: var

   include 'read_xml_scalar.inc'

end subroutine read_xml_double

subroutine read_xml_logical( info, tag, endtag, attribs, noattribs, data, nodata, &
                             var, has_var )
   logical, intent(inout)       :: var

   include 'read_xml_scalar.inc'

end subroutine read_xml_logical

subroutine read_xml_word( info, tag, endtag, attribs, noattribs, data, nodata, &
                          var, has_var )
   character(len=*), intent(inout)       :: var

   include 'read_xml_word.inc'

end subroutine read_xml_word

! read_xml_integer_array --
!    Routine to read a one-dimensional integer array from the parsed
!    ata
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question (error message only)
!    endtag      End tag found? (Dummy argument, actually)
!    attribs     Array of attributes and their values
!    noattribs   Number of attributes found
!    data        Array of strings, representing the data
!    nodata      Number of data strings
!    var         Variable to be filled
!    has_var     Has the variable been set?
!
subroutine read_xml_integer_array( info, tag, endtag, attribs, noattribs, data, &
                                   nodata, var, has_var )
   integer, dimension(:), pointer                :: var

   include 'read_xml_array.inc'

end subroutine read_xml_integer_array

! read_xml_line_array --
!    Routine to read an array of lines of text from the parsed data
!
! Arguments:
!    info        XML parser structure
!    tag         The tag in question (error message only)
!    attribs     Array of attributes and their values
!    noattribs   Number of attributes found
!    data        Array of strings, representing the data
!    nodata      Number of data strings
!    var         Variable to be filled
!    has_var     Has the variable been set?
!
subroutine read_xml_line_array( info, tag, endtag, attribs, noattribs, data, &
                                nodata, var, has_var )
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   character(len=*), dimension(:), pointer       :: var
   logical, intent(inout)                        :: has_var

   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: idx
   integer                                       :: idxv
   integer                                       :: ierr
   logical                                       :: started

   !
   ! The value can be stored in an attribute values="..." or in
   ! the data
   !
   has_var = .false.
   idx = xml_find_attrib( attribs, noattribs, 'values', buffer )
   if ( idx > 0 ) then
      allocate( var(1:1) )
      var(1) = buffer
      if ( buffer /= ' ' ) then
         has_var = .true.
      endif
   else
      idxv    = 0
      started = .false.
      do idx = 1,nodata
         if ( data(idx) /= ' ' .or. started ) then
            if ( .not. started ) then
               allocate( var(1:nodata-idx+1) )
               started = .true.
            endif
            idxv = idxv + 1
            var(idxv) = data(idx)
         endif
      enddo
      if ( started ) then
         has_var = .true.
      endif
   endif
end subroutine read_xml_line_array

! read_xml_real_array, ... --
!    See read_xml_integer_array for an explanation
!
subroutine read_xml_real_array( info, tag, endtag, attribs, noattribs, data, &
                                nodata, var, has_var )
   real, dimension(:), pointer :: var

   include 'read_xml_array.inc'

end subroutine read_xml_real_array

subroutine read_xml_double_array( info, tag, endtag, attribs, noattribs, data, &
                                  nodata, var, has_var )
   real(kind=kind(1.0d00)), dimension(:), pointer :: var

   include 'read_xml_array.inc'

end subroutine read_xml_double_array

subroutine read_xml_logical_array( info, tag, endtag, attribs, noattribs, data, &
                                   nodata, var, has_var )
   logical, dimension(:), pointer :: var

   include 'read_xml_array.inc'

end subroutine read_xml_logical_array

subroutine read_xml_word_array( info, tag, endtag, attribs, noattribs, data, &
                                nodata, var, has_var )
   character(len=*), dimension(:), pointer :: var

   include 'read_xml_array.inc'

end subroutine read_xml_word_array

! read_from_buffer_integers --
!    Routine to read all integers from a long string
!
! Arguments:
!    buffer      String containing the data
!    var         Variable to be filled
!    ierror      Error flag
!
subroutine read_from_buffer_integers( buffer, var, ierror )
   integer, dimension(:), pointer                :: var
   integer, dimension(:), pointer                :: work

   include 'read_from_buffer.inc'

end subroutine read_from_buffer_integers

! read_xml_from_buffer_reals, ... -
!    See read_xml_from_buffer_integers for an explanation
!
subroutine read_from_buffer_reals( buffer, var, ierror )
   real, dimension(:), pointer                :: var
   real, dimension(:), pointer                :: work

   include 'read_from_buffer.inc'

end subroutine read_from_buffer_reals

subroutine read_from_buffer_doubles( buffer, var, ierror )
   real(kind=kind(1.0d00)), dimension(:), pointer :: var
   real(kind=kind(1.0d00)), dimension(:), pointer :: work

   include 'read_from_buffer.inc'

end subroutine read_from_buffer_doubles

subroutine read_from_buffer_logicals( buffer, var, ierror )
   logical, dimension(:), pointer :: var
   logical, dimension(:), pointer :: work

   include 'read_from_buffer.inc'

end subroutine read_from_buffer_logicals

subroutine read_from_buffer_words( buffer, var, ierror )
   character(len=*), dimension(:), pointer :: var
   character(len=len(var)), dimension(:), pointer :: work

   include 'read_from_buffer.inc'

end subroutine read_from_buffer_words

! read_xml_word_1dim, ... -
!    Read an array of "words" (or ...) but from different elements
!
subroutine read_xml_integer_1dim( info, tag, endtag, attribs, noattribs, data, nodata, &
                                  var, has_var )
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   integer, dimension(:), pointer                :: var
   logical, intent(inout)                        :: has_var

   integer,dimension(:), pointer                 :: newvar
   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: newsize
   integer                                       :: ierr

   newsize = size(var) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = var
   deallocate( var )
   var => newvar

   call read_xml_integer( info, tag, endtag, attribs, noattribs, data, nodata, &
                          var(newsize), has_var )

end subroutine read_xml_integer_1dim

subroutine read_xml_real_1dim( info, tag, endtag, attribs, noattribs, data, nodata, &
                               var, has_var )
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   real, dimension(:), pointer                   :: var
   logical, intent(inout)                        :: has_var

   real, dimension(:), pointer                   :: newvar
   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: newsize
   integer                                       :: ierr

   newsize = size(var) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = var
   deallocate( var )
   var => newvar

   call read_xml_real( info, tag, endtag, attribs, noattribs, data, nodata, &
                       var(newsize), has_var )

end subroutine read_xml_real_1dim

subroutine read_xml_double_1dim( info, tag, endtag, attribs, noattribs, data, nodata, &
                                 var, has_var )
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   real(kind=kind(1.0d00)), dimension(:), pointer:: var
   logical, intent(inout)                        :: has_var

   real(kind=kind(1.0d00)), dimension(:), pointer:: newvar
   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: newsize
   integer                                       :: ierr

   newsize = size(var) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = var
   deallocate( var )
   var => newvar

   call read_xml_double( info, tag, endtag, attribs, noattribs, data, nodata, &
                       var(newsize), has_var )

end subroutine read_xml_double_1dim

subroutine read_xml_logical_1dim( info, tag, endtag, attribs, noattribs, data, nodata, &
                                  var, has_var )
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   logical, dimension(:), pointer                :: var
   logical, intent(inout)                        :: has_var

   logical, dimension(:), pointer                :: newvar
   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: newsize
   integer                                       :: ierr

   newsize = size(var) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = var
   deallocate( var )
   var => newvar

   call read_xml_logical( info, tag, endtag, attribs, noattribs, data, nodata, &
                          var(newsize), has_var )

end subroutine read_xml_logical_1dim

subroutine read_xml_word_1dim( info, tag, endtag, attribs, noattribs, data, nodata, &
                             var, has_var )
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   character(len=*), dimension(:), pointer       :: var
   logical, intent(inout)                        :: has_var

   character(len=len(var)),dimension(:), pointer :: newvar
   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: newsize
   integer                                       :: ierr

   newsize = size(var) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = var
   deallocate( var )
   var => newvar

   call read_xml_word( info, tag, endtag, attribs, noattribs, data, nodata, &
                       var(newsize), has_var )

end subroutine read_xml_word_1dim

subroutine read_xml_line_1dim( info, tag, endtag, attribs, noattribs, data, nodata, &
                             var, has_var )
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   character(len=*), dimension(:), pointer       :: var
   logical, intent(inout)                        :: has_var

   character(len=len(var)),dimension(:), pointer :: newvar
   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: newsize
   integer                                       :: ierr

   newsize = size(var) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = var
   deallocate( var )
   var => newvar

   call read_xml_line( info, tag, endtag, attribs, noattribs, data, nodata, &
                       var(newsize), has_var )

end subroutine read_xml_line_1dim


end module read_xml_primitives
