! Part of XML-Fortran library:
!
! $Id: read_xml_scalar.inc,v 1.3 2007/02/26 20:33:38 arjenmarkus Exp $
!
   type(XML_PARSE), intent(inout)               :: info
   character(len=*), intent(in)                 :: tag
   logical, intent(inout)                       :: endtag
   character(len=*), dimension(:,:), intent(in) :: attribs
   integer, intent(in)                          :: noattribs
   character(len=*), dimension(:), intent(in)   :: data
   integer, intent(in)                          :: nodata
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
   if ( idx .gt. 0 ) then
      read( buffer, *, iostat=ierr ) var
      has_var = .true.
   else
      do idx = 1,nodata
         if ( data(idx) .ne. ' ' ) then
            read( data(idx), *, iostat=ierr ) var
            has_var = .true.
            exit
         endif
      enddo
   endif

   if ( ierr .ne. 0 ) then
      write(*,*) 'Error reading variable - tag = ', trim(tag)
      has_var = .false.
   endif
