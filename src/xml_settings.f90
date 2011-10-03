module xml_data_settings_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   integer, private :: lurep_
   logical, private :: strict_

type xslibrary_xml
   character(len=250)                                :: path
end type xslibrary_xml

type criticality_xml
   integer                                         :: cycles
   integer                                         :: inactive
   integer                                         :: particles
end type criticality_xml

type source_xml
   character(len=10)                                :: type
   real(kind=kind(1.0d0)), dimension(:), pointer   :: coeffs => null()
end type source_xml

type cutoff_xml
   real(kind=kind(1.0d0))                          :: weight
   real(kind=kind(1.0d0))                          :: weight_avg
end type cutoff_xml
   type(xslibrary_xml)                             :: xslibrary
   type(criticality_xml)                           :: criticality
   integer                                         :: verbosity_
   type(source_xml)                                :: source_
   character(len=3), dimension(:), pointer         :: survival_ => null()
   type(cutoff_xml)                                :: cutoff_
contains
subroutine read_xml_type_xslibrary_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(xslibrary_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(xslibrary_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_xslibrary_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_xslibrary_xml_array

subroutine read_xml_type_xslibrary_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(xslibrary_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_path
   has_path                             = .false.
   call init_xml_type_xslibrary_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('path')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%path, has_path )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_path ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on path')
   endif
end subroutine read_xml_type_xslibrary_xml
subroutine init_xml_type_xslibrary_xml_array( dvar )
   type(xslibrary_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_xslibrary_xml_array
subroutine init_xml_type_xslibrary_xml(dvar)
   type(xslibrary_xml) :: dvar
end subroutine init_xml_type_xslibrary_xml
subroutine write_xml_type_xslibrary_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(xslibrary_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_xslibrary_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_xslibrary_xml_array

subroutine write_xml_type_xslibrary_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(xslibrary_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_word( info, 'path', indent+3, dvar%path)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_xslibrary_xml

subroutine read_xml_type_criticality_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(criticality_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(criticality_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_criticality_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_criticality_xml_array

subroutine read_xml_type_criticality_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(criticality_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_cycles
   logical                                         :: has_inactive
   logical                                         :: has_particles
   has_cycles                           = .false.
   has_inactive                         = .false.
   has_particles                        = .false.
   call init_xml_type_criticality_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('cycles')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%cycles, has_cycles )
      case('inactive')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%inactive, has_inactive )
      case('particles')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%particles, has_particles )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_cycles ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on cycles')
   endif
   if ( .not. has_inactive ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on inactive')
   endif
   if ( .not. has_particles ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on particles')
   endif
end subroutine read_xml_type_criticality_xml
subroutine init_xml_type_criticality_xml_array( dvar )
   type(criticality_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_criticality_xml_array
subroutine init_xml_type_criticality_xml(dvar)
   type(criticality_xml) :: dvar
end subroutine init_xml_type_criticality_xml
subroutine write_xml_type_criticality_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(criticality_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_criticality_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_criticality_xml_array

subroutine write_xml_type_criticality_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(criticality_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'cycles', indent+3, dvar%cycles)
   call write_to_xml_integer( info, 'inactive', indent+3, dvar%inactive)
   call write_to_xml_integer( info, 'particles', indent+3, dvar%particles)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_criticality_xml

subroutine read_xml_type_source_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(source_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(source_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_source_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_source_xml_array

subroutine read_xml_type_source_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(source_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_type
   logical                                         :: has_coeffs
   has_type                             = .false.
   has_coeffs                           = .false.
   call init_xml_type_source_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
      case('coeffs')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%coeffs, has_coeffs )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on type')
   endif
   if ( .not. has_coeffs ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on coeffs')
   endif
end subroutine read_xml_type_source_xml
subroutine init_xml_type_source_xml_array( dvar )
   type(source_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_source_xml_array
subroutine init_xml_type_source_xml(dvar)
   type(source_xml) :: dvar
end subroutine init_xml_type_source_xml
subroutine write_xml_type_source_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(source_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_source_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_source_xml_array

subroutine write_xml_type_source_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(source_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
   call write_to_xml_double_array( info, 'coeffs', indent+3, dvar%coeffs)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_source_xml

subroutine read_xml_type_cutoff_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(cutoff_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(cutoff_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_cutoff_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_cutoff_xml_array

subroutine read_xml_type_cutoff_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(cutoff_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_weight
   logical                                         :: has_weight_avg
   has_weight                           = .false.
   has_weight_avg                       = .false.
   call init_xml_type_cutoff_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('weight')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%weight, has_weight )
      case('weight_avg')
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%weight_avg, has_weight_avg )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
end subroutine read_xml_type_cutoff_xml
subroutine init_xml_type_cutoff_xml_array( dvar )
   type(cutoff_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_cutoff_xml_array
subroutine init_xml_type_cutoff_xml(dvar)
   type(cutoff_xml) :: dvar
   dvar%weight = 0.25
   dvar%weight_avg = 1.0
end subroutine init_xml_type_cutoff_xml
subroutine write_xml_type_cutoff_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(cutoff_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_cutoff_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_cutoff_xml_array

subroutine write_xml_type_cutoff_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(cutoff_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_double( info, 'weight', indent+3, dvar%weight)
   call write_to_xml_double( info, 'weight_avg', indent+3, dvar%weight_avg)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_cutoff_xml

subroutine read_xml_file_settings_t(fname, lurep, errout)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep
   logical, intent(out), optional         :: errout

   type(XML_PARSE)                        :: info
   logical                                :: error
   character(len=80)                      :: tag
   character(len=80)                      :: starttag
   logical                                :: endtag
   character(len=80), dimension(1:2,1:20) :: attribs
   integer                                :: noattribs
   character(len=200), dimension(1:100)   :: data
   integer                                :: nodata
   logical                                         :: has_xslibrary
   logical                                         :: has_criticality
   logical                                         :: has_verbosity_
   logical                                         :: has_source_
   logical                                         :: has_survival_
   logical                                         :: has_cutoff_
   has_xslibrary                        = .false.
   has_criticality                      = .false.
   has_verbosity_                       = .false.
   has_source_                          = .false.
   has_survival_                        = .false.
   allocate(survival_(0))
   has_cutoff_                          = .false.

   call init_xml_file_settings_t
   call xml_open( info, fname, .true. )
   call xml_options( info, report_errors=.false., ignore_whitespace=.true.)
   lurep_ = 0
   if ( present(lurep) ) then
      lurep_ = lurep
      call xml_options( info, report_lun=lurep )
   endif
   do
      call xml_get( info, starttag, endtag, attribs, noattribs, &
         data, nodata)
      if ( starttag /= '!--' ) exit
   enddo
   if ( starttag /= "settings" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "settings"')
      error = .true.
      call xml_close(info)
      return
   endif
   strict_ = .false.
   error = .false.
   do
      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
      if ( xml_error(info) ) then
         write(lurep_,*) 'Error reading input file!'
         error = .true.
         return
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('xslibrary')
         call read_xml_type_xslibrary_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            xslibrary, has_xslibrary )
      case('criticality')
         call read_xml_type_criticality_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            criticality, has_criticality )
      case('verbosity')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            verbosity_, has_verbosity_ )
      case('source')
         call read_xml_type_source_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            source_, has_source_ )
      case('survival_biasing')
         call read_xml_word_1dim( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            survival_, has_survival_ )
      case('cutoff')
         call read_xml_type_cutoff_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            cutoff_, has_cutoff_ )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_xslibrary ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on xslibrary')
   endif
   if ( .not. has_criticality ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on criticality')
   endif
   if ( .not. has_verbosity_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on verbosity_')
   endif
   if ( .not. has_source_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on source_')
   endif
   if ( .not. has_survival_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on survival_')
   endif
   if ( .not. has_cutoff_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on cutoff_')
   endif
   if ( present(errout) ) errout = error
end subroutine

subroutine write_xml_file_settings_t(fname, lurep)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep

   type(XML_PARSE)                        :: info
   integer                                :: indent = 0

   call xml_open( info, fname, .false. )
   call xml_options( info, report_errors=.true.)
   if ( present(lurep) ) then
       call xml_options( info, report_errors=.true.)
   endif
   write(info%lun,'(a)') &
      '<settings>'
   call write_xml_type_xslibrary_xml( info, 'xslibrary', indent+3, xslibrary)
   call write_xml_type_criticality_xml( info, 'criticality', indent+3, criticality)
   call write_to_xml_integer( info, 'verbosity', indent+3, verbosity_)
   call write_xml_type_source_xml( info, 'source', indent+3, source_)
   call write_to_xml_word_1dim( info, 'survival_biasing', indent+3, survival_)
   call write_xml_type_cutoff_xml( info, 'cutoff', indent+3, cutoff_)
   write(info%lun,'(a)') '</settings>'
   call xml_close(info)
end subroutine

subroutine init_xml_file_settings_t

end subroutine

end module
