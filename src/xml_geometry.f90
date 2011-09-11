module xml_data_geometry_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   integer, private :: lurep_
   logical, private :: strict_

type cell_xml
   integer                                         :: uid
   integer                                         :: universe
   integer                                         :: material
   integer                                         :: fill
   integer, dimension(:), pointer                  :: surfaces => null()
end type cell_xml

type surface_xml
   integer                                         :: uid
   character(len=15)                                :: type
   real(kind=kind(1.0d0)), dimension(:), pointer   :: coeffs => null()
   character(len=12)                                :: boundary
end type surface_xml

type lattice_xml
   integer                                         :: uid
   character(len=12)                                :: type
   integer, dimension(:), pointer                  :: dimension => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: origin => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: width => null()
   integer, dimension(:), pointer                  :: universes => null()
end type lattice_xml
   type(cell_xml), dimension(:), pointer           :: cell_ => null()
   type(surface_xml), dimension(:), pointer        :: surface_ => null()
   type(lattice_xml), dimension(:), pointer        :: lattice_ => null()
contains
subroutine read_xml_type_cell_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(cell_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(cell_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_cell_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_cell_xml_array

subroutine read_xml_type_cell_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(cell_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_uid
   logical                                         :: has_universe
   logical                                         :: has_material
   logical                                         :: has_fill
   logical                                         :: has_surfaces
   has_uid                              = .false.
   has_universe                         = .false.
   has_material                         = .false.
   has_fill                             = .false.
   has_surfaces                         = .false.
   call init_xml_type_cell_xml(dvar)
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
      case('uid')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%uid, has_uid )
      case('universe')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%universe, has_universe )
      case('material')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%material, has_material )
      case('fill')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%fill, has_fill )
      case('surfaces')
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%surfaces, has_surfaces )
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
   if ( .not. has_uid ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on uid')
   endif
   if ( .not. has_surfaces ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on surfaces')
   endif
end subroutine read_xml_type_cell_xml
subroutine init_xml_type_cell_xml_array( dvar )
   type(cell_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_cell_xml_array
subroutine init_xml_type_cell_xml(dvar)
   type(cell_xml) :: dvar
   dvar%universe = 0
   dvar%material = 0
   dvar%fill = 0
end subroutine init_xml_type_cell_xml
subroutine write_xml_type_cell_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(cell_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_cell_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_cell_xml_array

subroutine write_xml_type_cell_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(cell_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'uid', indent+3, dvar%uid)
   call write_to_xml_integer( info, 'universe', indent+3, dvar%universe)
   call write_to_xml_integer( info, 'material', indent+3, dvar%material)
   call write_to_xml_integer( info, 'fill', indent+3, dvar%fill)
   call write_to_xml_integer_array( info, 'surfaces', indent+3, dvar%surfaces)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_cell_xml

subroutine read_xml_type_surface_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(surface_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(surface_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_surface_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_surface_xml_array

subroutine read_xml_type_surface_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(surface_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_uid
   logical                                         :: has_type
   logical                                         :: has_coeffs
   logical                                         :: has_boundary
   has_uid                              = .false.
   has_type                             = .false.
   has_coeffs                           = .false.
   has_boundary                         = .false.
   call init_xml_type_surface_xml(dvar)
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
      case('uid')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%uid, has_uid )
      case('type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
      case('coeffs')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%coeffs, has_coeffs )
      case('boundary')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%boundary, has_boundary )
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
   if ( .not. has_uid ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on uid')
   endif
   if ( .not. has_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on type')
   endif
   if ( .not. has_coeffs ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on coeffs')
   endif
end subroutine read_xml_type_surface_xml
subroutine init_xml_type_surface_xml_array( dvar )
   type(surface_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_surface_xml_array
subroutine init_xml_type_surface_xml(dvar)
   type(surface_xml) :: dvar
   dvar%boundary = 'transmit'
end subroutine init_xml_type_surface_xml
subroutine write_xml_type_surface_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(surface_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_surface_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_surface_xml_array

subroutine write_xml_type_surface_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(surface_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'uid', indent+3, dvar%uid)
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
   call write_to_xml_double_array( info, 'coeffs', indent+3, dvar%coeffs)
   call write_to_xml_word( info, 'boundary', indent+3, dvar%boundary)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_surface_xml

subroutine read_xml_type_lattice_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(lattice_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(lattice_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_lattice_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_lattice_xml_array

subroutine read_xml_type_lattice_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(lattice_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_uid
   logical                                         :: has_type
   logical                                         :: has_dimension
   logical                                         :: has_origin
   logical                                         :: has_width
   logical                                         :: has_universes
   has_uid                              = .false.
   has_type                             = .false.
   has_dimension                        = .false.
   has_origin                           = .false.
   has_width                            = .false.
   has_universes                        = .false.
   call init_xml_type_lattice_xml(dvar)
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
      case('uid')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%uid, has_uid )
      case('type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
      case('dimension')
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%dimension, has_dimension )
      case('origin')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%origin, has_origin )
      case('width')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%width, has_width )
      case('universes')
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%universes, has_universes )
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
   if ( .not. has_uid ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on uid')
   endif
   if ( .not. has_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on type')
   endif
   if ( .not. has_dimension ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on dimension')
   endif
   if ( .not. has_origin ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on origin')
   endif
   if ( .not. has_width ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on width')
   endif
   if ( .not. has_universes ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on universes')
   endif
end subroutine read_xml_type_lattice_xml
subroutine init_xml_type_lattice_xml_array( dvar )
   type(lattice_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_lattice_xml_array
subroutine init_xml_type_lattice_xml(dvar)
   type(lattice_xml) :: dvar
end subroutine init_xml_type_lattice_xml
subroutine write_xml_type_lattice_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(lattice_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_lattice_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_lattice_xml_array

subroutine write_xml_type_lattice_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(lattice_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'uid', indent+3, dvar%uid)
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
   call write_to_xml_integer_array( info, 'dimension', indent+3, dvar%dimension)
   call write_to_xml_double_array( info, 'origin', indent+3, dvar%origin)
   call write_to_xml_double_array( info, 'width', indent+3, dvar%width)
   call write_to_xml_integer_array( info, 'universes', indent+3, dvar%universes)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_lattice_xml

subroutine read_xml_file_geometry_t(fname, lurep, errout)
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
   logical                                         :: has_cell_
   logical                                         :: has_surface_
   logical                                         :: has_lattice_
   has_cell_                            = .false.
   allocate(cell_(0))
   has_surface_                         = .false.
   allocate(surface_(0))
   has_lattice_                         = .false.
   allocate(lattice_(0))

   call init_xml_file_geometry_t
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
   if ( starttag /= "geometry" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "geometry"')
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
      case('cell')
         call read_xml_type_cell_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            cell_, has_cell_ )
      case('surface')
         call read_xml_type_surface_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            surface_, has_surface_ )
      case('lattice')
         call read_xml_type_lattice_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            lattice_, has_lattice_ )
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
   if ( .not. has_cell_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on cell_')
   endif
   if ( .not. has_surface_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on surface_')
   endif
   if ( .not. has_lattice_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on lattice_')
   endif
   if ( present(errout) ) errout = error
end subroutine

subroutine write_xml_file_geometry_t(fname, lurep)
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
      '<geometry>'
   call write_xml_type_cell_xml_array( info, 'cell', indent+3, cell_)
   call write_xml_type_surface_xml_array( info, 'surface', indent+3, surface_)
   call write_xml_type_lattice_xml_array( info, 'lattice', indent+3, lattice_)
   write(info%lun,'(a)') '</geometry>'
   call xml_close(info)
end subroutine

subroutine init_xml_file_geometry_t

end subroutine

end module
