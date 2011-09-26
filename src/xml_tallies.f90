module xml_data_tallies_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   integer, private :: lurep_
   logical, private :: strict_

type filter_xml
   character(len=250)                                :: cell
   character(len=250)                                :: surface
   character(len=250)                                :: universe
   character(len=250)                                :: material
   character(len=250)                                :: mesh
   character(len=250)                                :: cellborn
   character(len=250)                                :: energy
   character(len=250)                                :: energyout
end type filter_xml

type tally_xml
   integer                                         :: id
   type(filter_xml)                                :: filters
   character(len=250)                                :: macros
   character(len=250)                                :: reactions
   character(len=250)                                :: nuclides
end type tally_xml
   type(tally_xml), dimension(:), pointer          :: tally_ => null()
contains
subroutine read_xml_type_filter_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(filter_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(filter_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_filter_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_filter_xml_array

subroutine read_xml_type_filter_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(filter_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_cell
   logical                                         :: has_surface
   logical                                         :: has_universe
   logical                                         :: has_material
   logical                                         :: has_mesh
   logical                                         :: has_cellborn
   logical                                         :: has_energy
   logical                                         :: has_energyout
   has_cell                             = .false.
   has_surface                          = .false.
   has_universe                         = .false.
   has_material                         = .false.
   has_mesh                             = .false.
   has_cellborn                         = .false.
   has_energy                           = .false.
   has_energyout                        = .false.
   call init_xml_type_filter_xml(dvar)
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
      case('cell')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%cell, has_cell )
      case('surface')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%surface, has_surface )
      case('universe')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%universe, has_universe )
      case('material')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%material, has_material )
      case('mesh')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%mesh, has_mesh )
      case('cellborn')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%cellborn, has_cellborn )
      case('energy')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%energy, has_energy )
      case('energyout')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%energyout, has_energyout )
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
end subroutine read_xml_type_filter_xml
subroutine init_xml_type_filter_xml_array( dvar )
   type(filter_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_filter_xml_array
subroutine init_xml_type_filter_xml(dvar)
   type(filter_xml) :: dvar
   dvar%cell = ''
   dvar%surface = ''
   dvar%universe = ''
   dvar%material = ''
   dvar%mesh = ''
   dvar%cellborn = ''
   dvar%energy = ''
   dvar%energyout = ''
end subroutine init_xml_type_filter_xml
subroutine write_xml_type_filter_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(filter_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_filter_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_filter_xml_array

subroutine write_xml_type_filter_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(filter_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_word( info, 'cell', indent+3, dvar%cell)
   call write_to_xml_word( info, 'surface', indent+3, dvar%surface)
   call write_to_xml_word( info, 'universe', indent+3, dvar%universe)
   call write_to_xml_word( info, 'material', indent+3, dvar%material)
   call write_to_xml_word( info, 'mesh', indent+3, dvar%mesh)
   call write_to_xml_word( info, 'cellborn', indent+3, dvar%cellborn)
   call write_to_xml_word( info, 'energy', indent+3, dvar%energy)
   call write_to_xml_word( info, 'energyout', indent+3, dvar%energyout)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_filter_xml

subroutine read_xml_type_tally_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(tally_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_tally_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_tally_xml_array

subroutine read_xml_type_tally_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_id
   logical                                         :: has_filters
   logical                                         :: has_macros
   logical                                         :: has_reactions
   logical                                         :: has_nuclides
   has_id                               = .false.
   has_filters                          = .false.
   has_macros                           = .false.
   has_reactions                        = .false.
   has_nuclides                         = .false.
   call init_xml_type_tally_xml(dvar)
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
      case('id')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%id, has_id )
      case('filters')
         call read_xml_type_filter_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%filters, has_filters )
      case('macros')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%macros, has_macros )
      case('reactions')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%reactions, has_reactions )
      case('nuclides')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nuclides, has_nuclides )
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
   if ( .not. has_id ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on id')
   endif
   if ( .not. has_filters ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on filters')
   endif
end subroutine read_xml_type_tally_xml
subroutine init_xml_type_tally_xml_array( dvar )
   type(tally_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_tally_xml_array
subroutine init_xml_type_tally_xml(dvar)
   type(tally_xml) :: dvar
   dvar%macros = ''
   dvar%reactions = ''
   dvar%nuclides = ''
end subroutine init_xml_type_tally_xml
subroutine write_xml_type_tally_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_tally_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_tally_xml_array

subroutine write_xml_type_tally_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
   call write_xml_type_filter_xml( info, 'filters', indent+3, dvar%filters)
   call write_to_xml_word( info, 'macros', indent+3, dvar%macros)
   call write_to_xml_word( info, 'reactions', indent+3, dvar%reactions)
   call write_to_xml_word( info, 'nuclides', indent+3, dvar%nuclides)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_tally_xml

subroutine read_xml_file_tallies_t(fname, lurep, errout)
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
   logical                                         :: has_tally_
   has_tally_                           = .false.
   allocate(tally_(0))

   call init_xml_file_tallies_t
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
   if ( starttag /= "tallies" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "tallies"')
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
      case('tally')
         call read_xml_type_tally_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            tally_, has_tally_ )
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
   if ( .not. has_tally_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on tally_')
   endif
   if ( present(errout) ) errout = error
end subroutine

subroutine write_xml_file_tallies_t(fname, lurep)
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
      '<tallies>'
   call write_xml_type_tally_xml_array( info, 'tally', indent+3, tally_)
   write(info%lun,'(a)') '</tallies>'
   call xml_close(info)
end subroutine

subroutine init_xml_file_tallies_t

end subroutine

end module
