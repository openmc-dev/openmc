!===============================================================================
! XMLREADER:
!    Read an XML-file that contains a template for the data and the XML files
!    that will be used in a program and generate a module with reading routines.
!
! TODO:
! - Private routines!
! - Error for unknown data types
! - Length for character items
!===============================================================================

program xmlreader
   use XMLPARSE
   implicit none

   character(len=60)                      :: fname
   type(XML_PARSE)                        :: info

   character(len=80)                      :: tag
   character(len=80)                      :: starttag
   logical                                :: endtag
   character(len=80), dimension(1:2,1:20) :: attribs
   integer                                :: noattribs
   character(len=100000), dimension(1:1000)  :: data
   integer                                :: nodata
   integer                                :: i
   integer                                :: j
   integer, parameter                     :: notmps = 6 ! Number of temporary files needed!
   integer                                :: ludef
   integer                                :: ludeflt
   integer                                :: luinit
   integer                                :: lusubs
   integer                                :: luprolog
   integer                                :: lustart
   integer                                :: luloop
   integer                                :: luend
   integer                                :: lumain
   integer                                :: luwrite
   integer                                :: luwritp
   integer                                :: luwrv
   integer                                :: argc
   logical                                :: prolog_written
   logical                                :: comp
   logical                                :: error
   logical                                :: begin_loop = .true.
   logical                                :: begin_main_loop = .true.
   logical                                :: begin_component = .false.
   logical                                :: strict
   logical                                :: global_type
   logical                                :: dyn_strings

   character(len=32)                      :: root_name
   character(len=32)                      :: global_name
   character(len=32)                      :: typename
   character(len=32), dimension(1:20)     :: placeholder
   integer                                :: no_placeholders

   character(len=50)                      :: declare
   character(len=50), dimension(:,:), pointer :: types
   character(len=50), dimension(:,:), pointer :: new_types
   
   ! TODO: arrays of integers etc, in addition to integer-arrays - see: word
   integer, parameter                     :: notypes_predefined = 27
   character(len=50), dimension(1:4,1:notypes_predefined) :: predefined_types
   integer                                :: notypes = notypes_predefined
   data ((predefined_types(i,j) , i=1,4), j=1,notypes_predefined ) / &
'logical'       ,'   logical'                        , 'read_xml_logical',       'write_to_xml_logical', &
'logical-1dim'  ,'   logical, dimension(:), pointer' , 'read_xml_logical_1dim',  'write_to_xml_logical_1dim',  &
'logical-array' ,'   logical, dimension(:), pointer' , 'read_xml_logical_array', 'write_to_xml_logical_array', &
'logical-shape' ,'   logical, dimension(SHAPE)'      , 'read_xml_logical_array', 'write_to_xml_logical_array', &
'integer'       ,'   integer'                        , 'read_xml_integer',       'write_to_xml_integer', &
'integer-1dim'  ,'   integer, dimension(:), pointer' , 'read_xml_integer_1dim',  'write_to_xml_integer_1dim',  &
'integer-array' ,'   integer, dimension(:), pointer' , 'read_xml_integer_array', 'write_to_xml_integer_array', &
'integer-shape' ,'   integer, dimension(SHAPE)'      , 'read_xml_integer_array', 'write_to_xml_integer_array', &
'real'          ,'   real'                           , 'read_xml_real'   ,       'write_to_xml_real'   , &
'real-1dim'     ,'   real, dimension(:), pointer'    , 'read_xml_real_1dim',     'write_to_xml_real_1dim', &
'real-array'    ,'   real, dimension(:), pointer'    , 'read_xml_real_array',    'write_to_xml_real_array', &
'real-shape'    ,'   real, dimension(SHAPE)'         , 'read_xml_real_array',    'write_to_xml_real_array', &
'double'        ,'   real(kind=kind(1.0d0))'         , 'read_xml_double',        'write_to_xml_double',  &
'double-1dim'   ,'   real(kind=kind(1.0d0)), dimension(:), pointer' , &
                                                       'read_xml_double_1dim',   'write_to_xml_double_1dim', &
'double-array'  ,'   real(kind=kind(1.0d0)), dimension(:), pointer' , &
                                                       'read_xml_double_array',  'write_to_xml_double_array', &
'double-shape'  ,'   real(kind=kind(1.0d0)), dimension(SHAPE)'      , &
                                                       'read_xml_double_array',  'write_to_xml_double_array', &
'word'          ,'   character(len=?)'               , 'read_xml_word',          'write_to_xml_word',    &
'word-1dim'     ,'   character(len=?), dimension(:), pointer' , 'read_xml_word_1dim',  'write_to_xml_word_1dim',  &
'word-array'    ,'   character(len=?), dimension(:), pointer' , 'read_xml_word_array', 'write_to_xml_word_array', &
'word-shape'    ,'   character(len=?), dimension(SHAPE)'      , 'read_xml_word_array', 'write_to_xml_word_array', &
'line'          ,'   character(len=?)'               , 'read_xml_line',          'write_to_xml_line',    &
'line-1dim'     ,'   character(len=?), dimension(:), pointer'  , 'read_xml_line_1dim',  'write_to_xml_line_1dim',  &
'line-array'    ,'   character(len=?), dimension(:), pointer'  , 'read_xml_line_array', 'write_to_xml_line_array', &
'line-shape'    ,'   character(len=?), dimension(SHAPE)'       , 'read_xml_line_array', 'write_to_xml_line_array', &
'character'     ,'   character(len=?)'               , 'read_xml_line',                 'write_to_xml_line',    &
'character-array','   character(len=?), dimension(:), pointer' , 'read_xml_line_array', 'write_to_xml_line_array', &
'character-shape','   character(len=?), dimension(SHAPE)'      , 'read_xml_line_array', 'write_to_xml_line_array' /

   allocate( types(1:4,1:notypes) )
   types = predefined_types(:,1:notypes)

   !
   ! Read the global options file, if present
   !
   strict      = .false.
   global_type = .false.
   dyn_strings = .true.
   call get_global_options( attribs, noattribs, strict, global_type, global_name, &
                            root_name, dyn_strings )

   !
   ! Open the input file and read the name of the template.
   ! Load the template into a tree and then generate it all
   ! in stages
   !
   argc = COMMAND_ARGUMENT_COUNT()
   if (argc > 0) then
      call GET_COMMAND_ARGUMENT(1, fname)
   else
      open(UNIT=10, FILE='xmlreader.inp')
      read(UNIT=10, FMT='(a)') fname
      close(UNIT=10)
   end if

   open( 20, file = 'xmlreader.out' )
   call xml_options( info, report_lun = 20, report_details = .true. )

   prolog_written = .false.
   !
   ! Set the defaults
   !
   global_type = .false.
   global_name = fname
   root_name   = fname

   ludef    = 21
   lusubs   = 22
   luinit   = 23
   luwritp  = 24
   luwrv    = luwritp
   open( ludef,   file = trim(fname)//'.f90' )
   open( lusubs,  status = 'scratch' )
   open( luinit,  status = 'scratch' )
   open( luwritp, status = 'scratch' )
   call open_tmp_files( 31 )

   lumain   = luloop

   ! Read the template file and act as we go along:
   ! - write the declarations
   ! - write the main reading routine
   ! - write the individual reading routines
   ! - the root element is needed because of the definition of XML files,
   !   but we ignore it.
   !
   call xml_open( info, trim(fname)//'.xml', .true. )

   error = .false.
   comp  = .false.
   no_placeholders = 0

   call xml_get( info, starttag, endtag, attribs, noattribs, data, nodata )

   do
      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
      write(20,*) 'tag: ',tag
      write(20,*) 'attribs: ',noattribs
      if ( noattribs > 0 ) then
         write(20,'(4a)') ( '   ', trim(attribs(1,i)), ' = ', trim(attribs(2,i)), i=1,noattribs )
      endif
      write(20,*) 'data: ',nodata
      if ( nodata > 0 ) then
         write(20,'(3a)') ( '   >', trim(data(i)), '<', i=1,nodata )
      endif
      if ( xml_error(info) ) then
         write(*,*) 'Error reading template file!'
         !stop
         exit
      endif

      !
      ! When encountering the endtag, then close the
      ! current definition
      !
      if ( endtag .and. noattribs == 0 ) then
         select case ( tag )
         case ( 'typedef' )
            call close_typedef( begin_component )
         case ( 'placeholder' )
           !if ( comp ) then
           !   call close_placeholder
           !endif
           !comp = .false.
            call close_placeholder
            luwrv = luwritp
         case default
            !
            ! Have we found the end of the definition?
            !
            if ( tag == starttag ) then
               exit
            endif
         end select

         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif

      endif

      !
      ! Opening tags: dispatch on the actual tag
      !
      select case( tag )
      case( 'options' )
         call set_options( attribs, noattribs, strict, global_type, &
                           global_name, root_name, dyn_strings )
         if ( .not. prolog_written ) then
             prolog_written = .true.
             call write_prolog
         else
             write(20,*) 'Options element should be the first child of the root element',&
                         'Otherwise it has no effect!'
         endif

         if ( strict ) then
            write( lumain, '(a)' ) '   strict_ = .true.'
         else
            write( lumain, '(a)' ) '   strict_ = .false.'
         endif

      case( 'comment', '!--' )
         ! Do nothing

      case( 'placeholder' )
         if ( .not. prolog_written ) then
             prolog_written = .true.
             call write_prolog
         endif
         if ( begin_loop .or. begin_main_loop ) then
            begin_main_loop = .false.
            call add_begin_loop( .true., .false. )
         endif
         call add_placeholder(dyn_strings)
         begin_component = .false.
         luwrv = luwrite

      case( 'typedef' )
         if ( .not. prolog_written ) then
             prolog_written = .true.
             call write_prolog
         endif
         call add_typedef(dyn_strings)

      case( 'variable' )
         if ( .not. prolog_written ) then
             prolog_written = .true.
             call write_prolog
         endif
         if ( begin_loop .or. begin_main_loop ) then
            begin_main_loop = .false.
            call add_begin_loop( .true., begin_component )
         endif
         call add_variable( component=comp )

      case( 'component' )
         !
         ! Components of derived types are treated in much the
         ! same way as ordinary variables - with one syntactic
         ! difference
         !
         if ( .not. prolog_written ) then
             prolog_written = .true.
             call write_prolog
         endif
         if ( begin_loop ) then
            call add_begin_loop( .true., .true. )
         endif
         call add_variable( component=.true. )
         begin_component = .true.

      case default
         write(20,*) 'Unknown tag: ',trim(tag)
         write(20,*) '-- terminating the program!'
         stop
      end select
   end do

   !
   ! Now finish it all
   !
   write( luend, '(a)' ) &
      & '   if ( present(errout) ) errout = error', &
      & '   call xml_close(info)', &
      & 'end subroutine', &
      & ' '
   write( luwritp, '(a)' ) &
      & '   write(info%lun,''(a)'') ''</' // trim(root_name) // '>''', &
      & '   call xml_close(info)', &
      & 'end subroutine', &
      & ' '

   call append_files( luprolog )
   call merge_files

   write( ludef, '(/,a)' ) &
      & 'end subroutine', &
      & 'end module'

   if ( error ) then
       write(*,*) 'Errors found in the definition - please check!'
   endif
   stop
contains

!===============================================================================
! GET_GLOBAL_OPTIONS --
!    Routine to get the global options from the configuration file
! Arguments:
!    strict            Option to make the parser check for unknown tags
!    global_type       Option to generate an overall derived type
!    global_name       Name of that overall derived type (if requested)
!    root_name         Name of the root element of the XML file
!    dyn_strings       Whether to use dynamic strings or not
!===============================================================================

subroutine get_global_options( attribs, noattribs, strict, global_type, global_name, &
                               root_name, dyn_strings )
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   logical, intent(inout)                          :: strict
   logical, intent(inout)                          :: global_type
   character(len=*), intent(inout)                 :: global_name
   character(len=*), intent(inout)                 :: root_name
   logical, intent(inout)                          :: dyn_strings

   character(len=20)                               :: tag
   character(len=20), dimension(1)                 :: data
   integer                                         :: nodata
   logical                                         :: exists
   logical                                         :: endtag
   type(XML_PARSE)                                 :: info

   inquire( file = 'xmlreader.conf', exist = exists )
   if ( exists ) then
       call xml_open( info, 'xmlreader.conf', .true. )
       do while ( xml_ok(info) )
          call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
          if ( tag == 'xmlreader' ) then
              call set_options( attribs, noattribs, strict, global_type, &
                                global_name, root_name, dyn_strings )
          endif
       enddo
       call xml_close( info )
   endif
end subroutine get_global_options

!===============================================================================
! SET_OPTIONS --
!    Routine to set the options that influence the parser
! Arguments:
!    attribs           List of attributes
!    noattribs         Number of attributes
!    strict            Option to make the parser check for unknown tags
!    global_type       Option to generate an overall derived type
!    global_name       Name of that overall derived type (if requested)
!    root_name         Name of the root element of the XML file
!    dyn_strings       Whether to use dynamic strings or not
!===============================================================================

subroutine set_options( attribs, noattribs, strict, global_type, global_name, root_name, dyn_strings )
   character(len=*), dimension(:,:), intent(in) :: attribs
   integer, intent(in)                          :: noattribs
   logical, intent(inout)                       :: strict
   logical, intent(inout)                       :: global_type
   character(len=*), intent(inout)              :: global_name
   character(len=*), intent(inout)              :: root_name
   logical, intent(inout)                       :: dyn_strings

   integer :: i

   do i = 1,noattribs
      select case (attribs(1,i))
      case ('strict')
         if ( attribs(2,i) == 'yes' ) then
            strict = .true.
         else
            strict = .false.
         endif
      case ('globaltype')
         if ( attribs(2,i) == 'yes' ) then
            global_type = .true.
         else
            global_type = .false.
         endif
      case ('globalname')
         global_name = attribs(2,i)
      case ('rootname')
         root_name = attribs(2,i)
      case ('dynamicstrings')
         if ( attribs(2,i) == 'yes' ) then
            dyn_strings = .true.
         else
            dyn_strings = .false.
         endif
      case default
         write(20,*) 'Unknown option: ',trim(attribs(1,i)), ' - ignored'
      end select
   enddo
end subroutine set_options

!===============================================================================
! OPEN_TMP_FILES --
!    Routine to open the temporary files
! Arguments:
!    lufirst           First LU-number to use
!===============================================================================

subroutine open_tmp_files( lufirst )
   integer, intent(in) :: lufirst

   luprolog = lufirst
   lustart  = lufirst + 1
   luloop   = lufirst + 2
   luend    = lufirst + 3
   ludeflt  = lufirst + 4
   luwrite  = lufirst + 5
   open( luprolog, status = 'scratch' )
   open( lustart,  status = 'scratch' )
   open( luloop,   status = 'scratch' )
   open( luend,    status = 'scratch' )
   open( ludeflt,  status = 'scratch' )
   open( luwrite,  status = 'scratch' )
end subroutine open_tmp_files

!===============================================================================
! CLOSE_TMP_FILES --
!    Routine to close the temporary files
! Arguments:
!    None
!===============================================================================

subroutine close_tmp_files
   close( luprolog )
   close( lustart  )
   close( luloop   )
   close( luend    )
   close( ludeflt  )
   close( luwrite  )

   luprolog = luprolog - notmps
   lustart  = lustart  - notmps
   luloop   = luloop   - notmps
   luend    = luend    - notmps
   ludeflt  = ludeflt  - notmps
   luwrite  = luwrite  - notmps
end subroutine close_tmp_files

!===============================================================================
! APPEND_TMP_FILES --
!    Routine to append the contents of the temporary files
! Arguments:
!    lufirst           First LU-number to use
!===============================================================================

subroutine append_files( lufirst )
   integer, intent(in) :: lufirst

   integer             :: lu
   integer             :: io
   character(len=120)  :: line

   !
   ! If we have not written a subroutine yet, then
   ! now is the time to close the overall initialisation
   ! routine
   ! -- no longer needed
   !
   !if ( .not. contains ) then
   !   write( lusubs, '(a,/)' ) 'end subroutine'
   !   contains = .true.
   !endif

   !
   ! Copy the contents of the scratch files
   !
   do lu = lufirst,lufirst+notmps-1
      rewind( lu )
      do
         read( lu, '(a)', iostat=io ) line
         if ( io /= 0 ) exit
         write( lusubs, '(a)' ) trim(line)
      enddo
      rewind( lu )
   enddo
end subroutine append_files

!===============================================================================
! MERGE_FILES --
!    Routine to merge all temporary files into the definite file
! Arguments:
!    None
!===============================================================================

subroutine merge_files

   integer             :: io
   character(len=120)  :: line

   !
   ! Copy the contents of the "subroutines" file
   !
   rewind( lusubs )
   do
      read( lusubs, '(a)', iostat=io ) line
      if ( io /= 0 ) exit
      write( ludef, '(a)' ) trim(line)
   enddo

   !
   ! Copy the contents of the "write xml file subroutine" file
   !
   rewind( luwritp )
   do
      read( luwritp, '(a)', iostat=io ) line
      if ( io /= 0 ) exit
      write( ludef, '(a)' ) trim(line)
   enddo

   !
   ! Copy the contents of the "initialisation subroutine" file
   !
   rewind( luinit )
   do
      read( luinit, '(a)', iostat=io ) line
      if ( io /= 0 ) exit
      write( ludef, '(a)' ) trim(line)
   enddo
end subroutine merge_files

!===============================================================================
! WRITE_PROLOG --
!    Routine to write the beginning of the module
! Arguments:
!    None
!===============================================================================

subroutine write_prolog
   write( ludef, '(a)' ) &
  &   'module xml_data_' // trim(fname), &
  &   '   use READ_XML_PRIMITIVES', &
  &   '   use WRITE_XML_PRIMITIVES', &
  &   '   use XMLPARSE', &
  &   '   implicit none', &
  &   '   save', &
  &   '   integer, private :: lurep_', &
  &   '   logical, private :: strict_'

   write( luprolog, '(a)' ) &
  &   'subroutine read_xml_file_'//trim(fname)//'(fname, lurep, errout)' , &
  &   '   character(len=*), intent(in)           :: fname'     , &
  &   '   integer, intent(in), optional          :: lurep'     , &
  &   '   logical, intent(out), optional         :: errout'    , &
  &   '   '                                                    , &
  &   '   type(XML_PARSE)                        :: info'      , &
  &   '   logical                                :: error'     , &
  &   '   character(len=80)                      :: tag'       , &
  &   '   character(len=80)                      :: starttag'  , &
  &   '   logical                                :: endtag'    , &
  &   '   character(len=80), dimension(1:2,1:20) :: attribs'   , &
  &   '   integer                                :: noattribs' , &
  &   '   character(len=100000), dimension(1:1000)  :: data'      , &
  &   '   integer                                :: nodata'

   write( lusubs, '(a)' ) &
  &   'contains'
   write( luinit, '(a)' ) &
  &   'subroutine init_xml_file_'//trim(fname)

   write( luwritp, '(a)' ) &
  &   'subroutine write_xml_file_'//trim(fname)//'(fname, lurep)' , &
  &   '   character(len=*), intent(in)           :: fname'     , &
  &   '   integer, intent(in), optional          :: lurep'     , &
  &   '   '                                                    , &
  &   '   type(XML_PARSE)                        :: info'      , &
  &   '   integer                                :: indent = 0', &
  &   '   '                                                    , &
  &   '   call xml_open( info, fname, .false. )'               , &
  &   '   call xml_options( info, report_errors=.true.)'       , &
  &   '   if ( present(lurep) ) then'                          , &
  &   '       call xml_options( info, report_errors=.true.)'   , &
  &   '   endif'                                               , &
  &   '   write(info%lun,''(a)'') &'                           , &
  &   '      ''<' // trim(root_name) // '>'''

   write( lumain, '(a)' ) &
      &   '   ', &
      &   '   call init_xml_file_'//trim(fname), &
      &   '   call xml_open( info, fname, .true. )', &
      &   '   call xml_options( info, report_errors=.false., ignore_whitespace=.true.)', &
      &   '   lurep_ = 0', &
      &   '   if ( present(lurep) ) then', &
      &   '      lurep_ = lurep', &
      &   '      call xml_options( info, report_lun=lurep )', &
      &   '   endif', &
      &   '   do', &
      &   '      call xml_get( info, starttag, endtag, attribs, noattribs, &', &
      &   '         data, nodata)', &
      &   '      if ( starttag /= ''!--'' ) exit', &
      &   '   enddo', &
      &   '   if ( starttag /= "' // trim(root_name) // '" ) then', &
      &   '      call xml_report_errors( info, &', &
      &   '         ''XML-file should have root element "' // trim(root_name) // '"'')', &
      &   '      error = .true.', &
      &   '      call xml_close(info)', &
      &   '      return', &
      &   '   endif'

   call add_end_loop
end subroutine write_prolog

!===============================================================================
! ADD_BEGIN_LOOP --
!    Routine to write the start of the reading loop
! Arguments:
!    checktag        Whether code for checking the tag is required
!    component       Whether this is an ordinary variable or a component
!                    in a derived type
!===============================================================================

subroutine add_begin_loop( checktag, component )
   logical                                :: checktag
   logical                                :: component

   if ( component ) then
      write( luloop, '(a)' ) &
  &   '   call init_xml_type_'//trim(typename)//'(dvar)', &
  &   '   has_dvar = .true.'
   endif

   begin_loop = .false.

   if ( component ) then
      write( luloop, '(a)' ) &
     &   '   error  = .false.'     ,&
     &   '   att_   = 0'           ,&
     &   '   noatt_ = noattribs+1' ,&
     &   '   endtag_org = endtag'  ,&
     &   '   do',  &
     &   '      if ( nodata /= 0 ) then'           ,&
     &   '         noattribs = 0'                    ,&
     &   '         tag = starttag'                   ,&
     &   '      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then' ,&
     &   '         att_      = att_ + 1'             ,&
     &   '         if ( att_ <= noatt_-1 ) then'   ,&
     &   '            tag       = attribs(1,att_)'   ,&
     &   '            data(1)   = attribs(2,att_)'   ,&
     &   '            noattribs = 0'                 ,&
     &   '            nodata    = 1'                 ,&
     &   '            endtag    = .false.'           ,&
     &   '         else'                             ,&
     &   '            tag       = starttag'          ,&
     &   '            noattribs = 0'                 ,&
     &   '            nodata    = 0'                 ,&
     &   '            endtag    = .true.'            ,&
     &   '            cycle'                         ,&
     &   '         endif'                            ,&
     &   '      else', &
     &   '         if ( endtag_org ) then', &
     &   '            return', &
     &   '         else', &
     &   '            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )' ,&
     &   '            if ( xml_error(info) ) then'                 ,&
     &   '               write(lurep_,*) ''Error reading input file!''',&
     &   '               error = .true.'                           ,&
     &   '               return'                                   ,&
     &   '            endif'                                       ,&
     &   '         endif'                                          ,&
     &   '      endif'
   else
      write( luloop, '(a)' ) &
     &   '   error = .false.' ,&
     &   '   do',  &
     &   '      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )' ,&
     &   '      if ( xml_error(info) ) then'                       ,&
     &   '         write(lurep_,*) ''Error reading input file!'''  ,&
     &   '         error = .true.'                                 ,&
     &   '         return'                                         ,&
     &   '      endif'
   endif
   if ( checktag ) then
      write( luloop, '(a)' ) &
  &   '      if ( endtag .and. tag == starttag ) then'  ,&
  &   '         exit'                                     ,&
  &   '      endif'
   endif
   write( luloop, '(a)' ) &
  &   '      if ( endtag .and. noattribs == 0 ) then'   ,&
  &   '         if ( xml_ok(info) ) then'                 ,&
  &   '            cycle'                                 ,&
  &   '         else'                                     ,&
  &   '            exit'                                  ,&
  &   '         endif'                                    ,&
  &   '      endif'                                       ,&
  &   '      select case( tag )'
end subroutine add_begin_loop

!===============================================================================
! ADD_END_LOOP --
!    Routine to write the end of the reading loop
! Arguments:
!    None
!===============================================================================

subroutine add_end_loop

   write( luend, '(a)' ) &
  &   '      case (''comment'', ''!--'')' ,&
  &   '         ! Simply ignore', &
  &   '      case default' ,&
  &   '         if ( strict_ ) then', &
  &   '            error = .true.', &
  &   '            call xml_report_errors( info, &', &
  &   '               ''Unknown or wrongly placed tag: '' // trim(tag))',&
  &   '         endif'

   write( luend, '(a)' ) &
  &   '      end select' ,&
  &   '      nodata = 0' ,&
  &   '      if ( .not. xml_ok(info) ) exit' , &
  &   '   end do'
end subroutine add_end_loop

!===============================================================================
! ADD_VARIABLE --
!    Routine to write the definition of variables or components of
!    derived types
! Arguments:
!    component       Whether this is an ordinary variable or a component
!                    in a derived type
!===============================================================================

subroutine add_variable( component )
   logical                                :: component

   integer                                :: idx1
   integer                                :: idx2
   integer                                :: idx3
   integer                                :: idx4
   integer                                :: idx5
   integer                                :: idx6
   integer                                :: idx7
   integer                                :: k
   integer                                :: vdim

   character(len=32)                      :: varname
   character(len=40)                      :: varcomp
   character(len=40)                      :: varshape
   character(len=100)                     :: vardefault
   character(len=32)                      :: vartype
   character(len=32)                      :: vartag
   character(len=10)                      :: dim
   character(len=10)                      :: strlength
   character(len=32)                      :: initptr

   idx1 = xml_find_attrib( attribs, noattribs, 'name', varname )
   idx2 = xml_find_attrib( attribs, noattribs, 'type', vartype )
   strlength = "--"
   idx5 = xml_find_attrib( attribs, noattribs, 'length', strlength )
   if ( idx1 <= 0 ) then
      write( 20, * ) 'Variable/component found which has no name'
      error = .true.
   endif
   if ( idx2 <= 0 ) then
      write( 20, * ) 'Variable/component found which has no type - ',trim(varname)
      error = .true.
   else
      dim  = '--'
      idx7 = xml_find_attrib( attribs, noattribs, 'dimension', dim )
      idx6 = xml_find_attrib( attribs, noattribs, 'shape', varshape )
      if ( idx7 >= 1 ) then
         if ( dim == '1' ) then
            idx3 = xml_find_attrib( types, notypes, vartype, declare )
            if ( idx3 > notypes_predefined ) then
               vartype = trim(vartype) // '-array'
            else
               vartype = trim(vartype) // '-1dim'
            endif
         else
            error = .true.
            write(20,*) 'Dimension not supported: ',dim
         endif
      endif
      if ( idx6 >= 1 ) then
         vartype = trim(vartype) // '-shape'
         vdim    = 1
         if ( index(varshape, ',') > 0 ) then
             vdim = 2
         endif
      endif

      idx3 = xml_find_attrib( types, notypes, vartype, declare )
      if ( idx3 <= 0 ) then
         write( 20, * ) &
            'Variable/component with unknown type - ',trim(varname)
         error = .true.
      endif
   endif

   idx4 = xml_find_attrib( attribs, noattribs, 'default', vardefault )

   if ( component ) then
      varcomp = 'dvar%'//varname
   else
      varcomp = varname
   endif

   idx1 = xml_find_attrib( attribs, noattribs, 'tag', vartag )
   if ( idx1 < 1 ) then
      vartag = varname
   endif

   if ( .not. error ) then
      if ( index( declare, "pointer" ) > 0 ) then
         initptr = " => null()"
      else
         initptr = ""
      endif

      k = index( declare, 'SHAPE' )
      if ( k > 0 ) then
          declare = declare(1:k-1) // trim(varshape) // declare(k+5:)
      endif

      if ( index( declare, "?" ) <= 0 ) then
         write( ludef,    '(4a)' ) declare,    ' :: ', trim(varname), trim(initptr)
      else
         if ( strlength == "--" ) then
             strlength = "1" ! Hm, error is better?
         endif
         idx5 = index( declare, "?" )
         write( ludef,    '(6a)' ) declare(1:idx5-1), trim(strlength), declare(idx5+1:), &
            ' :: ', trim(varname), trim(initptr)
      endif

      if ( idx6 > 0 ) then
         k = index( types(2,idx3-1), '?' )
         if ( k <= 0 ) then
            write( luprolog, '(3a)' ) types(2,idx3-1), ' :: ', 'p_'//trim(varname)
         else
            write( luprolog, '(6a)' ) types(2,idx3-1)(1:k-1), trim(strlength), &
               types(2,idx3-1)(k+1:), ' :: ', 'p_'//trim(varname)
         endif
      endif
      write( luprolog, '(3a)' ) types(2,1), ' :: ', 'has_'//trim(varname)
      write( lustart,  '(3a)' ) '   has_', varname, ' = .false.'
      if ( dim /= '--' ) then
         write( lustart,  '(3a)' ) '   allocate(' // trim(varcomp), '(0))'
      endif
      write( luloop,   '(a)' ) '      case('''//trim(vartag)//''')'

      if ( idx6 <= 0 ) then
         write( luloop,   '(a)' ) &
            &'         call '//trim(types(3,idx3))//'( &', &
            &'            info, tag, endtag, attribs, noattribs, data, nodata, &',&
            &'            ' // trim(varcomp) // ', has_'//trim(varname) // ' )'
      else
         write( luloop,   '(a)' ) &
            &'         call '//trim(types(3,idx3))//'( &', &
            &'            info, tag, endtag, attribs, noattribs, data, nodata, &',&
            &'            p_' // trim(varname) // ', has_'//trim(varname) // ' )',&
            &'         if ( has_'//trim(varname) // ') then'
         if ( vdim == 1 ) then
             write( luloop,   '(a)' ) &
            &'            if ( size('//trim(varcomp)//') <= size(p_'//trim(varname)//') ) then', &
            &'               '//trim(varcomp) // ' = p_'//trim(varname)//'(1:size('//trim(varcomp)//'))', &
            &'            else', &
            &'               '//trim(varcomp) // '(1:size(p_'//trim(varname)//')) = p_'//trim(varname), &
            &'            endif'
         else
             write( luloop,   '(a)' ) &
            &'            if ( size(p_'//trim(varname)//') >= size('//trim(varcomp)//') ) then',&
            &'               '//trim(varcomp)//' = reshape(p_'//trim(varname)//', shape('//trim(varcomp)//'))',&
            &'            else',&
            &'               has_'//trim(varname)//' = .false.',&
            &'               call xml_report_errors(info, ''Incorrect number of values for '//trim(varname)//''')', &
            &'            endif'
         endif
         write( luloop,   '(a)' ) &
            &'            deallocate( p_'//trim(varname)//' )', &
            &'         endif'
      endif
      if ( idx4 <= 0 ) then
         write( luend,    '(a)' ) &
         &'   if ( .not. has_'//trim(varname)//' ) then'

         if ( component ) then
            write( luend,    '(a)' ) &
         &'      has_dvar = .false.'
         else
            write( luend,    '(a)' ) &
         &'      error = .true.'
         endif

         write( luend,    '(a)' ) &
         &'      call xml_report_errors(info, ''Missing data on '//trim(varname)//''')', &
         &'   endif'
      else
         !
         ! Note: the attribute value is supposed to have the quotes, if that
         ! is relevant for the variable's type
         !
         if ( component ) then
            write( ludeflt,    '(4a)' ) &
            &'   dvar%', trim(varname), ' = ', attribs(2,idx4)
         else
            write( luinit,   '(4a)' ) &
            &'   ', trim(varname), ' = ', attribs(2,idx4)
         endif
      endif
      !
      ! Write the component/variable
      !
      if ( component ) then
         write( luwrite,   '(4a)' ) &
         &'   call '//trim(types(4,idx3))//'(', &
         &' info, '''//trim(vartag)//''', indent+3, dvar%', trim(varname), ')'
      else
         write( luwrv,   '(4a)' ) &
         &'   call '//trim(types(4,idx3))//'(', &
         &' info, '''//trim(vartag)//''', indent+3, ', trim(varname), ')'
      endif
   endif

end subroutine add_variable

!===============================================================================
! ADD_TYPEDEF --
!    Routine to write the definition and other code for a derived type
! Arguments:
!    dyn_strings     Whether dynamic string lengths are allowed
!===============================================================================

subroutine add_typedef( dyn_strings )
   logical, intent(in)                    :: dyn_strings

   integer                                :: idx1
   integer                                :: idx2

   character(len=32)                      :: typetag

   idx1 = xml_find_attrib( attribs, noattribs, 'name', typename )
   if ( idx1 <= 0 ) then
      write( 20, * ) 'Type definition found which has no name'
      error = .true.
   endif

   !
   ! We need a new set of temporary files
   !
   call open_tmp_files( luprolog+notmps )

   idx2 = xml_find_attrib( attribs, noattribs, 'tag', typetag )
   if ( idx1 < 1 ) then
      typetag = typename
   endif

   if ( .not. error ) then
      write( ludef,    '(/,2a)' ) 'type ',trim(typename)
      write( luprolog, '(a)' ) &
  &   'subroutine read_xml_type_'//trim(typename)//'_array( &'      ,&
  &   '      info, tag, endtag, attribs, noattribs, data, nodata, &',&
  &   '      dvar, has_dvar )'                                      ,&
  &   '   type(XML_PARSE)                                 :: info'  ,&
  &   '   character(len=*), intent(inout)                 :: tag    ',&
  &   '   logical, intent(inout)                          :: endtag ',&
  &   '   character(len=*), dimension(:,:), intent(inout) :: attribs',&
  &   '   integer, intent(inout)                          :: noattribs',&
  &   '   character(len=*), dimension(:), intent(inout)   :: data   ',&
  &   '   integer, intent(inout)                          :: nodata ',&
  &   '   type('//trim(typename)//'), dimension(:), pointer :: dvar ',&
  &   '   logical, intent(inout)                       :: has_dvar  ',&
  &   '   '                                                          ,&
  &   '   integer                                      :: newsize   ',&
  &   '   type('//trim(typename)//'), dimension(:), pointer :: newvar',&
  &   '   '                                                          ,&
  &   '   newsize = size(dvar) + 1'                                  ,&
  &   '   allocate( newvar(1:newsize) )'                             ,&
  &   '   newvar(1:newsize-1) = dvar'                                ,&
  &   '   deallocate( dvar )'                                        ,&
  &   '   dvar => newvar'                                            ,&
  &   '   '                                                          ,&
  &   '   call read_xml_type_'//trim(typename)// &
  &         '( info, tag, endtag, attribs, noattribs, data, nodata, &',&
  &   '              dvar(newsize), has_dvar )'                      ,&
  &   'end subroutine read_xml_type_'//trim(typename)//'_array'      ,&
  &   '   '

      write( luwrite, '(a)' ) &
  &   'subroutine write_xml_type_'//trim(typename)//'_array( &'     ,&
  &   '      info, tag, indent, dvar )'                             ,&
  &   '   type(XML_PARSE)                                 :: info'  ,&
  &   '   character(len=*), intent(in)                    :: tag'   ,&
  &   '   integer                                         :: indent',&
  &   '   type('//trim(typename)//'), dimension(:)        :: dvar'  ,&
  &   '   integer                                         :: i'     ,&
  &   '   do i = 1,size(dvar)'                                      ,&
  &   '       call write_xml_type_'//trim(typename)// &
  &              '( info, tag, indent, dvar(i) )'                   ,&
  &   '   enddo'                                                    ,&
  &   'end subroutine write_xml_type_'//trim(typename)//'_array'    ,&
  &   '   ', &
  &   'subroutine write_xml_type_'//trim(typename)//'( &'           ,&
  &   '      info, tag, indent, dvar )'                             ,&
  &   '   type(XML_PARSE)                                 :: info'  ,&
  &   '   character(len=*), intent(in)                    :: tag'   ,&
  &   '   integer                                         :: indent',&
  &   '   type('//trim(typename)//')                      :: dvar'  ,&
  &   '   character(len=100)                              :: indentation' ,&
  &   '   indentation = '' '''                                      ,&
  &   '   write(info%lun, ''(4a)'' ) indentation(1:min(indent,100)),&',&
  &   '       ''<'',trim(tag), ''>'''

      write( luprolog, '(a)' ) &
  &   'subroutine read_xml_type_'//trim(typename)//&
  &         '( info, starttag, endtag, attribs, noattribs, data, nodata, &' ,&
  &   '              dvar, has_dvar )'                              ,&
  &   '   type(XML_PARSE)                                 :: info'   ,&
  &   '   character(len=*), intent(in)                    :: starttag',&
  &   '   logical, intent(inout)                          :: endtag  ',&
  &   '   character(len=*), dimension(:,:), intent(inout) :: attribs',&
  &   '   integer, intent(inout)                          :: noattribs',&
  &   '   character(len=*), dimension(:), intent(inout)   :: data   ',&
  &   '   integer, intent(inout)                          :: nodata ',&
  &   '   type('//trim(typename)//'), intent(inout)  :: dvar' ,&
  &   '   logical, intent(inout)                       :: has_dvar  ',&
  &   '   '                                                          ,&
  &   '   integer                                      :: att_      ',&
  &   '   integer                                      :: noatt_    ',&
  &   '   logical                                      :: error     ',&
  &   '   logical                                      :: endtag_org'
      if ( dyn_strings ) then
          write( luprolog, '(a)' ) &
  &   '   character(len=len(starttag))                 :: tag       '
      else
          write( luprolog, '(a)' ) &
  &   '   character(len=80)                            :: tag       '
      endif

      !
      ! Note: this may require a more sophisticated approach
      ! when the components of the type are also pointers ...
      !
      write( ludeflt, '(a)' ) &
  &   'subroutine init_xml_type_'//trim(typename)//'_array( dvar )  ',&
  &   '   type('//trim(typename)//'), dimension(:), pointer :: dvar ',&
  &   '   if ( associated( dvar ) ) then'                            ,&
  &   '      deallocate( dvar )'                                     ,&
  &   '   endif'                                                     ,&
  &   '   allocate( dvar(0) )'                                       ,&
  &   'end subroutine init_xml_type_'//trim(typename)//'_array'      ,&
  &   'subroutine init_xml_type_'//trim(typename)//'(dvar)'          ,&
  &   '   type('//trim(typename)//') :: dvar '

      begin_loop = .true.

      call add_end_loop

      !
      ! Add the names of the two new types to the list
      !
      allocate( new_types(1:4,1:notypes+3) )
      new_types(:,1:notypes) = types
      deallocate( types )
      types => new_types

      types(1,notypes+1) = typename
      types(2,notypes+1) = '   type('//trim(typename)//')'
      types(3,notypes+1) = 'read_xml_type_'//trim(typename)
      types(4,notypes+1) = 'write_xml_type_'//trim(typename)

      types(1,notypes+2) = trim(typename) // '-array'
      types(2,notypes+2) = '   type('//trim(typename)//'), dimension(:), pointer'
      types(3,notypes+2) = 'read_xml_type_'//trim(typename)//'_array'
      types(4,notypes+2) = 'write_xml_type_'//trim(typename)//'_array'

      types(1,notypes+3) = trim(typename) // '-shape'
      types(2,notypes+3) = '   type('//trim(typename)//'), dimension(SHAPE)'
      types(3,notypes+3) = 'read_xml_type_'//trim(typename)//'_array'
      types(4,notypes+3) = 'write_xml_type_'//trim(typename)//'_array'

      notypes = notypes + 3

   endif
end subroutine add_typedef

!===============================================================================
! CLOSE_TYPEDEF --
!    Routine to write the last code fragments for a derived type
! Arguments:
!    component       Turn off the "component" parameter
!===============================================================================

subroutine close_typedef( component )
   logical, intent(out) :: component

   component = .false.
   write( ludef, '(a)' ) 'end type '//trim(typename)
   write( luend, '(a)' ) &
  &   'end subroutine read_xml_type_'//trim(typename)
   write( ludeflt, '(a)' ) &
  &   'end subroutine init_xml_type_'//trim(typename)
   write( luwrite, '(a)' ) &
  &   '   write(info%lun,''(4a)'') indentation(1:min(indent,100)), &' ,&
  &   '       ''</'' //trim(tag) // ''>''', &
  &   'end subroutine write_xml_type_'//trim(typename) ,&
  &   ' '
   call append_files( luprolog )
   call close_tmp_files

end subroutine close_typedef

!===============================================================================
! ADD_PLACEHOLDER --
!    Routine to write the starting code fragments for a placeholder tag
! Arguments:
!    dyn_strings     Whether dynamic string lengths are allowed
!===============================================================================

subroutine add_placeholder( dyn_strings )
   logical, intent(in)                    :: dyn_strings

   integer                                :: idx1
   integer                                :: idx2

   character(len=32)                      :: tag
   character(len=20)                      :: optional

   idx1 = xml_find_attrib( attribs, noattribs, 'tag', tag )
   if ( idx1 <= 0 ) then
      write( 20, * ) 'Placeholder definition found which has no tag name'
      error = .true.
   endif

   optional = 'no'
   idx2 = xml_find_attrib( attribs, noattribs, 'optional', optional )

   if ( optional == 'yes' ) then
      if ( begin_loop ) then
         call add_begin_loop( .false., .false. )
      endif
      write( luloop, '(a)' ) &
         '      case('''//trim(tag)//''')',&
         '         ! Simply ignore the tag'
   else
      no_placeholders = no_placeholders + 1
      placeholder(no_placeholders) = tag

      write( luloop, '(a)' ) &
         '      case('''//trim(tag)//''')',&
        &'         call read_xml_place_'//trim(tag)//'( info, &', &
        &'            tag, attribs, noattribs, data, nodata )'
      comp = .false.

      !
      ! We need a new set of temporary files
      !
      call open_tmp_files( luprolog+notmps )

      !
      ! Write the first part of the routine
      ! NOTE:
      ! Will require an extra argument when collecting all variables
      ! in one derived type
      !
      write( luprolog, '(a)' ) &
  &   'subroutine read_xml_place_'//trim(tag)//&
  &         '( info, starttag, attribs, noattribs, data, nodata )'   ,&
  &   '   type(XML_PARSE)                                 :: info'   ,&
  &   '   character(len=*), intent(in)                    :: starttag',&
  &   '   character(len=*), dimension(:,:), intent(inout) :: attribs',&
  &   '   integer, intent(inout)                          :: noattribs',&
  &   '   character(len=*), dimension(:), intent(inout)   :: data   ',&
  &   '   integer, intent(inout)                          :: nodata ',&
  &   '   '                                                          ,&
  &   '   logical                                      :: error     ',&
  &   '   logical                                      :: endtag    '
      if ( dyn_strings ) then
          write( luprolog, '(a)' ) &
  &   '   character(len=len(starttag))                 :: tag       '
      else
          write( luprolog, '(a)' ) &
  &   '   character(len=80)                            :: tag       '
      endif

      begin_loop = .true.

      call add_end_loop

      write(luwrite,'(a)') &
        &'subroutine write_xml_place_'//trim(tag)//'( &'        ,&
        &'      info, indent )'                                 ,&
        &'   type(XML_PARSE)                          :: info'  ,&
        &'   integer                                  :: indent',&
        &'   character(len=100)                       :: indentation' ,&
        &'   indentation = '' '''                               ,&
        &'   write(info%lun, ''(4a)'' ) indentation(1:min(indent,100)),&',&
        &'       ''<'// trim(tag) // '>'''
      write(luwritp,'(a)') &
        &'   call write_xml_place_'//trim(tag)//'( info, indent+3 )'

   endif
end subroutine add_placeholder

!===============================================================================
! CLOSE_PLACEHOLDER --
!    Routine to write the last code fragments for a placeholder
! Arguments:
!    None
!===============================================================================

subroutine close_placeholder

   write( luend, '(a)' ) &
  &   'end subroutine read_xml_place_'//trim(placeholder(no_placeholders))

   write( luwrite, '(a)' ) &
  &   '   write(info%lun,''(4a)'') indentation(1:min(indent,100)), &' ,&
  &   '       ''</' //trim(placeholder(no_placeholders)) // '>''' ,&
  &   'end subroutine write_xml_place_'//trim(placeholder(no_placeholders)) ,&
  &   ' '

   call append_files( luprolog )
   call close_tmp_files

   no_placeholders = no_placeholders - 1

end subroutine close_placeholder

end program
