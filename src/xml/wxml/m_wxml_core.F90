module m_wxml_core

#ifndef DUMMYLIB
  use fox_m_fsys_abort_flush, only: pxfabort
  use fox_m_fsys_array_str, only: vs_str, str_vs, vs_str_alloc
  use fox_m_fsys_string, only: toLower
  use fox_m_utils_uri, only: URI, parseURI, destroyURI
  use m_common_attrs, only: dictionary_t, getLength, get_key, get_value, &
    hasKey, add_item_to_dict, init_dict, reset_dict, destroy_dict, &
    getWhitespaceHandling, sortAttrs
  use m_common_buffer, only: buffer_t, len, add_to_buffer, reset_buffer, &
    dump_buffer
  use m_common_charset, only: XML1_0, XML1_1, checkChars
  use m_common_element, only: parse_dtd_element, parse_dtd_attlist
  use m_common_elstack, only: elstack_t, len, get_top_elstack, pop_elstack, &
    is_empty, init_elstack, push_elstack, destroy_elstack
  use m_common_entities, only: existing_entity, is_unparsed_entity
  use m_common_error, only: FoX_warning_base, FoX_error_base, FoX_fatal_base, &
    error_stack, in_error, FoX_get_fatal_errors, FoX_get_fatal_warnings
  use m_common_io, only: get_unit
  use m_common_namecheck, only: checkEncName, checkName, checkQName, &
    checkCharacterEntityReference, checkPublicId, prefixOfQName, &
    localpartofQName, checkPEDef, checkPseudoAttValue, checkAttValue, checkNCName, &
    likeCharacterEntityReference, checkCharacterEntityReference
  use m_common_namespaces, only: namespaceDictionary, getnamespaceURI, &
  initnamespaceDictionary, addDefaultNS, destroyNamespaceDictionary, &
  addPrefixedNS, isPrefixInForce, checkNamespacesWriting, checkEndNamespaces
  use m_common_notations, only: add_notation, notation_exists
  use m_common_struct, only: xml_doc_state, init_xml_doc_state, destroy_xml_doc_state, &
    register_internal_PE, register_external_PE, register_internal_GE, register_external_GE
  use m_wxml_escape, only: escape_string
#ifdef PGF90
  use m_common_element, only : element_t
#endif
#endif

  implicit none
  private

#ifndef DUMMYLIB
  integer, parameter :: indent_inc = 2
  ! TOHW should we let this be set?

  !Output State Machines
  ! status wrt root element:
  integer, parameter :: WXML_STATE_1_JUST_OPENED = 0 
  !File is just opened, nothing written to it yet. 
  integer, parameter :: WXML_STATE_1_BEFORE_ROOT = 1
  !File has been opened, something has been written, but no root element yet.
  integer, parameter :: WXML_STATE_1_DURING_ROOT = 2
  !The root element has been opened but not closed
  integer, parameter :: WXML_STATE_1_AFTER_ROOT = 3
  !The root element has been opened but not closed

  ! status wrt tags:
  integer, parameter :: WXML_STATE_2_OUTSIDE_TAG = 0
  !We are not within a tag.
  integer, parameter :: WXML_STATE_2_INSIDE_PI = 1
  !We are inside a Processing Instruction tag
  integer, parameter :: WXML_STATE_2_INSIDE_ELEMENT = 2
  !We are inside an element tag.
  integer, parameter :: WXML_STATE_2_IN_CHARDATA = 3
  !We are inside deliberately-constructed text. (this is only necessary for preserve_whitespace)

  ! status wrt DTD
  integer, parameter :: WXML_STATE_3_BEFORE_DTD = 0
  ! No DTD has been encountered yet.
  integer, parameter :: WXML_STATE_3_DURING_DTD = 1
  ! Halfway throught outputting a DTD
  integer, parameter :: WXML_STATE_3_INSIDE_INTSUBSET = 2
  !We are inside the internal subset definition
  integer, parameter :: WXML_STATE_3_AFTER_DTD = 3
  ! Finished outputting a DTD
#endif


  type xmlf_t
    private
#ifdef DUMMYLIB
    integer :: i = 0
#else
    type(xml_doc_state) :: xds
    integer                   :: lun = -1
    type(buffer_t)            :: buffer
    type(elstack_t)           :: stack
    type(dictionary_t)        :: dict
    integer                   :: state_1 = -1
    integer                   :: state_2 = -1
    integer                   :: state_3 = -1
    logical                   :: minimize_overrun = .true.
    logical                   :: pretty_print = .false.
    logical                   :: canonical = .false.
    integer                   :: indent = 0
    character, pointer        :: name(:)
    logical                   :: namespace = .true.
    type(namespaceDictionary) :: nsDict
#endif
  end type xmlf_t

  public :: xmlf_t

  public :: xml_OpenFile
  public :: xml_NewElement
  public :: xml_EndElement
  public :: xml_Close
  public :: xml_AddXMLDeclaration
  public :: xml_AddXMLStylesheet
  public :: xml_AddXMLPI
  public :: xml_AddComment
  public :: xml_AddCharacters
  public :: xml_AddNewline
  public :: xml_AddEntityReference
  public :: xml_AddAttribute
  public :: xml_AddPseudoAttribute
  public :: xml_DeclareNamespace
  public :: xml_UnDeclareNamespace
  public :: xml_AddDOCTYPE
  public :: xml_AddParameterEntity
  public :: xml_AddInternalEntity
  public :: xml_AddExternalEntity
  public :: xml_AddNotation
  public :: xml_AddElementToDTD
  public :: xml_AddAttlistToDTD
  public :: xml_AddPEreferenceToDTD

  public :: xmlf_Name
  public :: xmlf_OpenTag

  public :: xmlf_SetPretty_print
  public :: xmlf_GetPretty_print

  interface xml_AddCharacters
    module procedure xml_AddCharacters_Ch
  end interface
  interface xml_AddAttribute
    module procedure xml_AddAttribute_Ch
  end interface
  interface xml_AddPseudoAttribute
    module procedure xml_AddPseudoAttribute_Ch
  end interface
 
#ifndef DUMMYLIB
  !overload error handlers to allow file info
  interface wxml_warning
    module procedure wxml_warning_xf, FoX_warning_base
  end interface
  interface wxml_error
    module procedure wxml_error_xf, FoX_error_base
  end interface
  interface wxml_fatal
    module procedure wxml_fatal_xf, FoX_fatal_base
  end interface

  ! Heuristic (approximate) target for justification of output
  ! only gets used for outputting attributes
  integer, parameter  :: COLUMNS = 80

  ! TOHW - This is the longest string that may be output without
  ! a newline. The buffer must not be larger than this, but its size 
  ! can be tuned for performance.
  !lowest value found so far is 4096, for NAG. We use 1024 just in case.
  integer, parameter  :: xml_recl = 1024
#endif

contains

  subroutine xml_OpenFile(filename, xf, unit, iostat, preserve_whitespace, &
    pretty_print, minimize_overrun, canonical, replace, addDecl, warning, &
    validate, namespace)
    character(len=*), intent(in)  :: filename
    type(xmlf_t), intent(inout)   :: xf
    integer, intent(in), optional :: unit
    integer, intent(out), optional :: iostat
    logical, intent(in), optional :: preserve_whitespace
    logical, intent(in), optional :: pretty_print
    logical, intent(in), optional :: minimize_overrun
    logical, intent(in), optional :: canonical
    logical, intent(in), optional :: replace
    logical, intent(in), optional :: addDecl
    logical, intent(in), optional :: warning
    logical, intent(in), optional :: validate
    logical, intent(in), optional :: namespace
    
#ifdef DUMMYLIB
    if (present(iostat)) iostat = 0
#else
    logical :: repl, decl
    integer :: iostat_

    if (xf%lun /= -1) &
      call wxml_fatal("Trying to reopen an already-open XML file")
    
    if (present(replace)) then
      repl = replace
    else
      repl = .true.
    endif
    if (present(addDecl)) then
      decl = addDecl
    else
      decl = .true.
    endif
    if (present(iostat)) iostat = 0

    allocate(xf%name(0))
    
    if (present(unit)) then
      if (unit==-1) then
        call get_unit(xf%lun,iostat_)
        if (iostat_ /= 0) then
          if (present(iostat)) iostat = iostat_
          return
        endif
      else
        xf%lun = unit
      endif
    else
      call get_unit(xf%lun,iostat_)
      if (iostat_ /= 0) then
        if (present(iostat)) iostat = iostat_
        return
      endif
    endif

    ! Use large I/O buffer in case the O.S./Compiler combination
    ! has hard-limits by default (i.e., NAGWare f95's 1024 byte limit)
    ! This is related to the maximum size of the buffer.
    ! TOHW - This is the longest string that may be output without
    ! a newline. The buffer must not be larger than this, but its size 
    ! can be tuned for performance.
    
    if (repl) then
      ! NAG insists on unnecessary duplication of iostat etc here
      if (present(iostat)) then
        open(unit=xf%lun, file=filename, form="formatted", status="replace", &
          action="write", recl=xml_recl, iostat=iostat)
      else
        open(unit=xf%lun, file=filename, form="formatted", status="replace", &
          action="write", recl=xml_recl)
      endif
    else 
      if (present(iostat)) then
        open(unit=xf%lun, file=filename, form="formatted", status="new", &
          action="write", recl=xml_recl, iostat=iostat)
      else
        open(unit=xf%lun, file=filename, form="formatted", status="new", &
          action="write", recl=xml_recl)
      endif
    endif

    call init_elstack(xf%stack)

    call init_dict(xf%dict)
    !NB it can make no difference which XML version we are using
    !until after we output the XML declaration. So we set it to
    !1.0 for the moment & reset below.
    ! Actually, this is done automatically in initializing xf%xds
    call init_xml_doc_state(xf%xds)
    xf%xds%documentURI => vs_str_alloc(filename)

    if (present(warning)) then
      xf%xds%warning = warning
    else
      xf%xds%warning = .false.
    endif
    if (present(validate)) then
      xf%xds%valid = validate
    else
      xf%xds%valid = .false.
    endif
    xf%state_1 = WXML_STATE_1_JUST_OPENED
    xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    xf%state_3 = WXML_STATE_3_BEFORE_DTD
    
    if (present(pretty_print)) then
      xf%pretty_print = pretty_print
    else
      xf%pretty_print = .true.
    endif
    if (present(minimize_overrun)) then
      xf%minimize_overrun = minimize_overrun
    else
      xf%minimize_overrun = .false.
    endif
    if (present(preserve_whitespace)) then
      xf%pretty_print = .not.preserve_whitespace
      xf%minimize_overrun = preserve_whitespace
    endif
    if (present(canonical)) then
      xf%canonical = canonical
    else
      xf%canonical = .false.
    endif
! FIXME interplay of above options
      
    xf%indent = 0
    
    if (decl) then
      call xml_AddXMLDeclaration(xf,encoding='UTF-8')
      ! which will reset the buffer itself
    else
      call reset_buffer(xf%buffer, xf%lun, xf%xds%xml_version)
    endif

    if (present(namespace)) then
      xf%namespace = namespace
    else
      xf%namespace = .true.
    endif
    if (xf%namespace) &
      call initNamespaceDictionary(xf%nsDict)
#endif
  end subroutine xml_OpenFile


  subroutine xml_AddXMLDeclaration(xf, version, encoding, standalone)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in), optional :: version
    character(len=*), intent(in), optional :: encoding
    logical, intent(in), optional :: standalone

#ifndef DUMMYLIB
    call check_xf(xf)
    ! Don't need to call checkChars on args, everything is checked
    ! fully below anyway.
    
    if (xf%state_1 /= WXML_STATE_1_JUST_OPENED) &
      call wxml_error("Tried to put XML declaration in wrong place")
    
    call reset_buffer(xf%buffer, xf%lun, xf%xds%xml_version)
    
    call xml_AddXMLPI(xf, "xml", xml=.true.)
    if (present(version)) then
      if (version =="1.0") then
        xf%xds%xml_version = XML1_0
        call xml_AddPseudoAttribute(xf, "version", version)
      elseif (version=="1.1") then
        xf%xds%xml_version = XML1_1
        call xml_AddPseudoAttribute(xf, "version", version)
      else
        call wxml_error("Invalid XML version.")
      endif
    else
      call xml_AddPseudoAttribute(xf, "version", "1.0")
      xf%xds%xml_version = XML1_0
    endif
    if (present(encoding)) then
      if (.not.checkEncName(encoding)) &
        call wxml_error("Invalid encoding name: "//encoding)
      if (encoding /= 'UTF-8' .and. encoding /= 'utf-8') &
        call wxml_warning(xf, "Non-default encoding specified: "//encoding)
      call xml_AddPseudoAttribute(xf, "encoding", encoding)
    endif
    if (present(standalone)) then
      xf%xds%standalone_declared = .true.
      xf%xds%standalone = standalone
      if (standalone) then
        call xml_AddPseudoAttribute(xf, "standalone", "yes")
      else
        call xml_AddPseudoAttribute(xf, "standalone", "no")
      endif
    endif
    call close_start_tag(xf)
    ! We have to close explicitly here to ensure nothing gets tied
    ! up in the XML declaration
    xf%state_1 = WXML_STATE_1_BEFORE_ROOT
#endif
  end subroutine xml_AddXMLDeclaration


  subroutine xml_AddDOCTYPE(xf, name, system, public)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: system, public

#ifndef DUMMYLIB
    type(URI), pointer :: URIref

    call check_xf(xf)
    
    if (xf%namespace) then
      if (.not.checkQName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Name in DTD "//name)
    endif
    
    if (present(system)) then
      URIref => parseURI(system)
      if (.not.associated(URIref)) call wxml_error("xml_AddDOCTYPE: Invalid SYSTEM URI")
      call destroyURI(URIref)
    endif
    if (present(public)) then
      if (.not.checkPublicId(public)) call wxml_error("xml_AddDOCTYPE: Invalid PUBLIC ID")
    endif
    
    if (present(public).and..not.present(system)) &
      call wxml_error('xml_AddDOCTYPE: PUBLIC supplied without SYSTEM for: '//name)
    
    ! By having an external ID we probably render this non-standalone (unless we've said that it is in the declaration)
    if (present(system).and..not.xf%xds%standalone_declared) &
      xf%xds%standalone=.false.
    
    call close_start_tag(xf)
    
    if (xf%state_1 /= WXML_STATE_1_BEFORE_ROOT) &
      call wxml_error("Tried to put XML DOCTYPE in wrong place: "//name)
    
    if (xf%state_3 /= WXML_STATE_3_BEFORE_DTD) then
      call wxml_error("Tried to output more than one DOCTYPE declaration: "//name)
    else
      xf%state_3 = WXML_STATE_3_DURING_DTD
    endif
    
    call add_eol(xf)
    call add_to_buffer("<!DOCTYPE "//name, xf%buffer, .false.)
    
    deallocate(xf%name)
    allocate(xf%name(len(name)))
    xf%name = vs_str(name)
    
    if (present(system)) then
      if (present(public)) then
        call add_to_buffer(" PUBLIC", xf%buffer, .false.)
        call add_to_buffer(" """//public//"""", xf%buffer, .true.)
      else
        call add_to_buffer(" SYSTEM", xf%buffer, .false.)
      endif
      if (scan(system, """")/=0) then
        call add_to_buffer(" '"//system//"'", xf%buffer, .true.)
      else
        call add_to_buffer(" """//system//"""", xf%buffer, .true.)
      endif
    endif
#endif
  end subroutine xml_AddDOCTYPE


  subroutine xml_AddParameterEntity(xf, name, PEdef, system, public)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: PEDef
    character(len=*), intent(in), optional :: system
    character(len=*), intent(in), optional :: public

#ifndef DUMMYLIB
    type(URI), pointer :: URIref
#ifdef PGF90
    type(URI), pointer :: nullURIref
#endif
    call check_xf(xf)
#ifdef PGF90
    nullURIref => null()
#endif
    if (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
         call wxml_error("Invalid Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Name in DTD "//name)
    endif

    if (present(PEDef)) then
      if (.not.checkChars(PEDef,xf%xds%xml_version)) call wxml_error("xml_AddParameterEntity: Invalid character in PEDef")
    endif

    if (present(system)) then
      URIref => parseURI(system)
      if (.not.associated(URIref)) call wxml_error("xml_AddParameterEntity: Invalid SYSTEM URI")
      call destroyURI(URIref)
    endif
    if (present(public)) then
      if (.not.checkPublicId(public)) call wxml_error("xml_AddParameterEntity: Invalid PUBLIC ID")
    endif

    ! By adding a parameter entity (internal or external) we make this
    ! a non-standalone document.
    if (.not.xf%xds%standalone_declared) &
      xf%xds%standalone = .false.

    if (xf%state_3 == WXML_STATE_3_DURING_DTD) then
      call add_to_buffer(" [", xf%buffer, .false.)
      xf%state_3 = WXML_STATE_3_INSIDE_INTSUBSET
    endif

    if (xf%state_3 /= WXML_STATE_3_INSIDE_INTSUBSET) &
      call wxml_fatal("Cannot define Parameter Entity here: "//name)
      
    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) then
      call close_start_tag(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    endif

    if (present(PEdef)) then
      if (present(system) .or. present(public)) &
        call wxml_fatal("Parameter entity "//name//" cannot have both a PEdef and an External ID")
    else
      if (.not.present(system)) &
        call wxml_fatal("Parameter entity "//name//" must have either a PEdef or an External ID")
    endif
    if (present(PEdef)) then
      if (.not.checkPEDef(PEDef, xf%xds%xml_version)) &
        call wxml_fatal("Parameter entity definition is invalid: "//PEDef)
      if (xf%xds%standalone) then
        if (.not.checkExistingRefs()) &
          call wxml_error("Tried to reference unregistered parameter entity")
      else
        if (.not.checkExistingRefs()) &
          call wxml_warning(xf, "Reference to unknown parameter entity")
      endif
#ifdef PGF90
      call register_internal_PE(xf%xds, name=name, text=PEdef, baseURI=nullURIref, wfc=.false.)
#else
      call register_internal_PE(xf%xds, name=name, text=PEdef, baseURI=null(), wfc=.false.)
#endif

    else
#ifdef PGF90
      call register_external_PE(xf%xds, name=name, systemId=system, &
        publicId=public, baseURI=nullURIref, wfc=.false.)
#else
      call register_external_PE(xf%xds, name=name, systemId=system, &
        publicId=public, baseURI=null(), wfc=.false.)
#endif
    endif

    call add_eol(xf)

    call add_to_buffer("<!ENTITY % "//name, xf%buffer, .false.) ! name can never contain whitespace
    if (present(PEdef)) then
      if (index(PEdef, """") > 0) then ! FIXME what if PEdef has both " and ' in it
        call add_to_buffer(" '"//PEdef//"'", xf%buffer, .true.)
      else
        call add_to_buffer(" """//PEdef//"""", xf%buffer, .true.)
      endif
        call add_to_buffer(">", xf%buffer, .false.)
    else
      if (present(public)) then
        call add_to_buffer(" PUBLIC", xf%buffer, .false.)
        call add_to_buffer(" """//public//"""", xf%buffer, .true.)
      else
        call add_to_buffer(" SYSTEM", xf%buffer, .false.)
      endif
      if (scan(system, """")/=0) then
        call add_to_buffer(" '"//system//"'", xf%buffer, .true.)
      else
        call add_to_buffer(" """//system//"""", xf%buffer, .true.)
      endif
      call add_to_buffer(">", xf%buffer)
    endif

  contains
    function checkExistingRefs() result(p)
      logical :: p
      
      integer :: i1, i2
      
      ! Here we assume we have syntactic well-formedness as
      ! checked by checkPEDef.
      
      p = .false.
      i1 = index(PEDef, '%')
      i2 = 0 
      do while (i1 > 0)
        i1 = i2 + i1
        i2 = index(PEDef(i1+1:),';')
        if (i2 == 0) return
        i2 = i1 + i2
        if (.not.existing_entity(xf%xds%PEList, PEDef(i1+1:i2-1))) &
          return
        i1 = index(PEDef(i2+1:), '%')
      enddo
      p = .true.
      
    end function checkExistingRefs
#endif
  end subroutine xml_AddParameterEntity


  subroutine xml_AddInternalEntity(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

#ifndef DUMMYLIB
#ifdef PGF90
    type(URI), pointer :: nullURI
    nullURI => null()
#endif
    call check_xf(xf)

    if (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
         call wxml_error("Invalid Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Name in DTD "//name)
    endif

    if (.not.checkChars(value, xf%xds%xml_version)) call wxml_error("xml_AddInternalEntity: Invalid character in value")

    if (xf%state_3 == WXML_STATE_3_DURING_DTD) then
      call add_to_buffer(" [", xf%buffer)
      xf%state_3 = WXML_STATE_3_INSIDE_INTSUBSET
    endif

    if (xf%state_3 /= WXML_STATE_3_INSIDE_INTSUBSET) &
      call wxml_fatal("Cannot define Entity here: "//name)
      
    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) then
      call close_start_tag(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    endif

    if (.not.checkName(name, xf%xds%xml_version)) &
      call wxml_error("xml_AddInternalEntity: Invalid Name: "//name)
#ifdef PGF90
    call register_internal_GE(xf%xds, name=name, text=value, baseURI=nullURI, wfc=.false.)
#else
    call register_internal_GE(xf%xds, name=name, text=value, baseURI=null(), wfc=.false.)
#endif

    call add_eol(xf)
    
    !FIXME - valid entity values?
    call add_to_buffer("<!ENTITY "//name//" ", xf%buffer, .false.) ! name cannot contain whitespace
    if (index(value, """") > 0) then
      call add_to_buffer("'"//value//"'>", xf%buffer, .true.)
    else
      call add_to_buffer(""""//value//""">", xf%buffer, .true.)
    endif
#endif
  end subroutine xml_AddInternalEntity


  subroutine xml_AddExternalEntity(xf, name, system, public, notation)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: system
    character(len=*), intent(in), optional :: public
    character(len=*), intent(in), optional :: notation

#ifndef DUMMYLIB
    type(URI), pointer :: URIref
#ifdef PGF90
    type(URI), pointer :: nullURI
    nullURI => null()
#endif
    call check_xf(xf)

    if (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
         call wxml_error("Invalid Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Name in DTD "//name)
    endif
    URIref => parseURI(system)
    if (.not.associated(URIref)) call wxml_error("xml_AddExternalEntity: Invalid SYSTEM URI")
    call destroyURI(URIref)
    if (present(public)) then
      if (.not.checkPublicId(public)) call wxml_error("xml_AddExternalEntity: Invalid PUBLIC ID")
    endif
    if (present(notation)) then
      if (xf%namespace) then
        if (.not.checkNCName(notation, xf%xds%xml_version)) &
          call wxml_error("Invalid Name in DTD "//name)
      else
        if (.not.checkName(notation, xf%xds%xml_version)) &
          call wxml_error("Invalid Name in DTD "//name)
      endif
    endif

    if (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
         call wxml_error("Invalid Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Name in DTD "//name)
    endif
    
    if (xf%state_3 == WXML_STATE_3_DURING_DTD) then
      call add_to_buffer(" [", xf%buffer, .false.)
      xf%state_3 = WXML_STATE_3_INSIDE_INTSUBSET
    endif

    if (xf%state_3 /= WXML_STATE_3_INSIDE_INTSUBSET) &
      call wxml_fatal("Cannot define Entity here: "//name)

    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) then
      call close_start_tag(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    endif

    ! Notation only needs checked if not already registered - done above.
#ifdef PGF90
    call register_external_GE(xf%xds, name=name, &
      systemID=system, publicId=public, notation=notation, &
      baseURI=nullURI, wfc=.false.)
#else
    call register_external_GE(xf%xds, name=name, &
      systemID=system, publicId=public, notation=notation, &
      baseURI=null(), wfc=.false.)
#endif
    
    call add_eol(xf)
    
    call add_to_buffer("<!ENTITY "//name, xf%buffer, .false.)
    if (present(public)) then
      call add_to_buffer(" PUBLIC", xf%buffer, .false.)
      call add_to_buffer(" """//public//"""", xf%buffer, .true.)
    else
      call add_to_buffer(" SYSTEM", xf%buffer, .false.)
    endif
    if (scan(system, """")/=0) then
      call add_to_buffer(" '"//system//"'", xf%buffer, .true.)
    else
      call add_to_buffer(" """//system//"""", xf%buffer, .true.)
    endif
    if (present(notation)) then
      call add_to_buffer(" NDATA "//notation, xf%buffer, .false.)
    endif
    call add_to_buffer(">", xf%buffer, .false.)
#endif
  end subroutine xml_AddExternalEntity


  subroutine xml_AddNotation(xf, name, system, public)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: system
    character(len=*), intent(in), optional :: public

#ifndef DUMMYLIB
    type(URI), pointer :: URIref
    call check_xf(xf)

    if (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
         call wxml_error("Invalid Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Name in DTD "//name)
    endif
    if (present(system)) then
      URIref => parseURI(system)
      if (.not.associated(URIref)) call wxml_error("xml_AddNotation: Invalid SYSTEM URI")
      call destroyURI(URIref)
    endif
    if (present(public)) then
      if (.not.checkPublicId(public)) call wxml_error("xml_AddNotation: Invalid PUBLIC ID")
    endif

    if (xf%state_3 == WXML_STATE_3_DURING_DTD) then
      call add_to_buffer(" [", xf%buffer, .false.)
      xf%state_3 = WXML_STATE_3_INSIDE_INTSUBSET
    endif

    if (xf%state_3 /= WXML_STATE_3_INSIDE_INTSUBSET) &
      call wxml_fatal("Cannot define Notation here: "//name)
    
    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) then
      call close_start_tag(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    endif

    if (notation_exists(xf%xds%nList, name)) &
      call wxml_error("Tried to create duplicate notation: "//name)
    
    call add_eol(xf)

    call add_notation(xf%xds%nList, name, system, public)
    call add_to_buffer("<!NOTATION "//name, xf%buffer, .false.)
    if (present(public)) then
      call add_to_buffer(" PUBLIC", xf%buffer, .false.)
      call add_to_buffer(" """//public//"""", xf%buffer, .true.)
    elseif (present(system)) then
      call add_to_buffer(" SYSTEM", xf%buffer, .false.)
    endif
    if (present(system)) then
      if (index(system, """") > 0) then
        call add_to_buffer(" '"//system//"'", xf%buffer, .true.)
      else
        call add_to_buffer(" """//system//"""", xf%buffer, .true.)
      endif
    endif
    call add_to_buffer(">", xf%buffer, .false.)
#endif
  end subroutine xml_AddNotation


  subroutine xml_AddElementToDTD(xf, name, declaration)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: declaration

#ifndef DUMMYLIB
    type(error_stack) :: stack
#ifdef PGF90
    type (element_t), pointer :: nullElement

    nullElement => null()
#endif
    call check_xf(xf)

    if (.not.checkChars(declaration,xf%xds%xml_version)) call wxml_error("xml_AddElementToDTD: Invalid character in declaration")

    if (xf%namespace) then
      if (.not.checkQName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Element Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Element Name in DTD "//name)
    endif
#ifdef PGF90
    call parse_dtd_element(declaration, xf%xds%xml_version, stack, nullElement, .true.)
#else
    call parse_dtd_element(declaration, xf%xds%xml_version, stack, null(), .true.)
#endif
    if (in_error(stack)) call wxml_error(xf, "Invalid ELEMENT declaration")
    
    if (xf%state_3 == WXML_STATE_3_DURING_DTD) then
      call add_to_buffer(" [", xf%buffer, .false.)
      xf%state_3 = WXML_STATE_3_INSIDE_INTSUBSET
    endif

    if (xf%state_3 /= WXML_STATE_3_INSIDE_INTSUBSET) &
      call wxml_fatal("Cannot write to DTD here: xml_AddElementToDTD")

    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) then
      call close_start_tag(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    endif

    call add_eol(xf)
    call add_to_buffer("<!ELEMENT "//name//" "//declaration//">", xf%buffer, .false.)
#endif
  end subroutine xml_AddElementToDTD


  subroutine xml_AddAttlistToDTD(xf, name, declaration)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: declaration

#ifndef DUMMYLIB
    type(error_stack) :: stack
#ifdef PGF90
    type (element_t), pointer :: nullElement

    nullElement => null()
#endif
    call check_xf(xf)

    if (.not.checkChars(declaration,xf%xds%xml_version)) call wxml_error("xml_AddAttListToDTD: Invalid character in declaration")

    if (xf%namespace) then
      if (.not.checkQName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Attribute Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Attribute Name in DTD "//name)
    endif

#ifdef PGF90
    call parse_dtd_attlist(declaration, xf%xds%xml_version, &
      validCheck=.false., namespaces=xf%namespace, stack=stack, &
      elem=nullElement, internal=.true.)
#else
    call parse_dtd_attlist(declaration, xf%xds%xml_version, &
      validCheck=.false., namespaces=xf%namespace, stack=stack, &
      elem=null(), internal=.true.)
#endif

    if (in_error(stack)) call wxml_error(xf, "Invalid ATTLIST declaration")

    if (xf%state_3 == WXML_STATE_3_DURING_DTD) then
      call add_to_buffer(" [", xf%buffer, .false.)
      xf%state_3 = WXML_STATE_3_INSIDE_INTSUBSET
    endif

    if (xf%state_3 /= WXML_STATE_3_INSIDE_INTSUBSET) &
      call wxml_fatal("Cannot write to DTD here: xml_AddAttlistToDTD")

    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) then
      call close_start_tag(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    endif

    call add_eol(xf)
    call add_to_buffer("<!ATTLIST "//name//" "//declaration//">", xf%buffer, .false.)
#endif
  end subroutine xml_AddAttlistToDTD
    

  subroutine xml_AddPEReferenceToDTD(xf, name)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name

#ifndef DUMMYLIB
    call check_xf(xf)

    if (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid PE Name in DTD "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid PE Name in DTD "//name)
    endif

    call wxml_warning(xf, "Adding PEReference to DTD. Cannot guarantee well-formedness")
    if (.not.existing_entity(xf%xds%PEList, name)) then
      if (.not.xf%xds%standalone) then
        call wxml_warning(xf, "Tried to reference possibly unregistered parameter entity in DTD: "//name)
      else
        call wxml_error("Tried to reference unregistered parameter entity in DTD "//name)
      endif
    else
      if (is_unparsed_entity(xf%xds%PEList, name)) &
        call wxml_error("Tried to reference unparsed parameter entity in DTD "//name)
    endif
    
    if (xf%state_3 == WXML_STATE_3_DURING_DTD) then
      call add_to_buffer(" [", xf%buffer, .false.)
      xf%state_3 = WXML_STATE_3_INSIDE_INTSUBSET
    endif

    if (xf%state_3 /= WXML_STATE_3_INSIDE_INTSUBSET) &
      call wxml_fatal("Cannot write to DTD here: xml_AddPEReferenceToDTD")

    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) then
      call close_start_tag(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    endif

    call add_eol(xf)
    call add_to_buffer("%"//name//";", xf%buffer, .false.)

#endif
  end subroutine xml_AddPEReferenceToDTD


  subroutine xml_AddXMLStylesheet(xf, href, type, title, media, charset, alternate)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in) :: href
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: media
    character(len=*), intent(in), optional :: charset
    logical,          intent(in), optional :: alternate

#ifndef DUMMYLIB
    call check_xf(xf)
    ! Don't bother checking name - all pseudoatts get checked anyway.
    
    if (xf%state_1 /= WXML_STATE_1_JUST_OPENED &
         .and. xf%state_1 /= WXML_STATE_1_BEFORE_ROOT) &
      call wxml_error("Cannot add stylesheet here: "//href)

    call close_start_tag(xf)
    
    call xml_AddXMLPI(xf, 'xml-stylesheet', xml=.true.)
    call xml_AddPseudoAttribute(xf, 'href', href)
    call xml_AddPseudoAttribute(xf, 'type', type)
    
    if (present(title)) call xml_AddPseudoAttribute(xf, 'title', title)
    if (present(media)) call xml_AddPseudoAttribute(xf, 'media', media)
    if (present(charset)) call xml_AddPseudoAttribute(xf, 'charset', charset)
    if (present(alternate)) then
      if (alternate) then
        call xml_AddPseudoAttribute(xf, 'alternate', 'yes')
      else
        call xml_AddPseudoAttribute(xf, 'alternate', 'no')
      endif
    endif
    if (xf%state_1 == WXML_STATE_1_JUST_OPENED) &
         xf%state_1 = WXML_STATE_1_BEFORE_ROOT 
    xf%state_2 = WXML_STATE_2_INSIDE_PI
#endif
  end subroutine xml_AddXMLStylesheet
  

  subroutine xml_AddXMLPI(xf, name, data, xml, ws_significant)
    type(xmlf_t), intent(inout)            :: xf
    character(len=*), intent(in)           :: name
    character(len=*), intent(in), optional :: data
    logical, intent(in), optional :: xml
    logical, intent(in), optional :: ws_significant

    logical :: xml_
#ifndef DUMMYLIB
    call check_xf(xf)

    if (present(xml)) then
      xml_ = xml
    else
      xml_ = .false.
    endif

    if (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid PI target "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid PI target "//name)
    endif
    if (.not.xml_) then
      if (len(name)==3.and.(toLower(name)=="xml")) &
        call wxml_error("Invalid PI target "//name)
    endif

    if (present(data)) then
      if (.not.checkChars(data,xf%xds%xml_version)) &
        call wxml_error("xml_AddXMLPI: Invalid character in data")
    endif

    select case (xf%state_1)
    case (WXML_STATE_1_JUST_OPENED) 
      xf%state_1 = WXML_STATE_1_BEFORE_ROOT
    case (WXML_STATE_1_DURING_ROOT)
      call close_start_tag(xf)
      if (xf%pretty_print) call add_eol(xf)
    case default
      call close_start_tag(xf)
      call add_eol(xf)
    end select
    call add_to_buffer("<?" // name, xf%buffer, .false.)
    if (present(data)) then
      if (len(data)>0) then
        if (index(data, '?>') > 0) &
          call wxml_error(xf, "Tried to output invalid PI data "//data)
        call add_to_buffer(' ', xf%buffer, .false.)
        call add_to_buffer(data//'?>', xf%buffer, ws_significant)
        ! state_2 is now OUTSIDE_TAG from close_start_tag
      else
        xf%state_2 = WXML_STATE_2_INSIDE_PI
        call reset_dict(xf%dict)
      endif
    else
      xf%state_2 = WXML_STATE_2_INSIDE_PI
      call reset_dict(xf%dict)
    endif
#endif
  end subroutine xml_AddXMLPI


  subroutine xml_AddComment(xf, comment, ws_significant)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in)  :: comment
    logical, intent(in), optional :: ws_significant

#ifndef DUMMYLIB
    call check_xf(xf)
    if (.not.checkChars(comment,xf%xds%xml_version)) call wxml_error("xml_AddComment: Invalid character in comment")
    
    select case (xf%state_1)
    case (WXML_STATE_1_JUST_OPENED) 
      xf%state_1 = WXML_STATE_1_BEFORE_ROOT
    case (WXML_STATE_1_DURING_ROOT)
      call close_start_tag(xf)
      if (xf%pretty_print.and.xf%state_2 == WXML_STATE_2_OUTSIDE_TAG) call add_eol(xf)
    case default
      call close_start_tag(xf)
      call add_eol(xf)
    end select

    if (index(comment,'--') > 0 .or. comment(len(comment):) == '-') &
         call wxml_error("Tried to output invalid comment "//comment)

    call add_to_buffer("<!--", xf%buffer, .false.)
    call add_to_buffer(comment, xf%buffer, ws_significant)
    call add_to_buffer("-->", xf%buffer, .false.)
#endif
  end subroutine xml_AddComment


  subroutine xml_NewElement(xf, name)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in)  :: name

#ifndef DUMMYLIB
    call check_xf(xf)
    
    if (xf%namespace) then
      if (.not.checkQName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Element Name "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Element Name "//name)
    endif
    
    select case (xf%state_1)
    case (WXML_STATE_1_JUST_OPENED, WXML_STATE_1_BEFORE_ROOT)
      if (xf%xds%valid) then
        if (size(xf%name)==0) then
          call wxml_error(xf, "No DTD specified for document")
        elseif (str_vs(xf%name) /= name) then
          call wxml_error(xf, "Root element name does not match DTD")
        endif
      endif
      call close_start_tag(xf)
      if (xf%state_3 /= WXML_STATE_3_BEFORE_DTD) then
        select case (xf%state_3)
        case (WXML_STATE_3_DURING_DTD)
          call add_to_buffer('>', xf%buffer, .false.)
          xf%state_3 = WXML_STATE_3_AFTER_DTD
        case (WXML_STATE_3_INSIDE_INTSUBSET)
          xf%state_3 = WXML_STATE_3_AFTER_DTD
          call add_eol(xf)
          call add_to_buffer(']>', xf%buffer, .false.)
        end select
      endif
      call add_eol(xf)
    case (WXML_STATE_1_DURING_ROOT)
      call close_start_tag(xf)
      if (xf%pretty_print) call add_eol(xf)
    case (WXML_STATE_1_AFTER_ROOT)
      call wxml_error(xf, "Two root elements: "//name)
    end select

    if (xf%namespace) then
      if (len(prefixOfQName(name)) > 0) then
        if (.not.isPrefixInForce(xf%nsDict, prefixOfQName(name))) &
          call wxml_error(xf, "Namespace prefix not registered: "//prefixOfQName(name))
      endif
    endif
    
    call push_elstack(xf%stack, name)
    call add_to_buffer("<"//name, xf%buffer, .false.)
    xf%state_2 = WXML_STATE_2_INSIDE_ELEMENT
    call reset_dict(xf%dict)
    xf%indent = xf%indent + indent_inc
    xf%state_1 = WXML_STATE_1_DURING_ROOT
#endif
  end subroutine xml_NewElement
  

  subroutine xml_AddCharacters_ch(xf, chars, parsed, ws_significant)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in)  :: chars
    logical, intent(in), optional :: parsed
    logical, intent(in), optional :: ws_significant

#ifndef DUMMYLIB
    logical :: pc

    call check_xf(xf)
    if (.not.checkChars(chars, xf%xds%xml_version)) call wxml_error("xml_AddCharacters: Invalid character in chars")
    
    if (xf%state_1 /= WXML_STATE_1_DURING_ROOT) &
         call wxml_fatal("Tried to add text section in wrong place: "//chars)
    
    if (present(parsed)) then
      pc = parsed
    else
      pc = .true.
    endif
    
    call close_start_tag(xf)

    if (pc) then
      call add_to_buffer(escape_string(chars, xf%xds%xml_version), xf%buffer, ws_significant)
    else
      ! FIXME what if we try and output two separate character events?
      ! need to keep track of this ...
      if (index(chars,']]>') > 0) &
           call wxml_fatal("Tried to output invalid CDATA: "//chars)
      call add_to_buffer("<![CDATA["//chars//"]]>", xf%buffer, ws_significant)
    endif
    
    xf%state_2 = WXML_STATE_2_IN_CHARDATA
#endif
  end subroutine xml_AddCharacters_Ch


  subroutine xml_AddNewline(xf)
    type(xmlf_t), intent(inout) :: xf

#ifndef DUMMYLIB
    call xml_AddCharacters(xf, "") ! To ensure we are in a text section
    call add_eol(xf)
#endif
  end subroutine xml_AddNewline

  
  subroutine xml_AddEntityReference(xf, name)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name

#ifndef DUMMYLIB
    call check_xf(xf)

    if (likeCharacterEntityReference(name)) then
      if (.not.checkCharacterEntityReference(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Character Entity Reference "//name)
    elseif (xf%namespace) then
      if (.not.checkNCName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Entity Name "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Entity Name "//name)
    endif

    call close_start_tag(xf)

    if (xf%state_2 /= WXML_STATE_2_OUTSIDE_TAG .and. &
      xf%state_2 /= WXML_STATE_2_IN_CHARDATA)         &
      call wxml_fatal("Tried to add entity reference in wrong place: "//name)

    if (.not.checkCharacterEntityReference(name, xf%xds%xml_version)) then
      !it's not just a unicode entity
      call wxml_warning(xf, "Entity reference added - document may not be well-formed")
      if (.not.existing_entity(xf%xds%entityList, name)) then
        if (xf%xds%standalone) then
          call wxml_error("Tried to reference unregistered entity")
        else
          call wxml_warning(xf, "Tried to reference unregistered entity")
        endif
      else
        if (is_unparsed_entity(xf%xds%entityList, name)) &
          call wxml_error("Tried to reference unparsed entity")
      endif
    endif

    call add_to_buffer('&'//name//';', xf%buffer, .false.)
    xf%state_2 = WXML_STATE_2_IN_CHARDATA
#endif
  end subroutine xml_AddEntityReference


  subroutine xml_AddAttribute_Ch(xf, name, value, escape, type, ws_significant)
    type(xmlf_t), intent(inout)             :: xf
    character(len=*), intent(in)            :: name
    character(len=*), intent(in)            :: value
    logical, intent(in), optional           :: escape
    character(len=*), intent(in), optional  :: type
    logical, intent(in), optional           :: ws_significant

#ifndef DUMMYLIB
    logical :: esc
    character, pointer :: type_(:)

    if (present(type)) then
      if (type/='CDATA'.and.type/='ID'.and.type/='IDREF'.and.type/='IDREFS'.and.type/='NMTOKEN'.and.type/='NMTOKENS' &
        .and.type/='ENTITY'.and.type/='ENTITIES'.and.type/='NOTATION') then
        call wxml_fatal("Invalid type in xml_AddAttribute: "//type)
      endif
      type_ => vs_str_alloc(type)
    else
      ! We assume CDATA, but need to worry about whether the caller cares about whitespace ...
      if (present(ws_significant)) then
        if (ws_significant) then
          type_ => vs_str_alloc('CDATA')
        else
          type_ => vs_str_alloc('CDANO') ! CDAta, whitespace Not significant
        endif
      else
        type_ => vs_str_alloc('CDAMB')   ! CDAta, whitespace MayBe significant
      endif
    endif

    call check_xf(xf)

    if (.not.checkChars(value, xf%xds%xml_version)) call wxml_error("xml_AddAttribute: Invalid character in value")

    if (xf%namespace) then
      if (.not.checkQName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Attribute Name "//name)
    else
      if (.not.checkName(name, xf%xds%xml_version)) &
        call wxml_error("Invalid Attribute Name "//name)
    endif

    if (present(escape)) then
      esc = escape
    else
      esc = .true.
    endif

    if (name=="xml:space") then
      ! The value can only be "default" or "preserve", by 2.10
      if (.not.esc) then
        if (value/="default".and.value/="preserve") & 
          call wxml_fatal("Invalid value for xml:space attrbute")
      endif
    endif

    ! FIXME when escape is false we should still do full verification
    ! where possible.
    ! Currently - minimal check: only extra allowed is character entity references.
    ! We check they exist, and are not unparsed.
    ! Ideally we would fully expand all entity references (at least for
    ! a standalone document where we can) and then
    ! match the resultant production against [XML]-3.3.1. This is
    ! initially too much work though, so we just check simple
    ! syntactic constraint.

    if (.not.esc) then
      if (.not.checkAttValue(value, xf%xds%xml_version)) &
        call wxml_error(xf, "Invalid attribute value: "//value)
      if (index(value, '&') > 0) then
        ! There are entity references
        ! They should exist (unless we are not standalone) and they must not be unparsed.
        if (.not.checkExistingRefsInAttValue()) then
          if (xf%xds%standalone) then
            call wxml_error(xf, "outputting unknown entity. Cannot guarantee validity.")
          else
            call wxml_warning(xf, "Warning: outputting unknown entity. Cannot guarantee validity.")
          endif
        endif
        if (.not.checkParsedRefsInAttValue()) &
          call wxml_error(xf, "Warning: outputting unknown entity. Cannot guarantee validity.")
      endif
    endif

    if (xf%state_2 /= WXML_STATE_2_INSIDE_ELEMENT) &
         call wxml_error(xf, "attributes outside element content: "//name)

    if (hasKey(xf%dict,name)) then
      call wxml_error(xf, "duplicate att name: "//name)
    elseif (xf%namespace) then
      if (hasKey(xf%dict, &
        getnamespaceURI(xf%nsDict,prefixOfQname(name)), localpartofQname(name))) then
        call wxml_error(xf, "duplicate att after namespace processing: "//name)
      endif
    endif

    if (xf%namespace) then
      if (len(prefixOfQName(name))>0) then
        if (prefixOfQName(name)/="xml".and.prefixOfQName(name)/="xmlns") then
          if (.not.isPrefixInForce(xf%nsDict, prefixOfQName(name))) &
            call wxml_error(xf, "namespace prefix not registered: "//prefixOfQName(name))
        endif
        if (esc) then
          call add_item_to_dict(xf%dict, localpartofQname(name), escape_string(value, xf%xds%xml_version), prefixOfQName(name), &
            getnamespaceURI(xf%nsDict,prefixOfQname(name)), type=str_vs(type_))
        else
          call add_item_to_dict(xf%dict, localpartofQname(name), value, prefixOfQName(name), &
            getnamespaceURI(xf%nsDict,prefixOfQName(name)), type=str_vs(type_))
        endif
      else
        if (esc) then
          call add_item_to_dict(xf%dict, name, escape_string(value, xf%xds%xml_version), type=str_vs(type_))
        else
          call add_item_to_dict(xf%dict, name, value, type=str_vs(type_))
        endif
      endif
    else
      if (esc) then
        call add_item_to_dict(xf%dict, name, escape_string(value, xf%xds%xml_version), type=str_vs(type_))
      else
        call add_item_to_dict(xf%dict, name, value, type=str_vs(type_))
      endif
    endif

    !FIXME need to deallocate this when we move to better error handling
    deallocate(type_)

  contains
    function checkExistingRefsInAttValue() result(p)
      logical :: p
      
      integer :: i1, i2
      
      ! Here we assume we have syntactic well-formedness as
      ! checked by checkAttValue.
      ! We also assume we do not have simply one entity as
      ! the contents - that is checked by checkAttValueEntity
      
      p = .false.
      i1 = index(value, '&')
      i2 = 0 
      do while (i1 > 0)
        i1 = i2 + i1
        i2 = index(value(i1+1:),';')
        if (i2 == 0) return
        i2 = i1 + i2
        if (.not.existing_entity(xf%xds%entityList, value(i1+1:i2-1)) .and. &
          .not.checkCharacterEntityReference(value(i1+1:i2-1), xf%xds%xml_version)) &
          return
        i1 = index(value(i2+1:), '&')
      enddo
      p = .true.
      
    end function checkExistingRefsInAttValue
    
    function checkParsedRefsInAttValue() result(p)
      logical :: p
      
      integer :: i1, i2
      
      ! Here we assume we have syntactic well-formedness as
      ! checked by checkAttValue.
      
      p = .false.
      i1 = index(value, '&')
      i2 = 0
      do while (i1 > 0)
        i1 = i1 + i2
        i2 = index(value(i1+1:),';')
        if (i2 == 0) return
        i2  = i1 + i2
        if (is_unparsed_entity(xf%xds%entityList, value(i1+1:i2-1))) &
          return
        i1 = index(value(i2+1:), '&')
      enddo
      p = .true.

    end function checkParsedRefsInAttValue
#endif
  end subroutine xml_AddAttribute_Ch


  subroutine xml_AddPseudoAttribute_Ch(xf, name, value, escape, ws_significant)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in)  :: name
    character(len=*), intent(in)  :: value
    logical, intent(in), optional :: escape
    logical, intent(in), optional :: ws_significant

#ifndef DUMMYLIB
    logical :: esc
    character(len=5) :: type

    call check_xf(xf)
    if (.not.checkChars(name, xf%xds%xml_version)) call wxml_error("xml_AddPseudoAttribute: Invalid character in name")
    if (.not.checkChars(value, xf%xds%xml_version)) call wxml_error("xml_AddPseudoAttribute: Invalid character in value")
    
    if (present(escape)) then
      esc = escape
    else
      esc = .true.
    endif
    if (present(ws_significant)) then
      if (ws_significant) then
        type='CDATA'
      else
        type='CDANO' ! CDAta, whitespace Not significant
      endif
    else
      type='CDAMB'   ! CDAta, whitespace MayBe significant
    endif

    if (index(value, '?>') > 0) &
        call wxml_error(xf, "Invalid pseudo-attribute value: "//value)
    if (.not.esc) then
      if (.not.checkPseudoAttValue(value, xf%xds%xml_version)) &
        call wxml_error(xf, "Invalid pseudo-attribute value: "//value)
    endif

    if (xf%state_2 /= WXML_STATE_2_INSIDE_PI) &
         call wxml_error("PI pseudo-attribute outside PI: "//name)

    ! This is mostly ad-hoc, pseudo-attribute names are not defined anywhere.
    if (.not.checkName(name, xf%xds%xml_version)) &
         call wxml_error("Invalid pseudo-attribute name: "//name)

    if (hasKey(xf%dict,name)) &
         call wxml_error(xf, "duplicate pseudo-attribute name: "//name)

    if (index(value, '?>') > 0) &
         call wxml_error(xf, "Invalid pseudo-attribute data: "//value)
    
    if (esc) then
      call add_item_to_dict(xf%dict, name, escape_string(value, xf%xds%xml_version), type=type)
    else
      call add_item_to_dict(xf%dict, name, value, type=type)
    endif
#endif
  end subroutine xml_AddPseudoAttribute_Ch


  subroutine xml_EndElement(xf, name)
    type(xmlf_t), intent(inout)             :: xf
    character(len=*), intent(in)            :: name

    character :: dummy
#ifndef DUMMYLIB
    call check_xf(xf)
    ! No point in doing checkChars, name is compared to stack anyway.

    if (len(xf%stack) == 0) &
      call wxml_fatal(xf,'Trying to close '//name//' but no tags are open.')

    if (get_top_elstack(xf%stack) /= name) &
      call wxml_fatal(xf, 'Trying to close '//name//' but '//get_top_elstack(xf%stack)// &
      ' is open. Either you have failed to open '//name//&
      ' or you have failed to close '//get_top_elstack(xf%stack)//'.') 
    xf%indent = xf%indent - indent_inc

    if (xf%state_2==WXML_STATE_2_INSIDE_ELEMENT) then
      if (xf%namespace) call checkNamespacesWriting(xf%dict, xf%nsDict, len(xf%stack))
      if (getLength(xf%dict) > 0) call write_attributes(xf)
      if (xf%minimize_overrun) call add_eol(xf)
    endif
    if (xf%state_2==WXML_STATE_2_INSIDE_ELEMENT.and..not.xf%canonical) then
      call add_to_buffer("/>", xf%buffer, .false.)
    else
      if (xf%state_2==WXML_STATE_2_INSIDE_ELEMENT) &
        call add_to_buffer('>', xf%buffer, .false.)
      if (xf%state_2==WXML_STATE_2_INSIDE_PI) &
        call close_start_tag(xf)
      if (xf%state_2==WXML_STATE_2_OUTSIDE_TAG.and.xf%pretty_print) &
        call add_eol(xf)
! XLF does a weird thing here, and if pop_elstack is called as an 
! argument to the add_to_buffer, it gets called twice. So we have to separate
! out get_top_... from pop_...
      call add_to_buffer("</" //get_top_elstack(xf%stack), xf%buffer, .false.)
      if (xf%minimize_overrun) call add_eol(xf)
      call add_to_buffer(">", xf%buffer, .false.)
    endif
    dummy = pop_elstack(xf%stack)

    if (xf%namespace) call checkEndNamespaces(xf%nsDict, len(xf%stack)+1)
    if (is_empty(xf%stack)) then
      xf%state_1 = WXML_STATE_1_AFTER_ROOT
    endif
    xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
#endif
  end subroutine xml_EndElement


  subroutine xml_DeclareNamespace(xf, nsURI, prefix, xml)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in) :: nsURI
    character(len=*), intent(in), optional :: prefix
    logical, intent(in), optional :: xml

#ifndef DUMMYLIB
    call check_xf(xf)
    if (.not.xf%namespace) call wxml_error("Cannot declare a namespace in a non-namespaced document")

    !if (.not.checkNCName(nsURI, xf%xds%xml_version)) call wxml_error("xml_DeclareNamespace: Invalid nsURI")
    if (present(prefix)) then
      if (.not.checkNCName(prefix, xf%xds%xml_version)) call wxml_error("xml_DeclareNamespace: Invalid prefix")
    endif

    if (xf%state_1 == WXML_STATE_1_AFTER_ROOT) &
      call wxml_error(xf, "adding namespace outside element content")

    if (len(nsURI) == 0) then
      if (present(prefix).and.xf%xds%xml_version==XML1_0) &
        call wxml_error(xf, "prefixed namespace with empty URI forbidden in XML 1.0")
    endif

    if (present(prefix)) then
      call addPrefixedNS(xf%nsDict, prefix, nsURI, len(xf%stack)+1, xf%xds, xml)
    else
      call addDefaultNS(xf%nsDict, nsURI, len(xf%stack)+1)
    endif
#endif
  end subroutine xml_DeclareNamespace


  subroutine xml_UndeclareNamespace(xf, prefix)
    type(xmlf_t), intent(inout)   :: xf
    character(len=*), intent(in), optional :: prefix

#ifndef DUMMYLIB
    call check_xf(xf)
    !No need to checkChars, prefix is checked against stack
    if (.not.xf%namespace) call wxml_error("Cannot declare a namespace in a non-namespaced document")

    if (present(prefix).and.xf%xds%xml_version==XML1_0) &
      call wxml_error("cannot undeclare prefixed namespaces in XML 1.0")
    
    if (xf%state_1 == WXML_STATE_1_AFTER_ROOT) &
      call wxml_error(xf, "Undeclaring namespace outside element content")
    
    if (present(prefix)) then
      call addPrefixedNS(xf%nsDict, prefix, "", len(xf%stack)+1, xf%xds)
    else
      call addDefaultNS(xf%nsDict, "", len(xf%stack)+1)
    endif
#endif
  end subroutine xml_UndeclareNamespace


  subroutine xml_Close(xf, empty)
    type(xmlf_t), intent(inout)   :: xf
    logical, optional :: empty

#ifndef DUMMYLIB
    logical :: empty_

    if (present(empty)) then
      empty_ = empty
    else
      empty_ = .false.
    endif

    if (xf%lun == -1) &
      call wxml_fatal('Tried to close XML file which is not open')

    if (xf%state_2 == WXML_STATE_2_INSIDE_PI) &
      call close_start_tag(xf)

    if (xf%state_3 /= WXML_STATE_3_BEFORE_DTD &
      .and. xf%state_3 /= WXML_STATE_3_AFTER_DTD) then
      select case (xf%state_3)
      case (WXML_STATE_3_DURING_DTD)
        call add_to_buffer('>', xf%buffer, .false.)
      case (WXML_STATE_3_INSIDE_INTSUBSET)
        call add_eol(xf)
        call add_to_buffer(']>', xf%buffer, .false.)
      end select
      xf%state_3 = WXML_STATE_3_AFTER_DTD
    endif
    
    do while (xf%state_1 == WXML_STATE_1_DURING_ROOT)
      call xml_EndElement(xf, get_top_elstack(xf%stack))
    enddo

    if (xf%state_1 /= WXML_STATE_1_AFTER_ROOT) then
      if (empty_) then
        call wxml_warning(xf, 'Invalid XML document produced: No root element')
      else
        call wxml_error(xf, 'Invalid XML document produced: No root element')
      endif
    endif
    
    call dump_buffer(xf%buffer)
    close(unit=xf%lun)
    xf%lun = -1

    call destroy_dict(xf%dict)
    call destroy_elstack(xf%stack)
    
    if (xf%namespace) &
      call destroyNamespaceDictionary(xf%nsDict)
    call destroy_xml_doc_state(xf%xds)
    
    deallocate(xf%name)
#endif
  end subroutine xml_Close

  subroutine xmlf_SetPretty_print(xf, new_value)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in)         :: new_value
#ifndef DUMMYLIB
    xf%pretty_print = new_value
#endif
  end subroutine xmlf_SetPretty_print

  pure function xmlf_GetPretty_print(xf) result(value)
    logical :: value
    type(xmlf_t), intent(in) :: xf
#ifdef DUMMYLIB
    value = .false.
#else
    value = xf%pretty_print
#endif
  end function xmlf_GetPretty_print

  pure function xmlf_name(xf) result(fn)
    type (xmlf_t), intent(in) :: xf
#ifdef DUMMYLIB
    character(len=1) :: fn
    fn = " "
#else
    character(len=size(xf%xds%documentURI)) :: fn
    fn = str_vs(xf%xds%documentURI)
#endif
  end function xmlf_name

#ifndef DUMMYLIB
  pure function xmlf_opentag_len(xf) result(n)
    type (xmlf_t), intent(in) :: xf
    integer :: n
    
    if (xf%lun == -1) then
      n = 0
    elseif (is_empty(xf%stack)) then
      n = 0
    else
      n = len(get_top_elstack(xf%stack))
    endif
  end function xmlf_opentag_len
#endif

  function xmlf_opentag(xf) result(fn)
    type (xmlf_t), intent(in) :: xf
#ifdef DUMMYLIB
    character(len=1) :: fn
    fn = " "
#else
    character(len=xmlf_opentag_len(xf)) :: fn
    
    if (xf%lun == -1) then
      fn = ''
    elseif (is_empty(xf%stack)) then
      fn = ''
    else
      fn = get_top_elstack(xf%stack)
    endif
#endif
  end function xmlf_opentag

#ifndef DUMMYLIB

  subroutine check_xf(xf)
    type(xmlf_t), intent(inout)   :: xf
    if (xf%lun == -1) &
      call wxml_fatal("Tried to manipulate an XML File which is not open")

  end subroutine check_xf


  subroutine add_eol(xf)
    type(xmlf_t), intent(inout)   :: xf
    
    integer :: indent_level
    
    ! In case we still have a zero-length stack, we must make
    ! sure indent_level is not less than zero.
    if (xf%state_3 == WXML_STATE_3_INSIDE_INTSUBSET) then
      indent_level = indent_inc
    else
      indent_level = xf%indent
    endif
    
    !We must flush here (rather than just adding an eol character)
    !since we don't know what the eol character is on this system.
    !Flushing with a linefeed will get it automatically, though.
    call dump_buffer(xf%buffer, lf=.true.)
    call reset_buffer(xf%buffer, xf%lun, xf%xds%xml_version)
    
    if (xf%pretty_print) &
      call add_to_buffer(repeat(' ',indent_level),xf%buffer, .false.)

  end subroutine add_eol


  subroutine close_start_tag(xf)
    type(xmlf_t), intent(inout)   :: xf

    select case (xf%state_2)
    case (WXML_STATE_2_INSIDE_ELEMENT)
      if (xf%namespace) call checkNamespacesWriting(xf%dict, xf%nsDict, len(xf%stack))
      if (getLength(xf%dict) > 0) call write_attributes(xf)
      if (xf%minimize_overrun) call add_eol(xf)
      call add_to_buffer('>', xf%buffer, .false.)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    case (WXML_STATE_2_INSIDE_PI)
      if (getLength(xf%dict) > 0) call write_attributes(xf)
      call add_to_buffer('?>', xf%buffer, .false.)
      if (xf%pretty_print.and.xf%state_3/=WXML_STATE_3_INSIDE_INTSUBSET) call add_eol(xf)
      xf%state_2 = WXML_STATE_2_OUTSIDE_TAG
    case (WXML_STATE_2_IN_CHARDATA)
      continue
    case (WXML_STATE_2_OUTSIDE_TAG)
      continue
    end select

  end subroutine close_start_tag


  subroutine write_attributes(xf)
    type(xmlf_t), intent(inout)   :: xf

    integer  :: i, j, size

    if (xf%state_2 /= WXML_STATE_2_INSIDE_PI .and. &
      xf%state_2 /= WXML_STATE_2_INSIDE_ELEMENT) &
      call wxml_fatal("Internal library error")

    if (xf%canonical) call sortAttrs(xf%dict)
    
    do i = 1, getLength(xf%dict)
      size = len(get_key(xf%dict, i)) + len(get_value(xf%dict, i)) + 4
      if (xf%minimize_overrun.and.(len(xf%buffer) + size) > COLUMNS) then
        call add_eol(xf)
      else
        call add_to_buffer(" ", xf%buffer, .false.)
      endif
      call add_to_buffer(get_key(xf%dict, i), xf%buffer, .false.)
      call add_to_buffer("=", xf%buffer, .false.)
      call add_to_buffer('"',xf%buffer, .false.)
      j = getWhiteSpaceHandling(xf%dict, i)
      if (j==0) then
        call add_to_buffer(get_value(xf%dict, i), xf%buffer, .true.)
      elseif (j==1) then
        call add_to_buffer(get_value(xf%dict, i), xf%buffer)
      else
        call add_to_buffer(get_value(xf%dict, i), xf%buffer, .false.)
      endif
      call add_to_buffer('"', xf%buffer, .false.)
    enddo

  end subroutine write_attributes
  
  subroutine wxml_warning_xf(xf, msg)
    ! Emit warning, but carry on.
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: msg
    
    if (FoX_get_fatal_warnings()) then
        write(6,'(a)') 'FoX warning made fatal'
        call wxml_fatal_xf(xf, msg)
    endif

    if (xf%xds%warning) then
      write(6,'(a)') 'WARNING(wxml) in writing to file ', xmlf_name(xf)
      write(6,'(a)')  msg
    endif
    
  end subroutine wxml_warning_xf
  
  
  subroutine wxml_error_xf(xf, msg)
    ! Emit error message, clean up file and stop.
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: msg
   
    if (FoX_get_fatal_errors()) then
        write(6,'(a)') 'FoX error made fatal'
        call wxml_fatal_xf(xf, msg)
    endif
 
    write(6,'(a)') 'ERROR(wxml) in writing to file ', xmlf_name(xf)
    write(6,'(a)')  msg
    
    !call xml_Close(xf)
    stop
    
  end subroutine wxml_error_xf
  
  
  subroutine wxml_fatal_xf(xf, msg)
    !Emit error message and abort with coredump. Does not try to
    !close file, so should be used from anything xml_Close might
    !itself call (to avoid infinite recursion!)
    
    type(xmlf_t), intent(in) :: xf
    character(len=*), intent(in) :: msg
    
    write(6,'(a)') 'ERROR(wxml) in writing to file ', xmlf_name(xf)
    write(6,'(a)')  msg
    
    call pxfabort()
    stop
    
  end subroutine wxml_fatal_xf

#endif

end module m_wxml_core
