module m_common_namespaces

#ifndef DUMMYLIB
  use fox_m_fsys_array_str, only: str_vs, vs_str, vs_str_alloc

  use fox_m_utils_uri, only: URI, parseURI, destroyURI, hasScheme
  use m_common_attrs, only: dictionary_t, get_key, get_value, remove_key, getLength, hasKey
  use m_common_attrs, only: set_nsURI, set_localName, get_prefix, add_item_to_dict
  use m_common_charset, only: XML1_0, XML1_1
  use m_common_error, only: FoX_error, FoX_warning, error_stack, add_error, in_error
  use m_common_namecheck, only: checkNCName
  use m_common_struct, only: xml_doc_state

  implicit none
  private

  character(len=*), parameter :: invalidNS = '::INVALID::'
  ! an invalid URI name to indicate a namespace error.

  type URIMapping
    character, dimension(:), pointer :: URI
    integer :: ix ! link back to node depth
  end type URIMapping
  !This is a single URI, and the node depth under which
  !its namespace applies.

  type prefixMapping
    character, dimension(:), pointer :: prefix
    type(URIMapping), dimension(:), pointer :: urilist
  end type prefixMapping
  !This is the mapping for a single prefix; with the
  !list of namespaces which are in force at various
  !depths

  type namespaceDictionary
    private
    type(URIMapping), dimension(:), pointer :: defaults
    type(prefixMapping), dimension(:), pointer :: prefixes
  end type namespaceDictionary
  !This is the full namespace dictionary; defaults is
  !the list of default namespaces in force; prefix a
  !list of all prefixes in force.

  public :: invalidNS

  public :: initNamespaceDictionary
  public :: destroyNamespaceDictionary
  public :: namespaceDictionary
  public :: checkNamespaces
  public :: checkNamespacesWriting
  public :: checkEndNamespaces
  public :: getnamespaceURI
  interface getnamespaceURI
     module procedure getURIofDefaultNS, getURIofPrefixedNS
  end interface
  public :: isPrefixInForce
  public :: isDefaultNSInForce
  public :: getNumberOfPrefixes
  public :: getPrefixByIndex

  public :: dumpnsdict !FIXME

  public :: addDefaultNS
  public :: removeDefaultNS
  public :: addPrefixedNS
  public :: removePrefixedNS

contains


  subroutine initNamespaceDictionary(nsDict)
    type(namespaceDictionary), intent(inout) :: nsDict

    !We need to properly initialize 0th elements
    !(which are never used) in order to provide
    !sensible behaviour when trying to manipulate
    !an empty dictionary.

    allocate(nsDict%defaults(0:0))
    allocate(nsDict%defaults(0)%URI(0))
    !The 0th element of the defaults NS is the empty namespace
    nsDict%defaults(0)%ix = -1
    
    allocate(nsDict%prefixes(0:0))
    allocate(nsDict%prefixes(0)%prefix(0))
    allocate(nsDict%prefixes(0)%urilist(0:0))
    allocate(nsDict%prefixes(0)%urilist(0)%URI(len(invalidNS)))
    nsDict%prefixes(0)%urilist(0)%URI = vs_str(invalidNS)
    nsDict%prefixes(0)%urilist(0)%ix = -1

  end subroutine initNamespaceDictionary


  subroutine destroyNamespaceDictionary(nsDict)
    type(namespaceDictionary), intent(inout) :: nsDict

    integer :: i, j

    do i = 0, ubound(nsDict%defaults,1)
       deallocate(nsDict%defaults(i)%URI)
    enddo
    deallocate(nsDict%defaults)
    do i = 0, ubound(nsDict%prefixes,1)
       do j = 0, ubound(nsDict%prefixes(i)%urilist,1)
          deallocate(nsDict%prefixes(i)%urilist(j)%URI)
       enddo
       deallocate(nsDict%prefixes(i)%prefix)
       deallocate(nsDict%prefixes(i)%urilist)
    enddo
    deallocate(nsDict%prefixes)
  end subroutine destroyNamespaceDictionary


  subroutine copyURIMapping(urilist1, urilist2, l_m)
    type(URIMapping), dimension(0:), intent(inout) :: urilist1
    type(URIMapping), dimension(0:), intent(inout) :: urilist2
    integer, intent(in):: l_m
    integer :: i

    if (ubound(urilist1,1) < l_m .or. ubound(urilist2,1) < l_m) then
       call FoX_error('Internal error in m_sax_namespaces:copyURIMapping')
    endif
    ! Now copy all defaults across (or rather - add pointers to them)
    do i = 0, l_m
       urilist2(i)%ix = urilist1(i)%ix
       urilist2(i)%URI => urilist1(i)%URI
    enddo

  end subroutine copyURIMapping
    
  
  subroutine addDefaultNS(nsDict, uri, ix, es)
    type(namespaceDictionary), intent(inout) :: nsDict
    character(len=*), intent(in) :: uri
    integer, intent(in) :: ix
    type(error_stack), intent(inout), optional :: es

    type(URIMapping), dimension(:), allocatable :: tempMap
    integer :: l_m, l_s

    if (uri=="http://www.w3.org/XML/1998/namespace") then
      if (present(es)) then
        call add_error(es, "Attempt to assign incorrect URI to prefix 'xml'")
      else
        call FoX_error("Attempt to assign incorrect URI to prefix 'xml'")
      endif
    elseif (uri=="http://www.w3.org/2000/xmlns/") then
      if (present(es)) then
        call add_error(es, "Attempt to assign prefix to xmlns namespace")
      else
        call FoX_error("Attempt to assign prefix to xmlns namespace")
      endif
    endif

    ! FIXME check URI is valid ...

    l_m = ubound(nsDict%defaults,1)
    allocate(tempMap(0:l_m))
    ! Now copy all defaults across ...
    call copyURIMapping(nsDict%defaults, tempMap, l_m)
    deallocate(nsDict%defaults)
    l_m = l_m + 1
    allocate(nsDict%defaults(0:l_m))
    !Now copy everything back ...
    call copyURIMapping(tempMap, nsDict%defaults, l_m-1)
    deallocate(tempMap)
    ! And finally, add the new default NS
    nsDict%defaults(l_m)%ix = ix
    l_s = len(uri)
    allocate(nsDict%defaults(l_m)%URI(l_s))
    nsDict%defaults(l_m)%URI = vs_str(uri)

  end subroutine addDefaultNS
  

  subroutine addPrefixedURI(nsPrefix, uri, ix)
    type(PrefixMapping), intent(inout) :: nsPrefix
    character, dimension(:), intent(in) :: uri
    integer, intent(in) :: ix

    type(URIMapping), dimension(:), allocatable :: tempMap
    integer :: l_m, l_s

    l_m = ubound(nsPrefix%urilist,1)
    allocate(tempMap(0:l_m))
    ! Now copy all across ...
    call copyURIMapping(nsPrefix%urilist, tempMap, l_m)
    deallocate(nsPrefix%urilist)
    l_m = l_m + 1
    allocate(nsPrefix%urilist(0:l_m))
    !Now copy everything back ...
    call copyURIMapping(tempMap, nsPrefix%urilist, l_m-1)
    deallocate(tempMap)
    ! And finally, add the new default NS
    nsPrefix%urilist(l_m)%ix = ix
    l_s = size(uri)
    allocate(nsPrefix%urilist(l_m)%URI(l_s))
    nsPrefix%urilist(l_m)%URI = uri

  end subroutine addPrefixedURI
   
  subroutine removeDefaultNS(nsDict)
    type(namespaceDictionary), intent(inout) :: nsDict

    type(URIMapping), dimension(:), allocatable :: tempMap
    integer :: l_m

    l_m = ubound(nsDict%defaults,1)
    allocate(tempMap(0:l_m-1))
    ! Now copy all defaults across ...
    call copyURIMapping(nsDict%defaults, tempMap, l_m-1)
    !And remove tail-end charlie
    deallocate(nsDict%defaults(l_m)%URI)
    deallocate(nsDict%defaults)
    l_m = l_m - 1
    allocate(nsDict%defaults(0:l_m))
    !Now copy everything back ...
    call copyURIMapping(tempMap, nsDict%defaults, l_m)
    deallocate(tempMap)

  end subroutine removeDefaultNS

  subroutine removePrefixedURI(nsPrefix)
    type(PrefixMapping), intent(inout) :: nsPrefix

    type(URIMapping), dimension(:), allocatable :: tempMap
    integer :: l_m

    l_m = ubound(nsPrefix%urilist,1)
    allocate(tempMap(0:l_m-1))
    ! Now copy all defaults across ...
    call copyURIMapping(nsPrefix%urilist, tempMap, l_m-1)
    !And remove tail-end charlie
    deallocate(nsPrefix%urilist(l_m)%URI)
    deallocate(nsPrefix%urilist)
    l_m = l_m - 1
    allocate(nsPrefix%urilist(0:l_m))
    !Now copy everything back ...
    call copyURIMapping(tempMap, nsPrefix%urilist, l_m)
    deallocate(tempMap)

  end subroutine removePrefixedURI

  subroutine addPrefixedNS(nsDict, prefix, URI, ix, xds, xml, es)
    type(namespaceDictionary), intent(inout) :: nsDict
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: uri
    integer, intent(in) :: ix
    type(xml_doc_state), intent(in) :: xds
    logical, intent(in), optional :: xml
    type(error_stack), intent(inout), optional :: es
    
    integer :: l_p, p_i, i
    logical :: xml_

    if (present(xml)) then
      xml_ = xml
    else
      xml_ = .false.
    endif

    if (prefix=='xml' .and. &
      URI/='http://www.w3.org/XML/1998/namespace') then
      if (present(es)) then
        call add_error(es, "Attempt to assign incorrect URI to prefix 'xml'")
      else
        call FoX_error("Attempt to assign incorrect URI to prefix 'xml'")
      endif
    elseif (prefix/='xml' .and. &
      URI=='http://www.w3.org/XML/1998/namespace') then
      if (present(es)) then
        call add_error(es, "Attempt to assign incorrect prefix to XML namespace")
      else
        call FoX_error("Attempt to assign incorrect prefix to XML namespace")
      endif
    elseif (prefix == 'xmlns') then
      if (present(es)) then
        call add_error(es, "Attempt to declare 'xmlns' prefix")
      else
        call FoX_error("Attempt to declare 'xmlns' prefix")
      endif
    elseif (URI=="http://www.w3.org/2000/xmlns/") then
      if (present(es)) then
        call add_error(es, "Attempt to assign prefix to xmlns namespace")
      else
        call FoX_error("Attempt to assign prefix to xmlns namespace")
      endif
    elseif (len(prefix) > 2) then
      if ((verify(prefix(1:1), 'xX') == 0) &
        .and. (verify(prefix(2:2), 'mM') == 0) &
        .and. (verify(prefix(3:3), 'lL') == 0)) then
        if (.not.xml_) then
          ! FIXME need working warning infrastructure
          !if (present(es)) then
          !  call add_error(es, "Attempt to declare reserved prefix: "//prefix)
          !else
          call FoX_warning("Attempt to declare reserved prefix: "//prefix)
          !endif
        endif
      endif
    endif

    if (.not.checkNCName(prefix, xds%xml_version)) &
      call FoX_error("Attempt to declare invalid prefix: "//prefix)

    ! FIXME check URI is valid

    l_p = ubound(nsDict%prefixes, 1)
    
    p_i = 0
    do i = 1, l_p
       if (str_vs(nsDict%prefixes(i)%prefix) == prefix) then
          p_i = i
          exit
       endif
    enddo

    if (p_i == 0) then
       call addPrefix(nsDict, vs_str(prefix))
       p_i = l_p + 1
    endif
 
    call addPrefixedURI(nsDict%prefixes(p_i), vs_str(URI), ix)

  end subroutine addPrefixedNS

  subroutine removePrefixedNS(nsDict, prefix)
    type(namespaceDictionary), intent(inout) :: nsDict
    character, dimension(:), intent(in) :: prefix
    integer :: l_p, p_i, i
    l_p = ubound(nsDict%prefixes, 1)

    p_i = 0
    do i = 1, l_p
      if (str_vs(nsDict%prefixes(i)%prefix) == str_vs(prefix)) then
        p_i = i
        exit
      endif
    enddo

    if (p_i /= 0) then
      call removePrefixedURI(nsDict%prefixes(p_i))
      if (ubound(nsDict%prefixes(p_i)%urilist,1) == 0) then
        !that was the last mapping for that prefix
        call removePrefix(nsDict, p_i)
      endif
    else
      call FoX_error('Internal error in m_sax_namespaces:removePrefixedNS')
    endif
    
  end subroutine removePrefixedNS

  subroutine addPrefix(nsDict, prefix)
    type(namespaceDictionary), intent(inout) :: nsDict
    character, dimension(:), intent(in) :: prefix
    integer :: l_p

    type(prefixMapping), dimension(:), pointer :: tempPrefixMap

    integer :: i

    !Add a new prefix to the namespace dictionary.
    !Unfortunately this involves copying the entire
    !prefixes dictionary to a temporary structure, then
    !reallocating the prefixes dictionary to be one
    !longer, then copying everything back:

    l_p = ubound(nsDict%prefixes, 1)
    allocate(tempPrefixMap(0:l_p))

    !for each current prefix, append everything to temporary structure
    do i = 0, l_p
       tempPrefixMap(i)%prefix => nsDict%prefixes(i)%prefix
       tempPrefixMap(i)%urilist => nsDict%prefixes(i)%urilist
    enddo
    deallocate(nsDict%prefixes)
    !extend prefix dictionary by one ...
    l_p = l_p + 1
    allocate(nsDict%prefixes(0:l_p))
    !and copy back ...
    do i = 0, l_p-1
       nsDict%prefixes(i)%prefix => tempPrefixMap(i)%prefix
       nsDict%prefixes(i)%urilist => tempPrefixMap(i)%urilist
    enddo
    deallocate(tempPrefixMap)

    allocate(nsDict%prefixes(l_p)%prefix(size(prefix)))
    nsDict%prefixes(l_p)%prefix = prefix
    allocate(nsDict%prefixes(l_p)%urilist(0:0))
    allocate(nsDict%prefixes(l_p)%urilist(0)%URI(len(invalidNS)))
    nsDict%prefixes(l_p)%urilist(0)%URI = vs_str(invalidNS)
    nsDict%prefixes(l_p)%urilist(0)%ix = -1
    
  end subroutine addPrefix

  subroutine removePrefix(nsDict, i_p)
    type(namespaceDictionary), intent(inout) :: nsDict
    integer, intent(in) :: i_p
    integer :: l_p

    type(prefixMapping), dimension(:), pointer :: tempPrefixMap

    integer :: i

    !Remove a prefix from the namespace dictionary.
    !Unfortunately this involves copying the entire
    !prefixes dictionary to a temporary structure, then
    !reallocating the prefixes dictionary to be one
    !shorter, then copying everything back:

    l_p = ubound(nsDict%prefixes, 1)
    allocate(tempPrefixMap(0:l_p-1))

    !for each current prefix, append everything to temporary structure
    do i = 0, i_p-1
       tempPrefixMap(i)%prefix => nsDict%prefixes(i)%prefix
       tempPrefixMap(i)%urilist => nsDict%prefixes(i)%urilist
    enddo
    deallocate(nsDict%prefixes(i_p)%urilist(0)%URI)
    deallocate(nsDict%prefixes(i_p)%urilist)
    deallocate(nsDict%prefixes(i_p)%prefix)
    !this subroutine will only get called if the urilist is already
    !empty, so no need to deallocate it.
    do i = i_p+1, l_p
       tempPrefixMap(i-1)%prefix => nsDict%prefixes(i)%prefix
       tempPrefixMap(i-1)%urilist => nsDict%prefixes(i)%urilist
    enddo
    deallocate(nsDict%prefixes)
    !shorten prefix dictionary by one ...
    l_p = l_p - 1
    allocate(nsDict%prefixes(0:l_p))
    !and copy back ...
    do i = 0, l_p
       nsDict%prefixes(i)%prefix => tempPrefixMap(i)%prefix
       nsDict%prefixes(i)%urilist => tempPrefixMap(i)%urilist
    enddo
    deallocate(tempPrefixMap)

  end subroutine removePrefix


  subroutine checkNamespaces(atts, nsDict, ix, xds, namespace_prefixes, xmlns_uris, es, &
    partial, start_prefix_handler, end_prefix_handler)
    type(dictionary_t), intent(inout) :: atts
    type(namespaceDictionary), intent(inout) :: nsDict
    integer, intent(in) :: ix ! depth of nesting of current element.
    type(xml_doc_state), intent(in) :: xds
    logical, intent(in) :: namespace_prefixes, xmlns_uris
    type(error_stack), intent(inout) :: es
    logical, intent(in) :: partial ! if so, don't try and resolve anything except xml & xmlns
    optional :: start_prefix_handler, end_prefix_handler

    interface
      subroutine start_prefix_handler(namespaceURI, prefix)
        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: prefix
      end subroutine start_prefix_handler
      subroutine end_prefix_handler(prefix)
        character(len=*), intent(in) :: prefix
      end subroutine end_prefix_handler
    end interface

    character(len=6) :: xmlns
    character, dimension(:), pointer :: QName, URIstring
    integer :: i, n
    type(URI), pointer :: URIref
    !Check for namespaces; *and* remove xmlns references from 
    !the attributes dictionary.

    ! we can't do a simple loop across the attributes,
    ! because we need to remove some as we go along ...
    i = 1
    do while (i <= getLength(atts))
      xmlns = get_key(atts, i)
      if (xmlns == 'xmlns ') then
        !Default namespace is being set
        URIstring => vs_str_alloc(get_value(atts, i))
        if (str_vs(URIstring)=="") then
          ! Empty nsURI on default namespace has same effect in 1.0 and 1.1
          if (present(end_prefix_handler)) &
            call end_prefix_handler("")
          call addDefaultNS(nsDict, invalidNS, ix)
          deallocate(URIstring)
        else
          URIref => parseURI(str_vs(URIstring))
          if (.not.associated(URIref)) then
            call add_error(es, "Invalid URI: "//str_vs(URIstring))
            deallocate(URIstring)
            return
          elseif (.not.hasScheme(URIref)) then
            call add_error(es, "Relative namespace in URI deprecated: "//str_vs(URIstring))
            deallocate(URIstring)
            call destroyURI(URIref)
            return
          endif
          call destroyURI(URIref)
          if (present(start_prefix_handler)) &
            call start_prefix_handler(str_vs(URIstring), "")
          call addDefaultNS(nsDict, str_vs(URIstring), ix)
          deallocate(URIstring)
        endif
        if (namespace_prefixes) then
          i = i + 1
        else
          call remove_key(atts, i)
        endif
      elseif (xmlns == 'xmlns:') then
        !Prefixed namespace is being set
        QName => vs_str_alloc(get_key(atts, i))
        URIstring => vs_str_alloc(get_value(atts, i))
        if (str_vs(URIstring)=="") then
          if (xds%xml_version==XML1_0) then
            call add_error(es, "Empty nsURI is invalid in XML 1.0")
            deallocate(URIstring)
            deallocate(QName)
            return
          elseif (xds%xml_version==XML1_1) then
            call addPrefixedNS(nsDict, str_vs(QName(7:)), invalidNS, ix, xds, es=es)
            if (in_error(es)) then
              deallocate(URIstring)
              deallocate(QName)
              return
            elseif (present(end_prefix_handler)) then
              call end_prefix_handler(str_vs(QName(7:)))
            endif
            deallocate(URIstring)
            deallocate(QName)
          endif
        else
          URIref => parseURI(str_vs(URIstring))
          if (.not.associated(URIref)) then
            call add_error(es, "Invalid URI: "//str_vs(URIstring))
            deallocate(URIstring)
            deallocate(QName)
            return
          elseif (.not.hasScheme(URIref)) then
            call add_error(es, "Relative namespace in URI deprecated: "//str_vs(URIstring))
            deallocate(URIstring)
            deallocate(QName)
            call destroyURI(URIref)
            return
          endif
          call destroyURI(URIref)
          call addPrefixedNS(nsDict, str_vs(QName(7:)), str_vs(URIstring), ix, xds, es=es)
          if (in_error(es)) then
            deallocate(URIstring)
            deallocate(QName)
            return
          elseif (present(start_prefix_handler)) then
            call start_prefix_handler(str_vs(URIstring), str_vs(QName(7:)))
          endif
          deallocate(URIstring)
          deallocate(QName)
        endif
        if (namespace_prefixes) then
          i = i + 1
        else
          call remove_key(atts, i)
        endif
      else
        ! we only increment if we haven't removed a key
        i = i + 1
      endif
    enddo

    ! having done that, now resolve all attribute namespaces:
    do i = 1, getLength(atts)
      QName => vs_str_alloc(get_key(atts,i))
      n = index(str_vs(QName), ":")
      if (n > 0) then
        if (str_vs(QName(1:n-1))=="xmlns") then
          ! FIXME but this can be controlled by SAX configuration xmlns-uris
          if (xmlns_uris) then
            call set_nsURI(atts, i, "http://www.w3.org/2000/xmlns/")
          else
            call set_nsURI(atts, i, "")
          endif
        else
          if (str_vs(QName(1:n-1))=="xml") then
            call set_nsURI(atts, i, "http://www.w3.org/XML/1998/namespace")
          elseif (getnamespaceURI(nsDict, str_vs(QName(1:n-1)))==invalidNS) then
            ! Sometimes we don't want to worry about unbound prefixes,
            ! eg if we are in the middle of parsing an entity.
            if (.not.partial) then
              call add_error(es, "Unbound namespace prefix")
              deallocate(QName)
              return
            else
              call set_nsURI(atts, i, "")
            endif
          else
            call set_nsURI(atts, i, getnamespaceURI(nsDict, str_vs(QName(1:n-1))))
          endif
        endif
      else
        if (xmlns_uris.and.str_vs(QName)=="xmlns") then
          call set_nsURI(atts, i, "http://www.w3.org/2000/xmlns/")
        else
          call set_nsURI(atts, i, "") ! no such thing as a default namespace on attributes
        endif
      endif
      ! Check for duplicates
      if (hasKey(atts, getnamespaceURI(nsDict, str_vs(QName(1:n-1))), str_vs(QName(n+1:)))) then
        call add_error(es, "Duplicate attribute names after namespace processing")
        deallocate(QName)
        return
      endif
      call set_localName(atts, i, QName(n+1:))
      deallocate(QName)
    enddo

  end subroutine checkNamespaces

  
  subroutine checkNamespacesWriting(atts, nsdict, ix)
    type(dictionary_t), intent(inout) :: atts
    type(namespaceDictionary), intent(inout) :: nsDict
    integer, intent(in) :: ix
    ! Read through a list of attributes, check with currently
    ! active namespaces & add any necessary declarations

    integer :: i, i_p, l_d, l_ps, n

    n = getLength(atts) ! we need the length before we fiddle with it

    !Does the default NS need added?
    l_d = ubound(nsDict%defaults,1)
    if (nsDict%defaults(l_d)%ix == ix) then
      !It's not been registered yet:
      call add_item_to_dict(atts, "xmlns", &
           str_vs(nsDict%defaults(l_d)%URI), type="CDATA")
    endif

    !next, add any overdue prefixed NS's in the same way:
    ! there should only ever be one. More would be an error,
    ! but the check should have been done earlier.
    do i_p = 0, ubound(nsDict%prefixes, 1)
      l_ps = ubound(nsDict%prefixes(i_p)%urilist,1)
      if (nsDict%prefixes(i_p)%urilist(l_ps)%ix == ix) then
        call add_item_to_dict(atts, &
             "xmlns:"//str_vs(nsDict%prefixes(i_p)%prefix), &
             str_vs(nsDict%prefixes(i_p)%urilist(l_ps)%URI), &
             type="CDATA")
      endif
    enddo
    

    !Finally, we may have some we've added for attribute QNames
    ! have to get those too:
    do i = 1, getLength(atts)
      ! get prefix, and identify the relevant NS mapping
      i_p = getPrefixIndex(nsDict, get_prefix(atts, i))
      l_ps = ubound(nsDict%prefixes(i_p)%urilist,1)
      !If the index is greater than what it should be:
      if (nsDict%prefixes(i_p)%urilist(l_ps)%ix > ix) then
        !we only just added this, so we need to declare it
        call add_item_to_dict(atts, "xmlns:"//get_prefix(atts, i), &
             str_vs(nsDict%prefixes(i_p)%urilist(l_ps)%URI), &
             type="CDATA")
        !Reset the index to the right value:
        nsDict%prefixes(i_p)%urilist(l_ps)%ix = ix
      endif
    enddo

  end subroutine checkNamespacesWriting


  subroutine checkEndNamespaces(nsDict, ix, end_prefix_handler)
    type(namespaceDictionary), intent(inout) :: nsDict
    integer, intent(in) :: ix

    optional :: end_prefix_handler
    
    interface
      subroutine end_prefix_handler(prefix)
        character(len=*), intent(in) :: prefix
      end subroutine end_prefix_handler
    end interface

    integer :: l_d, l_p, l_ps, i
    character, pointer :: prefix(:)

    !It will only ever be the final elements in the list which
    ! might have expired.

    l_d = ubound(nsDict%defaults,1)
    do while (nsDict%defaults(l_d)%ix == ix)
      if (present(end_prefix_handler)) &
        call end_prefix_handler("")
      call removeDefaultNS(nsDict)
      l_d = ubound(nsDict%defaults,1)
    enddo

    l_p = ubound(nsDict%prefixes, 1)
    i = 1
    do while (i <= l_p)
      l_ps = ubound(nsDict%prefixes(l_p)%urilist,1)
      if (nsDict%prefixes(i)%urilist(l_ps)%ix == ix) then
        if (present(end_prefix_handler)) &
          call end_prefix_handler(str_vs(nsDict%prefixes(i)%prefix))
        ! We have to assign this pointer explicitly, otherwise the next call
        ! aliases its arguments illegally.
        prefix =>  nsDict%prefixes(i)%prefix
        call removePrefixedNS(nsDict, prefix)
        if (l_p > ubound(nsDict%prefixes, 1)) then
          ! we just removed the last reference to that prefix,
          ! so our list of prefixes has shrunk - update the running total.
          ! and go to the next prefix, which is at the same index
          l_p = l_p - 1
          cycle
        endif
      endif
      i = i + 1
    enddo

  end subroutine checkEndNamespaces


  subroutine dumpnsdict(nsdict)
    type(namespaceDictionary), intent(in) :: nsdict
    integer :: i, j
    write(*,'(a)')'* default namespaces *'

    do i = 1, ubound(nsdict%defaults, 1)
      write(*,'(i0,a)') nsdict%defaults(i)%ix, str_vs(nsdict%defaults(i)%URI)
    enddo
    write(*,'(a)') '* Prefixed namespaces *'
    do i = 1, ubound(nsdict%prefixes, 1)
       write(*,'(2a)') '* prefix: ', str_vs(nsdict%prefixes(i)%prefix)
       do j = 1, ubound(nsdict%prefixes(i)%urilist, 1)
          write(*,'(i0,a)') nsdict%prefixes(i)%urilist(j)%ix, str_vs(nsdict%prefixes(i)%urilist(j)%URI)
       enddo
    enddo

  end subroutine dumpnsdict


  pure function getURIofDefaultNS(nsDict) result(uri)
    type(namespaceDictionary), intent(in) :: nsDict
    character(len=size(nsDict%defaults(ubound(nsDict%defaults,1))%URI)) :: URI
    
    integer :: l_d
    l_d = ubound(nsDict%defaults,1)
    uri = str_vs(nsDict%defaults(l_d)%URI)
  end function getURIofDefaultNS


  pure function isPrefixInForce(nsDict, prefix) result(force)
    type(namespaceDictionary), intent(in) :: nsDict
    character(len=*), intent(in) :: prefix
    logical :: force
    integer :: i, l_s

    force = .false.
    do i = 1, ubound(nsDict%prefixes, 1)
       if (str_vs(nsDict%prefixes(i)%prefix) == prefix) then
         l_s = ubound(nsDict%prefixes(i)%urilist, 1)
         force = (size(nsdict%prefixes(i)%urilist(l_s)%URI) > 0)
         exit
       endif
    enddo

  end function isPrefixInForce


  pure function isDefaultNSInForce(nsDict) result(force)
    type(namespaceDictionary), intent(in) :: nsDict
    logical :: force
    integer :: l_s

    force = .false.
    l_s = ubound(nsDict%defaults, 1)
    if (l_s > 0) &
      force = (size(nsdict%defaults(l_s)%URI) > 0)

  end function isDefaultNSInForce


  pure function getPrefixIndex(nsDict, prefix) result(p)
    type(namespaceDictionary), intent(in) :: nsDict
    character(len=*), intent(in) :: prefix
    integer :: p
    
    integer :: i
    p = 0
    do i = 1, ubound(nsDict%prefixes, 1)
      if (str_vs(nsDict%prefixes(i)%prefix) == prefix) then
           p = i
           exit
       endif
    enddo
  end function getPrefixIndex

  
  function getNumberOfPrefixes(nsDict) result(n)
    type(namespaceDictionary), intent(in) :: nsDict
    integer :: n
    n = ubound(nsDict%prefixes, 1)
  end function getNumberOfPrefixes


  function getPrefixByIndex(nsDict, i) result(c)
    type(namespaceDictionary), intent(in) :: nsDict
    integer, intent(in) :: i
    character(len=size(nsDict%prefixes(i)%prefix)) :: c

    c = str_vs(nsDict%prefixes(i)%prefix)
  end function getPrefixByIndex


  pure function getURIofPrefixedNS(nsDict, prefix) result(uri)
    type(namespaceDictionary), intent(in) :: nsDict
    character(len=*), intent(in) :: prefix
    character(len=size( &
              nsDict%prefixes( &
           getPrefixIndex(nsDict,prefix) &
                             ) &
                     %urilist( &
           ubound(nsDict%prefixes(getPrefixIndex(nsDict,prefix))%urilist, 1) &
                             ) & 
                      %uri)) :: URI
    integer :: p_i, l_m
    p_i = getPrefixIndex(nsDict, prefix)
    l_m = ubound(nsDict%prefixes(p_i)%urilist, 1)
    uri = str_vs(nsDict%prefixes(p_i)%urilist(l_m)%URI)

  end function getURIofPrefixedNS

#endif
end module m_common_namespaces
