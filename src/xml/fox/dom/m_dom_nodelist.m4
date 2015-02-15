TOHW_m_dom_publics(`

  public :: item
  public :: append
  public :: pop_nl
  public :: remove_nl
  public :: destroyNodeList
  
  interface append
    module procedure append_nl
  end interface
  
  interface item
    module procedure item_nl
  end interface

  interface getLength
    module procedure getLength_nl
  end interface getLength
')`'dnl
dnl
TOHW_m_dom_contents(`

  TOHW_function(item_nl, (list, index), np)
    type(NodeList), pointer :: list
    integer, intent(in) :: index
    type(Node), pointer :: np

    if (.not.associated(list)) then
      TOHW_m_dom_throw_error(FoX_LIST_IS_NULL)
    endif

    if (index>=0.and.index<list%length)  then
      np => list%nodes(index+1)%this
    else
      np => null()
    endif

  end function item_nl

  subroutine append_nl(list, arg)
    type(NodeList), intent(inout) :: list
    type(Node), pointer :: arg

    type(ListNode), pointer :: temp_nl(:)
    integer :: i

    if (.not.associated(list%nodes)) then
      allocate(list%nodes(1))
      list%nodes(1)%this => arg
      list%length = 1
    else
      temp_nl => list%nodes
      allocate(list%nodes(size(temp_nl)+1))
      do i = 1, size(temp_nl)
        list%nodes(i)%this => temp_nl(i)%this
      enddo
      deallocate(temp_nl)
      list%nodes(size(list%nodes))%this => arg
      list%length = size(list%nodes)
    endif
    
  end subroutine append_nl

  TOHW_function(pop_nl, (list), np)
    type(NodeList), pointer :: list
    type(Node), pointer :: np

    type(ListNode), pointer :: temp_nl(:)
    integer :: i

    if (list%length==0) then
      TOHW_m_dom_throw_error(FoX_INTERNAL_ERROR)
    endif

    np => list%nodes(size(list%nodes))%this

    if (list%length==1) then
      deallocate(list%nodes)
      list%length = 0
    else
      temp_nl => list%nodes
      allocate(list%nodes(size(temp_nl)-1))
      do i = 1, size(temp_nl)-1
        list%nodes(i)%this => temp_nl(i)%this
      enddo
      deallocate(temp_nl)
      list%length = size(list%nodes)
    endif
    
  end function pop_nl


  TOHW_function(remove_nl, (nl, index), np)
    type(NodeList), intent(inout) :: nl
    integer, intent(in) :: index
    type(Node), pointer :: np

    type(ListNode), pointer :: temp_nl(:)

    integer :: i

    if (index>nl%length) then
      TOHW_m_dom_throw_error(FoX_INTERNAL_ERROR)
    endif

    np => nl%nodes(index)%this
    temp_nl => nl%nodes
    allocate(nl%nodes(size(temp_nl)-1))
    nl%length = nl%length - 1 
    do i = 1, index - 1
      nl%nodes(i)%this => temp_nl(i)%this
    enddo
    do i = index, nl%length
      nl%nodes(i)%this => temp_nl(i+1)%this
    enddo
    deallocate(temp_nl)

  end function remove_nl


  subroutine remove_node_nl(nl, np)
    type(NodeList), intent(inout) :: nl
    type(Node), pointer :: np

    integer :: i

    do i = 1, nl%length
      if (associated(nl%nodes(i)%this, np)) exit
    enddo
    np => remove_nl(nl, i)

  end subroutine remove_node_nl


  TOHW_function(getLength_nl, (nl), n)
    type(NodeList), pointer :: nl
    integer :: n

    if (.not.associated(nl)) then
      TOHW_m_dom_throw_error(FoX_LIST_IS_NULL)
    endif

    n = size(nl%nodes)
  end function getLength_nl

  subroutine destroyNodeList(nl)
    type(NodeList), pointer :: nl
    
    if (associated(nl%nodes)) deallocate(nl%nodes)
    if (associated(nl%nodeName)) deallocate(nl%nodeName)
    if (associated(nl%localName)) deallocate(nl%localName)
    if (associated(nl%namespaceURI)) deallocate(nl%namespaceURI)
    deallocate(nl)
  end subroutine destroyNodeList

  subroutine updateNodeLists(doc)
    ! When triggered, update all nodelists
    type(Node), pointer :: doc

    type(NodeList), pointer :: nl, nl_orig
    type(NodeListPtr), pointer :: temp_nll(:)
    integer :: i, i_t

    if (.not.getGCstate(doc)) return
    if (.not.doc%docExtras%liveNodeLists) return
    if (.not.associated(doc%docExtras%nodelists)) return

    ! We point the old list of nodelists to temp_nll, then recalculate 
    ! them all (which repopulates nodelists)
    temp_nll => doc%docExtras%nodelists
    i_t = size(temp_nll)
    allocate(doc%docExtras%nodelists(0))
    do i = 1, i_t
      nl_orig => temp_nll(i)%this
      !
      ! Although all nodes should be searched whatever the result,
      ! we should only do the appropriate sort of search for this
      ! list - according to namespaces or not.
      !
      if (associated(nl_orig%nodeName)) then 
        ! this was made by getElementsByTagName
        nl => getElementsByTagName(nl_orig%element, str_vs(nl_orig%nodeName))
      elseif (associated(nl_orig%namespaceURI)) then 
        ! this was made by getElementsByTagNameNS
        nl => getElementsByTagNameNS(nl_orig%element, &
          str_vs(nl_orig%localName), str_vs(nl_orig%namespaceURI))
      endif
    enddo
    ! We dont care about the nodelists weve calculated now
    nullify(nl)

    deallocate(temp_nll)    

  end subroutine updateNodeLists

')`'dnl
