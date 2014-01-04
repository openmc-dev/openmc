module m_common_content_model

#ifndef DUMMYLIB
  ! Allow validating the content model of an XML document

  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc, vs_vs_alloc
  implicit none
  private

  integer, parameter :: OP_NULL = 0
  integer, parameter :: OP_EMPTY = 1
  integer, parameter :: OP_ANY = 2
  integer, parameter :: OP_MIXED = 3
  integer, parameter :: OP_NAME = 4
  integer, parameter :: OP_CHOICE = 5
  integer, parameter :: OP_SEQ = 6

  integer, parameter :: REP_NULL = 0
  integer, parameter :: REP_ONCE = 1
  integer, parameter :: REP_QUESTION_MARK = 2
  integer, parameter :: REP_ASTERISK = 3
  
  type content_particle_t
    character, pointer :: name(:) => null()
    integer :: operator = OP_NULL
    integer :: repeater = REP_NULL
    type(content_particle_t), pointer :: nextSibling => null()
    type(content_particle_t), pointer :: parent => null()
    type(content_particle_t), pointer :: firstChild => null()
  end type content_particle_t

  public :: content_particle_t

  public :: newCP
  public :: transformCPPlus
  public :: checkCP
  public :: checkCPToEnd
  public :: elementContentCP
  public :: emptyContentCP
  public :: destroyCPtree
  public :: dumpCPtree

  public :: OP_NULL, OP_NAME, OP_MIXED, OP_CHOICE, OP_SEQ
  public :: REP_QUESTION_MARK, REP_ASTERISK

contains

  function newCP(empty, any, name, repeat) result(cp)
    logical, intent(in), optional :: empty
    logical, intent(in), optional :: any
    character(len=*), intent(in), optional :: name
    character, intent(in), optional :: repeat
    type(content_particle_t), pointer :: cp

    allocate(cp)
    if (present(empty)) then
      cp%operator = OP_EMPTY
    elseif (present(any)) then
      cp%operator = OP_ANY
    elseif (present(name)) then
      cp%operator = OP_NAME
      cp%name => vs_str_alloc(name)
    else
      cp%operator = OP_SEQ
    endif
    if (present(repeat)) then
      select case (repeat)
      case("?")
        cp%repeater = REP_QUESTION_MARK
      case("*")
        cp%repeater = REP_ASTERISK
      end select
    endif

  end function newCP

  function copyCP(cp) result(cp_out)
    type(content_particle_t), pointer :: cp
    type(content_particle_t), pointer :: cp_out

    allocate(cp_out)
    if (associated(cp%name)) cp_out%name => vs_vs_alloc(cp%name)
    cp_out%operator = cp%operator
    cp_out%repeater = cp%repeater
  end function copyCP

  function copyCPtree(cp) result(cp_out)
    type(content_particle_t), pointer :: cp
    type(content_particle_t), pointer :: cp_out

    type(content_particle_t), pointer :: tcp, tcp_out, tcpn_out, tcpp_out
    logical :: done

    tcp => cp
    cp_out => copyCP(cp)
    tcp_out => cp_out
    done = .false.
    do while (associated(tcp_out))
      if (.not.done) then
        do while (associated(tcp%firstChild))
          tcp => tcp%firstChild
          tcpn_out => copyCP(tcp)
          tcp_out%firstChild => tcpn_out
          tcpn_out%parent => tcp_out
          tcp_out => tcpn_out
        enddo
      endif
      tcpp_out => tcp_out%parent
      if (associated(tcp%nextSibling)) then
        done = .false.
        tcp => tcp%nextSibling
        tcpn_out => copyCP(tcp)
        tcp_out%nextSibling => tcpn_out
        tcpn_out%parent => tcpp_out
        tcp_out => tcpn_out
      else
        done = .true.
        tcp => tcp%parent
        tcp_out => tcp_out%parent
      endif
    enddo
  end function copyCPtree

  subroutine transformCPPlus(cp)
    type(content_particle_t), pointer :: cp

    type(content_particle_t), pointer :: tcp, cp_new

    ! Make copy of cp, and graft children on
    cp_new => copyCP(cp)
    cp_new%firstChild => cp%firstChild

    ! Reset children's parents ...
    tcp => cp%firstChild
    do while (associated(tcp))
      tcp%parent => cp_new
      tcp => tcp%nextSibling
    enddo

    ! Clear cp & make it an SEQ
    if (associated(cp%name)) deallocate(cp%name)
    cp%operator = OP_SEQ

    ! Append our copied cp to the now-an-SEQ
    cp%firstChild => cp_new
    cp_new%parent => cp

    ! Copy it for a sibling, and make the sibling a *
    cp_new%nextSibling => copyCPtree(cp_new)
    cp_new%nextSibling%parent => cp
    cp_new%nextSibling%repeater = REP_ASTERISK

  end subroutine transformCPPlus

  function checkCP(cp, name) result(p)
    type(content_particle_t), pointer :: cp
    character(len=*), intent(in) :: name
    logical :: p

    type(content_particle_t), pointer :: tcp

    ! for EMPTY, ANY or MIXED, cp never moves.
    ! for element content, we move the pointer as we
    ! move through the regex.

    ! If the regex includes ambiguous content, we are
    ! a bit screwed. But the document is in error if so.
    ! (and we are not required to diagnose errors.)

    p = .false.
    if (.not.associated(cp)) return

    select case(cp%operator)
    case (OP_EMPTY)
      continue ! anything fails
    case (OP_ANY)
      p = .true.
    case (OP_MIXED)
      tcp => cp%firstChild
      do while (associated(tcp))
        if (name==str_vs(tcp%name)) then
          p = .true.
          exit
        endif
        tcp => tcp%nextSibling
      enddo
    case default
      do
        if (.not.associated(cp)) exit
        select case (cp%operator)
        case (OP_NAME)
          p = (name==str_vs(cp%name))
          if (p) then
            tcp => nextCPAfterMatch(cp)
            cp => tcp
            exit
          else
            tcp => nextCPAfterFail(cp)
            cp => tcp
          endif
        case (OP_CHOICE, OP_SEQ)
          cp => cp%firstChild
        end select
      end do
    end select

  end function checkCP

  function nextCPaftermatch(cp) result(cp_next)
    type (content_particle_t), pointer :: cp
    type (content_particle_t), pointer :: cp_next

    type (content_particle_t), pointer :: tcp

    cp_next => cp

    do
      if (cp_next%repeater==REP_ASTERISK) exit
      tcp => cp_next%parent
      if (associated(tcp)) then
        if (tcp%operator==OP_CHOICE) then
          ! siblings are uninteresting, we've matched this CHOICE
          cp_next => tcp
          ! Move up & try the whole thing again on the parent CHOICE
        elseif (tcp%operator==OP_SEQ) then
          ! we do care about siblings, move onto next one
          if (associated(cp_next%nextSibling)) then
            cp_next => cp_next%nextSibling
            ! thatll do, itll be the next thing to try
            exit
          else
            ! No sibling, move up a level
            cp_next => tcp
            ! and try again
          endif
        endif
      else
        ! We've got to the top already.
        cp_next => tcp
        exit
      endif
    enddo
        
  end function nextCPaftermatch

  function nextCPafterfail(cp) result(cp_next)
    type(content_particle_t), pointer :: cp
    type(content_particle_t), pointer :: cp_next

    type(content_particle_t), pointer :: tcp
    logical :: match

    match = .false.
    cp_next => cp
    do
      tcp => cp_next%parent
      if (associated(tcp)) then
        if (tcp%operator==OP_CHOICE) then
          ! we care about siblings, lets try the next one
          if (associated(cp_next%nextSibling)) then
            cp_next => cp_next%nextSibling
            ! super, lets go back and try that
            exit
          else ! weve failed to match any cp in this CHOICE ...
            cp_next => tcp
            ! go up a level and see if theres another legitimate choice
          endif
        elseif (tcp%operator==OP_SEQ) then
          if ((match.or.cp_next%repeater/=REP_NULL) &
            .and.associated(cp_next%nextSibling)) then
            ! we were allowed to fail to match, try sibling
            cp_next => cp_next%nextSibling
            exit
          elseif (cp_next%repeater/=REP_NULL) then
            match = .true.
            ! The last item was optional, so weve matched at this level
            cp_next => tcp
          elseif (associated(tcp%firstChild, cp_next)) then
            ! we havent matched - but we hadnt started, Maybe it was ok
            ! not to match because we are nested inside an optional thingy
            cp_next => tcp
          else
            ! We were not allowed to fail there,
            ! there is no legitimate next choice.
            cp_next => null()
            exit
          endif
        endif
      else
        ! weve got all the way to the top without
        ! finding a new cp to try. But if this top-level
        ! cp is ASTERISK'ed we can try it agin
        cp_next => null()
        exit
      endif
    enddo
    
  end function nextCPafterfail

  function checkCPToEnd(cp) result(p)
    type(content_particle_t), pointer :: cp
    logical :: p

    type(content_particle_t), pointer :: tcp

    if (associated(cp)) then
      select case(cp%operator)
      case (OP_EMPTY, OP_ANY, OP_MIXED)
        p = .true.
      case default
        tcp => nextCPMustMatch(cp)
        p = .not.associated(tcp)
      end select
    else
      p = .true.
    endif
  end function checkCPToEnd

  function nextCPMustMatch(cp) result(cp_next)
    type(content_particle_t), pointer :: cp
    type(content_particle_t), pointer :: cp_next

    type(content_particle_t), pointer :: tcp

    if (.not.associated(cp)) return
    if (.not.associated(cp%parent)) then
      ! we havent started exploring this one.
      ! get the first starting position
      cp_next => cp
      do while (cp_next%repeater==REP_NULL)
        if (associated(cp_next%firstChild)) then
          cp_next => cp_next%firstChild
        else
          exit
        endif
      enddo
    else
      cp_next => cp
    endif
    if (cp_next%repeater==REP_NULL) return
    do
      tcp => cp_next%parent
      if (associated(tcp)) then
        if (tcp%operator==OP_CHOICE) then
          ! its matched by the optional one we are on, go up a level
          cp_next => tcp
        elseif (tcp%operator==OP_SEQ) then
          ! check all siblings for any compulsory ones
          do while (associated(cp_next%nextSibling))
            cp_next => cp_next%nextSibling
            if (cp_next%repeater==REP_NULL) return
          enddo
          ! all were optional, go up a level
          cp_next => tcp
        endif
      else
        ! weve got all the way to the top without
        ! finding a new cp to try
        cp_next => tcp
        exit
      endif
    enddo
    
  end function nextCPMustMatch

  function elementContentCP(cp) result(p)
    type(content_particle_t), pointer :: cp
    logical :: p

    if (associated(cp)) then
      select case (cp%operator)
      case (OP_EMPTY, OP_ANY, OP_MIXED)
        p = .false.
      case default
        p = .true.
      end select
    else
      p = .true.
    endif

  end function elementContentCP

  function emptyContentCP(cp) result(p)
    type(content_particle_t), pointer :: cp
    logical :: p

    if (associated(cp)) then
      p = cp%operator==OP_EMPTY
    else
      p = .false.
    endif

  end function emptyContentCP

  subroutine destroyCP(cp)
    type(content_particle_t), pointer :: cp

    if (associated(cp%name)) deallocate(cp%name)
    deallocate(cp)
  end subroutine destroyCP

  subroutine destroyCPtree(cp)
    type(content_particle_t), pointer :: cp

    type(content_particle_t), pointer :: current, tcp

    current => cp
    do
      do while (associated(current%firstChild))
        current => current%firstChild
      enddo
      if (associated(current, cp)) exit
      tcp => current
      if (associated(current%nextSibling)) then
        current => current%nextSibling
        call destroyCP(tcp)
      else
        current => current%parent
        call destroyCP(tcp)
        current%firstChild => null()
      endif
    enddo
    call destroyCP(cp)

  end subroutine destroyCPtree

  subroutine dumpCP(cp)
    type(content_particle_t), pointer :: cp

    select case(cp%operator)
    case (OP_EMPTY)
      write(*,'(a)', advance="no") "EMPTY"
    case (OP_ANY)
      write(*,'(a)', advance="no") "ANY"
    case (OP_MIXED)
      write(*,'(a)', advance="no") "MIXED"
    case (OP_NAME)
      write(*,'(a)', advance="no") str_vs(cp%name)
    case (OP_CHOICE)
      write(*,'(a)', advance="no") "CHOICE"
    case (OP_SEQ)
      write(*,'(a)', advance="no") "SEQ"
    end select
    select case(cp%repeater)
    case (REP_QUESTION_MARK)
      write(*,'(a)', advance="no") "?"
    case (REP_ASTERISK)
      write(*,'(a)', advance="no") "*"
    end select
    write(*,*)
  end subroutine dumpCP

  subroutine dumpCPtree(cp)
    type(content_particle_t), pointer :: cp

    type(content_particle_t), pointer :: current

    integer :: i
    logical :: done
    i = 0
    current => cp
    done = .false.
    call dumpCP(current)
    do
      if (.not.done) then
        do while (associated(current%firstChild))
          i = i + 2
          current => current%firstChild
          write(*,'(a)', advance="no") repeat(" ",i)
          call dumpCP(current)
        enddo
      endif
      if (associated(current, cp)) exit
      if (associated(current%nextSibling)) then
        done = .false.
        current => current%nextSibling
        write(*,'(a)', advance="no") repeat(" ",i)
        call dumpCP(current)
      else
        done = .true.
        i = i - 2
        current => current%parent
      endif
    enddo

  end subroutine dumpCPtree

#endif

end module m_common_content_model
