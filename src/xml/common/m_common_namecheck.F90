module m_common_namecheck

#ifndef DUMMYLIB
  ! These are basically a collection of what would be regular
  ! expressions in a more sensible language.
  ! The only external dependency should be knowing how these
  ! regular expressions may differ between XML-1.0 and 1.1 (which
  ! is only in the areas of
  ! 1: allowing character entity references to control characters
  ! 2: More characters allowed in Names (but this only affects
  !    unicode-aware programs, so is only skeleton here)

  use fox_m_fsys_format, only: str_to_int_10, str_to_int_16, operator(//)
  use fox_m_fsys_string, only: toLower
  use m_common_charset, only: isLegalCharRef, isNCNameChar, &
    isInitialNCNameChar, isInitialNameChar, isNameChar, isRepCharRef

  implicit none
  private

  character(len=*), parameter :: lowerCase = "abcdefghijklmnopqrstuvwxyz"
  character(len=*), parameter :: upperCase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(len=*), parameter :: letters = lowerCase//upperCase
  character(len=*), parameter :: digits = "0123456789"
  character(len=*), parameter :: hexdigits = "0123456789abcdefABCDEF"
  character(len=*), parameter :: NameChars = lowerCase//upperCase//digits//".-_:"

  public :: checkName
  public :: checkNames
  public :: checkQName
  public :: checkQNames
  public :: checkNmtoken
  public :: checkNmtokens
  public :: checkNCName
  public :: checkNCNames
  public :: checkEncName
  public :: checkPITarget
  public :: checkPublicId
  public :: checkPEDef
  public :: checkPseudoAttValue
  public :: checkAttValue
  public :: checkCharacterEntityReference
  public :: checkRepCharEntityReference
  public :: likeCharacterEntityReference

  public :: prefixOfQName
  public :: localpartOfQName

contains

  pure function checkEncName(name) result(good)
    ![81]   	EncName	   ::=   	[A-Za-z] ([A-Za-z0-9._] | '-')*
    character(len=*), intent(in) :: name
    logical :: good

    integer :: n
    
    n = len(name)
    good = (n > 0)
    if (good) good = (scan(name(1:1), letters) /= 0)
    if (good .and. n > 1) &
         good = (verify(name(2:), letters//digits//'.-_') == 0)
  end function checkEncName


  function checkPITarget(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the XML requirements for a NAME
    ! Is not fully compliant; ignores UTF issues.

    good = checkName(name, xv) &
      .and.toLower(name)/="xml"

  end function checkPITarget


  pure function checkName(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the XML requirements for a NAME
    ! Is not fully compliant; ignores UTF issues.

    good = (len(name) > 0)
    if (.not.good) return
    if (good) good = isInitialNameChar(name(1:1), xv)
    if (.not.good.or.len(name)==1) return
    good = isNameChar(name(2:), xv)
       
  end function checkName

  pure function checkNames(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the production for NAMES

    integer :: i, j

    good = (len(name) > 0)
    if (.not.good) return
    i = verify(name, " ")
    if (i==0) then
      good = .false.
      return
    endif
    j = scan(name(i:), " ")
    if (j==0) then
      j = len(name)
    else
      j = i + j - 2
    endif
    do
      good = checkName(name(i:j), xv)
      if (.not.good) return
      i = j + 1
      j = verify(name(i:), " ")
      if (j==0) exit
      i = i + j - 1
      j = scan(name(i:), " ")
      if (j==0) then
        j = len(name)
      else
        j = i + j - 2
      endif
    enddo
       
  end function checkNames


  pure function checkQName(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the XML requirements for a NAME
    ! Is not fully compliant; ignores UTF issues.

    integer :: n

    n = index(name, ':')
    if (n == 0) then
      good = checkNCName(name, xv)
    else
      good = (checkNCName(name(:n-1), xv) .and. checkNCName(name(n+1:), xv))
    endif
  end function checkQName


  pure function checkQNames(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the production for NAMES

    integer :: i, j

    good = (len(name) > 0)
    if (.not.good) return
    i = verify(name, " ")
    if (i==0) then
      good = .false.
      return
    endif
    j = scan(name(i:), " ")
    if (j==0) then
      j = len(name)
    else
      j = i + j - 2
    endif
    do
      good = checkQName(name(i:j), xv)
      if (.not.good) return
      i = j + 1
      j = verify(name(i:), " ")
      if (j==0) exit
      i = i + j - 1
      j = scan(name(i:), " ")
      if (j==0) then
        j = len(name)
      else
        j = i + j - 2
      endif
    enddo
       
  end function checkQNames


  pure function checkNCName(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the XML requirements for an NCNAME
    ! Is not fully compliant; ignores UTF issues.

    good = (len(name)/=0)
    if (.not.good) return
    good = isInitialNCNameChar(name(1:1), xv)
    if (.not.good.or.len(name)==1) return
    good = isNCNameChar(name(2:), xv)
       
  end function checkNCName


  pure function checkNCNames(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the production for NAMES

    integer :: i, j

    good = (len(name) > 0)
    if (.not.good) return
    i = verify(name, " ")
    if (i==0) then
      good = .false.
      return
    endif
    j = scan(name(i:), " ")
    if (j==0) then
      j = len(name)
    else
      j = i + j - 2
    endif
    do
      good = checkNCName(name(i:j), xv)
      if (.not.good) return
      i = j + 1
      j = verify(name(i:), " ")
      if (j==0) exit
      i = i + j - 1
      j = scan(name(i:), " ")
      if (j==0) then
        j = len(name)
      else
        j = i + j - 2
      endif
    enddo
       
  end function checkNCNames


  pure function checkNmtoken(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the XML requirements for an NCNAME
    ! Is not fully compliant; ignores UTF issues.

    good = isNameChar(name, xv)
       
  end function checkNmtoken


  pure function checkNmtokens(name, xv) result(good)
    character(len=*), intent(in) :: name
    integer, intent(in) :: xv
    logical :: good
    ! Validates a string against the XML requirements for an NCNAME
    ! Is not fully compliant; ignores UTF issues.

    integer :: i, j

    good = (len(name) > 0)
    if (.not.good) return
    i = verify(name, " ")
    if (i==0) then
      good = .false.
      return
    endif
    j = scan(name(i:), " ")
    if (j==0) then
      j = len(name)
    else
      j = i + j - 2
    endif
    do
      good = isNameChar(name(i:j), xv)
      if (.not.good) return
      i = j + 1
      j = verify(name(i:), " ")
      if (j==0) exit
      i = i + j - 1
      j = scan(name(i:), " ")
      if (j==0) then
        j = len(name)
      else
        j = i + j - 2
      endif
    enddo
       
  end function checkNmtokens


  function checkPublicId(value) result(good)
    character(len=*), intent(in) :: value
    logical :: good
    character(len=*), parameter :: PubIdChars = &
      " "//achar(10)//achar(13)//lowerCase//upperCase//digits//"-'()+,./:=?;!*#@$_%"

    good = (verify(value, PubIdChars)==0) 
  end function checkPublicId


  function checkPEDef(value, xv) result(p)
    character(len=*), intent(in) :: value
    integer, intent(in) :: xv
    logical :: p

    integer :: i1, i2

    p = .false.
    if (scan(value, '%&')==0) then
      p = .true.
    elseif (scan(value, '"')==0) then
      i1 = scan(value, '%&')
      i2 = 0
      do while (i1>0)
        i1 = i2 + i1
        i2 = index(value(i1+1:),';')
        if (i2==0) return
        i2 = i1 + i2
        if (value(i1:i1)=='&') then
          if (.not.checkName(value(i1+1:i2-1), xv) .and. &
            .not.checkCharacterEntityReference(value(i1+1:i2-1), xv)) return
        else
          if (.not.checkName(value(i1+1:i2-1), xv)) &
            return
        endif
        i1 = scan(value(i2+1:), '%&')
      enddo
      p = .true.
    endif
  end function checkPEDef

  function checkPseudoAttValue(value, xv) result(p)
    character(len=*), intent(in) :: value
    integer, intent(in) :: xv
    logical :: p

    integer :: i1, i2
!fixme can we have entrefs in PIs?
    p = .false.
    if (scan(value, '"<&')==0) then
      p = .true.
    elseif (index(value, '&') > 0) then
      i1 = index(value, '&')
      i2 = 0
      do while (i1 > 0)
        i1 = i2 + i1
        i2 = index(value(i1+1:),';')
        if (i2==0) return
        i2 = i1 + i2
        if (value(i1+1:i2-1) /= 'amp' .and. &
          value(i1+1:i2-1) /= 'lt' .and. &
          value(i1+1:i2-1) /= 'gt' .and. &
          value(i1+1:i2-1) /= 'quot' .and. &
          value(i1+1:i2-1) /= 'apos' .and. &
          .not.checkCharacterEntityReference(value(i1+1:i2-1), xv)) &
          return
        i1 = index(value(i2+1:), '&')
      enddo
      p = .true.
    endif
  end function checkPseudoAttValue

  function checkAttValue(value, xv) result(p)
    character(len=*), intent(in) :: value
    integer, intent(in) :: xv
    logical :: p
    ! This function is called basically to make sure
    ! that any attribute value looks like one. It will
    ! not flag up out-of-range character entities, and
    ! is a very weak check. Only used from xml_AddAttribute
    ! when escaping is off.
    integer :: i1, i2

    p = .false.
    if (scan(value, '"<&'//"'")==0) then
      p = .true.
    elseif (index(value, '&') > 0) then
      i1 = index(value, '&')
      i2 = 0
      do while (i1 > 0)
        i1 = i1 + i2 + 1
        i2 = scan(value(i1+1:),';')
        if (i2 == 0) return
        i2 = i1 + i2
        if (.not.checkName(value(i1+1:i2-1), xv) .and. &
          .not.likeCharacterEntityReference(value(i1+1:i2-1))) then
          print*, value(i1+1:i2-1), " ", &
            likeCharacterEntityReference(value(i1+1:i2-1))
          return
        endif
        i1 = index(value(i2+1:), '&')
      enddo
      p = .true.
    endif
  end function checkAttValue

  
  function likeCharacterEntityReference(code) result(good)
    character(len=*), intent(in) :: code
    logical :: good

    good = .false.
    if (len(code) > 0) then
      if (code(1:1) == "#") then
        if (code(2:2) == "x") then
          if (len(code) > 2) then
            good = (verify(code(3:), hexdigits) == 0)
          endif
        else
          good = (verify(code(2:), digits) == 0)
        endif
      endif
    endif

  end function likeCharacterEntityReference

  function checkCharacterEntityReference(code, xv) result(good)
    character(len=*), intent(in) :: code
    integer, intent(in) :: xv
    logical :: good

    integer :: i

    good = .false.
    if (len(code) > 0) then
      if (code(1:1) == "#") then
        if (code(2:2) == "x") then
          if (len(code) > 2) then
            good = (verify(code(3:), hexdigits) == 0)
            if (good) then
              i = str_to_int_16(code(3:))
            endif
          endif
        else
          good = (verify(code(2:), digits) == 0)
          if (good) then
            i = str_to_int_10(code(2:))
          endif
        endif
      endif
    endif
    if (good) good = isLegalCharRef(i, xv)

  end function checkCharacterEntityReference

  function checkRepCharEntityReference(code, xv) result(good)
    character(len=*), intent(in) :: code
    integer, intent(in) :: xv
    logical :: good
    
    ! Is this a reference to a character we can actually represent
    ! in memory? ie without unicode, US-ASCII only.

    integer :: i

    good = .false.
    if (len(code) > 0) then
      if (code(1:1) == "#") then
        if (code(2:2) == "x") then
          if (len(code) > 2) then
            good = (verify(code(3:), hexdigits) == 0)
            if (good) then
              i = str_to_int_16(code(3:))
            endif
          endif
        else
          good = (verify(code(2:), digits) == 0)
          if (good) then
            i = str_to_int_10(code(2:))
          endif
        endif
      endif
    endif
    if (good) good = isRepCharRef(i, xv)

  end function checkRepCharEntityReference
  

  pure function prefixOfQName(qname) result(prefix)
    character(len=*), intent(in) :: qname
    character(len=max(index(qname, ':')-1,0)) :: prefix

    prefix = qname ! automatic truncation
  end function prefixOfQName

  
  pure function localpartOfQname(qname) result(localpart)
    character(len=*), intent(in) :: qname
    character(len=max(len(qname)-index(qname,':'),0)) ::localpart

    localpart = qname(index(qname,':')+1:)
  end function localpartOfQname

#endif
end module m_common_namecheck
