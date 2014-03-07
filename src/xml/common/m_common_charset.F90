module m_common_charset

#ifndef DUMMYLIB
  ! Written to use ASCII charset only. Full UNICODE would
  ! take much more work and need a proper unicode library.

  use fox_m_fsys_string, only: toLower

  implicit none
  private

!!$  character(len=1), parameter :: ASCII = &
!!$achar(0)//achar(1)//achar(2)//achar(3)//achar(4)//achar(5)//achar(6)//achar(7)//achar(8)//achar(9)//&
!!$achar(10)//achar(11)//achar(12)//achar(13)//achar(14)//achar(15)//achar(16)//achar(17)//achar(18)//achar(19)//&
!!$achar(20)//achar(21)//achar(22)//achar(23)//achar(24)//achar(25)//achar(26)//achar(27)//achar(28)//achar(29)//&
!!$achar(30)//achar(31)//achar(32)//achar(33)//achar(34)//achar(35)//achar(36)//achar(37)//achar(38)//achar(39)//&
!!$achar(40)//achar(41)//achar(42)//achar(43)//achar(44)//achar(45)//achar(46)//achar(47)//achar(48)//achar(49)//&
!!$achar(50)//achar(51)//achar(52)//achar(53)//achar(54)//achar(55)//achar(56)//achar(57)//achar(58)//achar(59)//&
!!$achar(60)//achar(61)//achar(62)//achar(63)//achar(64)//achar(65)//achar(66)//achar(67)//achar(68)//achar(69)//&
!!$achar(70)//achar(71)//achar(72)//achar(73)//achar(74)//achar(75)//achar(76)//achar(77)//achar(78)//achar(79)//&
!!$achar(80)//achar(81)//achar(82)//achar(83)//achar(84)//achar(85)//achar(86)//achar(87)//achar(88)//achar(89)//&
!!$achar(90)//achar(91)//achar(92)//achar(93)//achar(94)//achar(95)//achar(96)//achar(97)//achar(98)//achar(99)//&
!!$achar(100)//achar(101)//achar(102)//achar(103)//achar(104)//achar(105)//achar(106)//achar(107)//achar(108)//achar(109)//&
!!$achar(110)//achar(111)//achar(112)//achar(113)//achar(114)//achar(115)//achar(116)//achar(117)//achar(118)//achar(119)//&
!!$achar(120)//achar(121)//achar(122)//achar(123)//achar(124)//achar(125)//achar(126)//achar(127)

  character(len=1), parameter :: SPACE           = achar(32)
  character(len=1), parameter :: NEWLINE         = achar(10)
  character(len=1), parameter :: CARRIAGE_RETURN = achar(13)
  character(len=1), parameter :: TAB             = achar(9)

  character(len=*), parameter :: whitespace = SPACE//NEWLINE//CARRIAGE_RETURN//TAB

  character(len=*), parameter :: lowerCase = "abcdefghijklmnopqrstuvwxyz"
  character(len=*), parameter :: upperCase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(len=*), parameter :: digits = "0123456789"
  character(len=*), parameter :: hexdigits = "0123456789abcdefABCDEF"
  character(len=*), parameter :: InitialNCNameChars = lowerCase//upperCase//"_"
  character(len=*), parameter :: NCNameChars = InitialNCNameChars//digits//".-"
  character(len=*), parameter :: InitialNameChars = InitialNCNameChars//":"
  character(len=*), parameter :: NameChars = NCNameChars//":"

  character(len=*), parameter :: PubIdChars = NameChars//whitespace//"'()+,/=?;!*#@$%"
  character(len=*), parameter :: validchars = &
    whitespace//"!""#$%&'()*+,-./"//digits// &
    ":;<=>?@"//upperCase//"[\]^_`"//lowerCase//"{|}~"
  ! these are all the standard ASCII printable characters: whitespace + (33-126)
  ! which are the only characters we can guarantee to know how to handle properly.

  integer, parameter :: XML1_0 = 10      ! NB 0x7F was legal in XML-1.0, but illegal in XML-1.1

  integer, parameter :: XML1_1 = 11

  character(len=*), parameter :: XML1_0_ILLEGALCHARS = achar(0)
  character(len=*), parameter :: XML1_1_ILLEGALCHARS = NameChars

  character(len=*), parameter :: XML1_0_INITIALNAMECHARS = InitialNameChars
  character(len=*), parameter :: XML1_1_INITIALNAMECHARS = InitialNameChars
  character(len=*), parameter :: XML1_0_NAMECHARS = NameChars
  character(len=*), parameter :: XML1_1_NAMECHARS = NameChars

  character(len=*), parameter :: XML1_0_INITIALNCNAMECHARS = InitialNCNameChars
  character(len=*), parameter :: XML1_1_INITIALNCNAMECHARS = InitialNCNameChars
  character(len=*), parameter :: XML1_0_NCNAMECHARS = NCNameChars
  character(len=*), parameter :: XML1_1_NCNAMECHARS = NCNameChars


  character(len=*), parameter :: XML_WHITESPACE = whitespace
  character(len=*), parameter :: XML_INITIALENCODINGCHARS = lowerCase//upperCase
  character(len=*), parameter :: XML_ENCODINGCHARS = lowerCase//upperCase//digits//'._-'

  public :: validchars
  public :: whitespace
  public :: uppercase

  public :: digits
  public :: hexdigits

  public :: XML1_0
  public :: XML1_1
  public :: XML1_0_NAMECHARS
  public :: XML1_1_NAMECHARS
  public :: XML1_0_INITIALNAMECHARS
  public :: XML1_1_INITIALNAMECHARS
  public :: XML_WHITESPACE
  public :: XML_INITIALENCODINGCHARS
  public :: XML_ENCODINGCHARS

  public :: isLegalChar
  public :: isLegalCharRef
  public :: isRepCharRef
  public :: isInitialNameChar
  public :: isNameChar
  public :: isInitialNCNameChar
  public :: isNCNameChar
  public :: isXML1_0_NameChar
  public :: isXML1_1_NameChar
  public :: checkChars
  public :: isUSASCII
  public :: allowed_encoding

contains

  pure function isLegalChar(c, ascii_p, xml_version) result(p)
    character, intent(in) :: c
    ! really we should check the encoding here & be more intelligent
    ! for now we worry only about is it ascii or not.
    logical, intent(in) :: ascii_p
    integer, intent(in) :: xml_version
    logical :: p 
    ! Is this character legal as a source character in the document?
    integer :: i
    i = iachar(c)
    if (i<0) then
      p = .false.
      return
    elseif (i>127) then
      p = .not.ascii_p
      return
      ! ie if we are ASCII, then >127 is definitely illegal.
      ! otherwise maybe it's ok
    endif
    select case(xml_version)
    case (XML1_0)
      p = (i==9.or.i==10.or.i==13.or.(i>31.and.i<128))
    case (XML1_1)
      p = (i==9.or.i==10.or.i==13.or.(i>31.and.i<127))
      ! NB 0x7F was legal in XML-1.0, but illegal in XML-1.1
    end select
  end function isLegalChar

  pure function isLegalCharRef(i, xml_version) result(p)
    integer, intent(in) :: i
    integer, intent(in) :: xml_version
    logical :: p 

    ! Is Unicode character #i legal as a character reference?

    if (xml_version==XML1_0) then
      p = (i==9).or.(i==10).or.(i==13).or.(i>31.and.i<55296).or.(i>57343.and.i<65534).or.(i>65535.and.i<1114112)
    elseif (xml_version==XML1_1) then
      p = (i>0.and.i<55296).or.(i>57343.and.i<65534).or.(i>65535.and.i<1114112)
      ! XML 1.1 made all control characters legal as character references.
    end if
  end function isLegalCharRef

  pure function isRepCharRef(i, xml_version) result(p)
    integer, intent(in) :: i
    integer, intent(in) :: xml_version
    logical :: p 

    ! Is Unicode character #i legal and representable here?

    if (xml_version==XML1_0) then
      p = (i==9).or.(i==10).or.(i==13).or.(i>31.and.i<128)
    elseif (xml_version==XML1_1) then
      p = (i>0.and.i<128)
      ! XML 1.1 made all control characters legal as character references.
    end if
  end function isRepCharRef

  pure function isInitialNameChar(c, xml_version) result(p)
    character, intent(in) :: c
    integer, intent(in) :: xml_version
    logical :: p

    select case(xml_version)
    case (XML1_0)
      p = (verify(c, XML1_0_INITIALNAMECHARS)==0)
    case (XML1_1)
      p = (verify(c, XML1_1_INITIALNAMECHARS)==0)
    end select

  end function isInitialNameChar

  pure function isNameChar(c, xml_version) result(p)
    character(len=*), intent(in) :: c
    integer, intent(in) :: xml_version
    logical :: p

    select case(xml_version)
    case (XML1_0)
      p = (verify(c, XML1_0_NAMECHARS)==0)
    case (XML1_1)
      p = (verify(c, XML1_1_NAMECHARS)==0)
    end select

  end function isNameChar

  pure function isInitialNCNameChar(c, xml_version) result(p)
    character, intent(in) :: c
    integer, intent(in) :: xml_version
    logical :: p

    select case(xml_version)
    case (XML1_0)
      p = (verify(c, XML1_0_INITIALNCNAMECHARS)==0)
    case (XML1_1)
      p = (verify(c, XML1_1_INITIALNCNAMECHARS)==0)
    end select
  end function isInitialNCNameChar

  pure function isNCNameChar(c, xml_version) result(p)
    character(len=*), intent(in) :: c
    integer, intent(in) :: xml_version
    logical :: p

    select case(xml_version)
    case (XML1_0)
      p = (verify(c, XML1_0_NCNAMECHARS)==0)
    case (XML1_1)
      p = (verify(c, XML1_1_NCNAMECHARS)==0)
    end select
  end function isNCNameChar

  function isXML1_0_NameChar(c) result(p)
    character, intent(in) :: c
    logical :: p

    p = (verify(c, XML1_0_NAMECHARS)==0)

  end function isXML1_0_NameChar

  function isXML1_1_NameChar(c) result(p)
    character, intent(in) :: c
    logical :: p

    p = (verify(c, XML1_1_NAMECHARS)==0)

  end function isXML1_1_NameChar

  pure function checkChars(value, xv) result(p)
    character(len=*), intent(in) :: value
    integer, intent(in) :: xv
    logical :: p

    ! This checks if value only contains values
    ! legal to appear (escaped or unescaped) 
    ! according to whichever XML version is in force.
    integer :: i

    p = .true.
    do i = 1, len(value)
      if (xv == XML1_0) then
        select case(iachar(value(i:i)))
        case (0,8)
          p = .false.
        case (11,12)
          p = .false.
        end select
      else
        if (iachar(value(i:i))==0) p =.false.
      endif
    enddo
  end function checkChars

  function isUSASCII(encoding) result(p)
    character(len=*), intent(in) :: encoding
    logical :: p

    character(len=len(encoding)) :: enc
    enc = toLower(encoding)
    p = (enc=="ansi_x3.4-1968" &
      .or. enc=="ansi_x3.4-1986" &
      .or. enc=="iso_646.irv:1991" &
      .or. enc=="ascii" &
      .or. enc=="iso646-us" &
      .or. enc=="us-ascii" &
      .or. enc=="us" &
      .or. enc=="ibm367" &
      .or. enc=="cp367" &
      .or. enc=="csascii")

  end function isUSASCII

  function allowed_encoding(encoding) result(p)
    character(len=*), intent(in) :: encoding
    logical :: p

    character(len=len(encoding)) :: enc
    logical :: utf8, usascii, iso88591, iso88592, iso88593, iso88594, &
      iso88595, iso88596, iso88597, iso88598, iso88599, iso885910,  &
      iso885913, iso885914, iso885915, iso885916

    enc = toLower(encoding)

    ! From http://www.iana.org/assignments/character-sets
    ! We can only reliably do US-ASCII (the below is mostly
    ! a list of synonyms for US-ASCII) but we also accept
    ! UTF-8 as a practicality. We bail out if any non-ASCII
    ! characters are used later on.
    utf8 = (enc=="utf-8")

    usascii = (enc=="ansi_x3.4-1968" &
      .or. enc=="ansi_x3.4-1986" &
      .or. enc=="iso_646.irv:1991" &
      .or. enc=="ascii" &
      .or. enc=="iso646-us" &
      .or. enc=="us-ascii" &
      .or. enc=="us" &
      .or. enc=="ibm367" &
      .or. enc=="cp367" &
      .or. enc=="csascii")
! As of FoX 4.0, we accept ISO-8859-??, also as practicality
! since we know it is identical to ASCII as far as 0x7F

    iso88591 = (enc =="iso_8859-1:1987" &
      .or. enc=="iso-ir-100" &
      .or. enc=="iso_8859-1" &
      .or. enc=="iso-8859-1" &
      .or. enc=="latin1" &
      .or. enc=="l1" &
      .or. enc=="ibm819" &
      .or. enc=="cp819" &
      .or. enc=="csisolatin1")

    iso88592 = (enc=="iso_8859-2:1987" &
      .or. enc=="iso-ir-101" &
      .or. enc=="iso_8859-2" &
      .or. enc=="iso-8859-2" &
      .or. enc=="latin2" &
      .or. enc=="l2" &
      .or. enc=="csisolatin2")

    iso88593 = (enc=="iso_8859-3:1988" &
      .or. enc=="iso-ir-109" &
      .or. enc=="iso_8859-3" &
      .or. enc=="iso-8859-3" &
      .or. enc=="latin3" &
      .or. enc=="l3" &
      .or. enc=="csisolatin3")

    iso88594 = (enc=="iso_8859-4:1988" &
      .or. enc=="iso-ir-110" &
      .or. enc=="iso_8859-4" &
      .or. enc=="iso-8859-4" &
      .or. enc=="latin4" &
      .or. enc=="l4" &
      .or. enc=="csisolatin4")

    iso88595 = (enc=="iso_8859-5:1988" &
      .or. enc=="iso-ir-144" &
      .or. enc=="iso_8859-5" &
      .or. enc=="iso-8859-5" &
      .or. enc=="cyrillic" &
      .or. enc=="csisolatincyrillic")

    iso88596 = (enc=="iso_8859-6:1987" &
      .or. enc=="iso-ir-127" &
      .or. enc=="iso_8859-6" &
      .or. enc=="iso-8859-6" &
      .or. enc=="ecma-114" &
      .or. enc=="asmo-708" &
      .or. enc=="arabic" &
      .or. enc=="csisolatinarabic")

    iso88597 = (enc=="iso_8859-7:1987" &
      .or. enc=="iso-ir-126" &
      .or. enc=="iso_8859-7" &
      .or. enc=="iso-8859-7" &
      .or. enc=="elot_928" &
      .or. enc=="ecma-118" &
      .or. enc=="greek" &
      .or. enc=="greek8" &
      .or. enc=="csisolatingreek")

    iso88598 = (enc=="iso_8859-8:1988" &
      .or. enc=="iso-ir-138" &
      .or. enc=="iso_8859-8" &
      .or. enc=="iso-8859-8" &
      .or. enc=="hebrew" &
      .or. enc=="csisolatinhebrew")

    iso88599 = (enc=="iso_8859-9:1989" &
      .or. enc=="iso-ir-148" &
      .or. enc=="iso_8859-9" &
      .or. enc=="iso-8859-9" &
      .or. enc=="latin5" &
      .or. enc=="l5" &
      .or. enc=="csisolatin5")

    iso885910 = (enc=="iso-8859-10" &
      .or. enc=="iso-ir-157" &
      .or. enc=="l6" &
      .or. enc=="iso_8859-10:1992" &
      .or. enc=="csisolatin6" &
      .or. enc=="latin6")

! ISO 6937 replaces $ sign with currency sign.
! JIS-X0201 has Yen instead of backslash, macron instead of tilde
! 16, 17, 18, 19 - Japanese encoding we can't use.
! BS 4730 replaces hash with UK pound sign, and tilde to macron
! 21, 22, 23, 24, 25, 26 - other variants of iso646, similar but not identical
      
!      iso10646utf1 = (enc=="iso-10646-utf-1") ! FIXME check
!      iso656basic1983 = (enc=="iso_646.basic:1983" &
!        .or. enc=="csiso646basic1983") ! FIXME check
! INVARIANT - almost but not quite a subset of ASCII
!      iso646irv = (enc=="iso_646.irv:1983" &
!        .or. enc=="iso-ir-2" &
!        .or. enc=="irv")
! 31, 32, 33, 34 - NATS scandinavian, different from ASCII
! 35 - another iso646 variant
! 36, 37, 38 Korean shifted/multibyte
! 39, 40 Japanese shifted/multibyte
! 41, 42, JIS (iso646inv 7 bits)
! 43 another iso646 variantt
! 44, 45 greek variants
! 46 another iso646 variant
! 47 greek
! 48 cyrillic ascii relationship unknown
! 49 JIS again
! 50 similar not identical
! 51, 52, 53 not identical
! 54 see 48
! 55 see 47
! 56 another iso646 variant
! 57 chinese
! ... to be continued

      iso885913 = (enc=="iso-8859-13")
      iso885914 = (enc=="iso-8859-14" &
      .or. enc=="iso-ir-199" &
      .or. enc=="iso_8859-14:1998" &
      .or. enc=="iso_8849-14" &
      .or. enc=="iso_latin8" &
      .or. enc=="iso-celtic" &
      .or. enc=="l8")
      iso885915 = (enc=="iso-8859-15" &
        .or. enc=="iso-8859-15" &
        .or. enc=="latin-9")
      iso885916 = (enc=="iso-8859-16" &
        .or. enc=="iso-ir226" &
        .or. enc=="iso_8859-16:2001" &
        .or. enc=="iso_8859-16" &
        .or. enc=="latin10" &
        .or. enc=="l10")

      p = utf8.or.usascii.or.iso88591.or.iso88592.or.iso88593 &
        .or.iso88594.or.iso88595.or.iso88596.or.iso88597 &
        .or.iso88598.or.iso88599.or.iso885910.or.iso885913 &
        .or.iso885914.or.iso885915.or.iso885916

  end function allowed_encoding

#endif
end module m_common_charset
