module fox_m_utils_uri
#ifndef DUMMYLIB

  ! Manipulate URIs and URI references a la RFC 2396
  ! NB: ...
  ! Forbidden (ASCII control) characters are not handled correctly
  ! checking of reg names (not hosts) is done wrongly
  ! checking of ipv6/X is untested

  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc, vs_vs_alloc
  use fox_m_fsys_format, only: str_to_int_10, str_to_int_16, str
  use fox_m_fsys_string, only: toLower

  implicit none
  private

  type path_segment
    character, pointer :: s(:) => null()
  end type path_segment
#endif

  type URI
    private
#ifndef DUMMYLIB
    character, pointer :: scheme(:) => null()
    character, pointer :: authority(:) => null()
    character, pointer :: userinfo(:) => null()
    character, pointer :: host(:) => null()
    integer :: port = -1
    character, pointer :: path(:) => null()
    type(path_segment), pointer :: segments(:) => null()
    character, pointer :: query(:) => null()
    character, pointer :: fragment(:) => null()
#else
    integer :: i
#endif
  end type URI

#ifndef DUMMYLIB
  character(len=*), parameter :: lowalpha = "abcdefghijklmnopqrstuvwxyz"
  character(len=*), parameter :: upalpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(len=*), parameter :: alpha = lowalpha//upalpha
  character(len=*), parameter :: digit = "0123456789"
  character(len=*), parameter :: hexdigit = "0123456789abcdefABCDEF"
  character(len=*), parameter :: alphanum = alpha//digit
  character(len=*), parameter :: unreserved = alphanum//"-._~"
  character(len=*), parameter :: gen_delims  = ":/?#[]@"
  character(len=*), parameter :: sub_delims  = "!$&'()*+,;="
  character(len=*), parameter :: reserved = gen_delims//sub_delims
  character(len=*), parameter :: pchar = unreserved//":@&=+$,"
  character(len=*), parameter :: uric_no_slash = unreserved//";?:@&=+$,"
  character(len=*), parameter :: uric = unreserved//reserved
  character(len=*), parameter :: unwise = "{}|\^[]`"
#endif 

  public :: URI
  public :: parseURI
  public :: expressURI
  public :: isAbsoluteURI
  public :: rebaseURI
  public :: copyURI
  public :: destroyURI

  public :: hasScheme
  public :: getScheme
  public :: hasAuthority
  public :: getAuthority
  public :: hasUserinfo
  public :: getUserinfo
  public :: hasHost
  public :: getHost
  public :: hasPort
  public :: getPort
  public :: getPath
  public :: hasQuery
  public :: getQuery
  public :: hasFragment
  public :: getFragment

#ifndef DUMMYLIB
  public :: dumpURI
#endif

contains

#ifndef DUMMYLIB
  function unEscape_alloc(s) result(c)
    character(len=*), intent(in) :: s
    character, pointer :: c(:)

    integer :: i, j, n
    character(len(s)) :: t

    c => null()

    i = 1
    j = 0
    do while (i<=len(s))
      j = j + 1
      if (s(i:i)=="%") then
        if (i+2>len(s)) return
        if (verify(s(i+1:i+2), hexdigit)/=0) return
        n = str_to_int_16(s(i+1:i+2))
        t(j:j) = achar(n)
        i = i + 3
      else
        t(j:j) = s(i:i)
        i = i + 1
      endif
    enddo

    c => vs_str_alloc(t(:j))
  end function unEscape_alloc

  function verifyWithPctEncoding(s, chars) result(p)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: chars
    logical :: p

    integer :: i

    p = .false.
    i = 1
    do while (i<=len(s))
      if (s(i:i)=="%") then
        if (i+2>len(s)) return
        if (verify(s(i+1:i+2), hexdigit)>0) return
        i = i + 3
      else
        if (verify(s(i:i),chars)>0) return
        i = i + 1
      endif
    enddo
    p = .true.
  end function verifyWithPctEncoding

  pure function pctEncode_len(s, chars) result(n)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: chars
    integer :: n

    integer :: i
    n = 0
    do i = 1, len(s)
      n = n + 1
      if (verify(s(i:i), unwise)==0.or.verify(s(i:i), chars)>0) n = n + 2
    enddo

  end function pctEncode_len

  function pctEncode(s, chars) result(ps)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: chars
    character(len=pctEncode_len(s, chars)) :: ps

    integer :: i, n

    n = 1
    do i = 1, len(s)
      if (verify(s(i:i), unwise)==0.or.verify(s(i:i), chars)>0) then
        ps(n:n+2) = "%"//str(iachar(s(i:i)), "x2")
        n = n + 3
      else
        ps(n:n) = s(i:i)
        n = n + 1
      endif
    enddo

  end function pctEncode

  function checkOpaquePart(part) result(p)
    character(len=*), intent(in) :: part
    logical :: p

    if (len(part)>0) then
      p = verify(part(1:1), uric_no_slash)==0
      if (p.and.len(part)>1) &
        p = verify(part(1:1), uric)==0
    endif
  end function checkOpaquePart
    

  function checkScheme(scheme) result(p)
    character(len=*), intent(in) :: scheme
    logical :: p

    p = len(scheme)>0 
    if (p) then
      p = verify(scheme(1:1), lowalpha//upalpha)==0
      if (p.and.len(scheme)>1) then
        p = verify(scheme(2:), alphanum//"+-.")==0
      endif
    endif
  end function checkScheme

  function checkIpvX(host) result(p)
    character(len=*), intent(in) :: host
    logical :: p

    integer :: i, n1, n2

    p = (len(host)>5).and.(host(1:1)=="[".and.host(len(host):len(host))=="]")

    if (p) then

      ! Try IPvFuture:
      p = (verify(host(2:2),"Vv")==0 &
        .and.verify(host(3:3),hexdigit)==0 &
        .and.host(4:4)=="." &
        .and.verify(host(3:3),unreserved//sub_delims//":")==0)

      if (.not.p) then ! is it IPv6?
        n1 = 0
        do i = 1, 4
          n2 = index(host(n1+1:), ":")
          if (n2==0.or.n2>6) return
          n2 = n2 + n1
          if (verify(host(n1+1:n2-1),hexdigit)>0) return
          n1 = n2
        enddo
        n2 = index(host(n1+1:), ":")
        if (n2==0) then
          ! this must be ipv4 format
          do i = 1, 3
            n2 = index(host(n1+1:), ".")
            if (n2==0) return
            n2 = n2 + n1
            if (verify(host(n1+1:n2-1),digit)>0) return
            if (str_to_int_10(host(n1+1:n2-1))>255) return
            n1 = n2
          enddo
          ! Now there must be 3 or less digits followed by ]
          n2 = len(host)-1
          if (verify(host(n1+1:n2-1),digit)>0) return
          if (str_to_int_10(host(n1+1:n2-1))>255) return 
        elseif (n2<6) then
          n2 = n2 + n1
          if (verify(host(n1+1:n2-1),hexdigit)>0) return
          ! Now there must be 4 or less digits followed by ]
          n1 = n2
          n2 = len(host)
          if (n2-n1>4) return
          if (verify(host(n1+1:n2-1),hexdigit)>0) return
        endif
        p = .true.
      endif
    endif
  end function checkIpvX

  function checkHost(host) result(p)
    character(len=*), intent(in) :: host
    logical :: p

    p = checkIpvX(host)
    if (.not.p) &
      p = verifyWithPctEncoding(host, unreserved//sub_delims)

  end function checkHost

  function checkAuthority(authority, userinfo, host, port) result(p)
    character(len=*), intent(in) :: authority
    character, pointer :: userinfo(:), host(:)
    integer :: port
    logical :: p

    integer :: i1, i2

    p = .true.
    if (len(authority)==0) return

    i1 = index(authority, "@")
    if (i1>0) then
      i2 = index(authority(i1+1:), ":")
    else
      i2 = index(authority, ":")
    endif
    if (i1==0) then
      userinfo => null()
    else
      p = verifyWithPctEncoding(authority(:i1-1), unreserved//sub_delims//":")
      if (p) userinfo => unEscape_alloc(authority(:i1-1))
    endif
    if (i2==0) then
      i2 = len(authority)+1
    else
      i2 = i1 + i2
      p = p.and.verify(authority(i2+1:), digit)==0
      if (p) port = str_to_int_10(authority(i2+1:))
    endif
    p = p.and.checkHost(authority(i1+1:i2-1))
    if (p) then
      host => vs_str_alloc(authority(i1+1:i2-1))
    else
      if (associated(userinfo)) deallocate(userinfo)
    end if

  end function checkAuthority
  
  function checkPathSegment(segment) result(p)
    character(len=*), intent(in) :: segment
    logical :: p

    integer :: i1

    i1 = index(segment, ";")
    if (i1>0) then
      p = verifyWithPctEncoding(segment(:i1-1), pchar) &
        .and.verifyWithPctEncoding(segment(i1+1:), pchar)
    else
      p = verifyWithPctEncoding(segment, unreserved//pchar)
    endif
  end function checkPathSegment

  function checkNonOpaquePath(path, segments) result(p)
    character(len=*), intent(in) :: path
    type(path_segment), pointer :: segments(:)
    logical :: p

    integer :: i, i1, i2
    type(path_segment), pointer :: temp(:)
    
    p = .true.
    i1 = index(path, "/")
    if (i1==1) then
      allocate(segments(1))
      segments(1)%s => vs_str_alloc("/")
    else
      allocate(segments(0))
      i1 = 0
    endif

    do
      i2 = index(path(i1+1:), "/")
      if (i2==0) then
        i2 = len(path)
      else
        i2 = i1 + i2
      endif
      if (checkPathSegment(path(i1+1:i2-1))) then
        allocate(temp(size(segments)+1))
        do i = 1, size(segments)
          temp(i)%s => segments(i)%s
        enddo
        temp(i)%s => unEscape_alloc(path(i1+1:i2))
        deallocate(segments)
        segments => temp
      else
        do i = 1, size(segments)
          deallocate(segments(i)%s)
        enddo
        deallocate(segments)
        p = .false.
        return
      endif
      if (i2==len(path)) exit
      i1 = i2
    end do
  end function checkNonOpaquePath

  function checkPath(path, segments) result(p)
    character(len=*), intent(in) :: path
    type(path_segment), pointer :: segments(:)
    logical :: p

    p = checkNonOpaquePath(path, segments)
    if (.not.p) then
      p = checkOpaquePart(path)
      if (p) allocate(segments(0))
    endif

  end function checkPath

  function checkQuery(query) result(p)
    character(len=*), intent(in) :: query
    logical :: p

    p = verifyWithPctEncoding(query, uric)
  end function checkQuery

  function checkFragment(fragment) result(p)
    character(len=*), intent(in) :: fragment
    logical :: p

    p = verifyWithPctEncoding(fragment, uric)
  end function checkFragment
#endif

  function parseURI(URIstring) result(u)
    character(len=*), intent(in) :: URIstring
    type(URI), pointer :: u
#ifndef DUMMYLIB
    character, pointer, dimension(:) :: scheme, authority, &
      userinfo, host, path, query, fragment
    integer :: port
    type(path_segment), pointer :: segments(:)
    integer :: i1, i2, i3, i4
    logical :: p

#endif
    u => null()
#ifndef DUMMYLIB

    scheme => null()
    authority => null()
    userinfo => null()
    host => null()
    port = -1
    path => null()
    segments => null()
    query => null()
    fragment => null()

    i1 = index(URIstring, ":")
    if (i1>0) then 
      p = checkScheme(URIstring(:i1-1))
      if (p) then
        scheme => vs_str_alloc(toLower(URIstring(:i1-1)))
      else
        i1 = 0
      endif
    endif
    ! if either i1==0 or the scheme doesn't validate, there is no scheme..
    if (len(URIstring)>=i1+3) then
      if (URIstring(i1+1:i1+2)=="//") then
        i2 = scan(URIstring(i1+3:), "/#?")
        if (i2==0) then
          i2 = len(URIstring) + 1
        else
          i2 = i1 + i2 + 2
        endif
        p = checkAuthority(URIstring(i1+3:i2-1), userinfo, host, port)
        if (.not.p) then
          call cleanUp
          return
        endif
        authority => vs_str_alloc(URIstring(i1+3:i2-1))
      else
        i2 = i1 + 1
      endif
    else
      i2 = i1 + 1
    endif

    if (i2>len(URIstring)) then
      path => vs_str_alloc("")
      allocate(segments(1))
      segments(1)%s => vs_str_alloc("")
      call produceResult
      return
    endif

    i3 = scan(URIstring(i2:),"#?")
    if (i3==0) then
      i3 = len(URIstring) + 1
    else
      i3 = i2 + i3 - 1
    endif
    p = checkPath(URIstring(i2:i3-1), segments)
    if (.not.p) then
      call cleanUp
      return
    endif
    path => unEscape_alloc(URIstring(i2:i3-1))

    if (i3>len(URIstring)) then
      call produceResult
      return
    endif

    if (URIstring(i3:i3)=="?") then
      i4 = index(URIstring(i3+1:), "#")
      if (i4==0) then
        i4 = len(URIstring) + 1
      else
        i4 = i3 + i4
      endif
      p = checkQuery(URIstring(i3+1:i4-1))
      if (.not.p) then
        call cleanUp
        return
      endif
      query => vs_str_alloc(URIstring(i3+1:i4-1))
    else
      i4 = i3
    endif

    if (i4>len(URIstring)) then
      call produceResult
      return
    endif

    p = checkFragment(URIstring(i4+1:))
    if (.not.p) then
      call cleanUp
      return
    endif
    fragment => vs_str_alloc(URIstring(i4+1:))
    call produceResult

    contains
      subroutine cleanUp
        integer :: i
        if (associated(scheme)) deallocate(scheme)
        if (associated(authority)) deallocate(authority)
        if (associated(userinfo)) deallocate(userinfo)
        if (associated(host)) deallocate(host)
        if (associated(path)) deallocate(path)
        if (associated(query)) deallocate(query)
        if (associated(fragment)) deallocate(fragment)
        if (associated(segments)) then
          do i = 1, size(segments)
            deallocate(segments(i)%s)
          enddo
          deallocate(segments)
        endif
      end subroutine cleanUp
      subroutine produceResult
        allocate(u)
        u%scheme => scheme
        u%authority => authority
        u%userinfo => userinfo
        u%host => host
        u%port = port
        u%path => path
        u%segments => segments
        u%query => query
        u%fragment => fragment
        u%segments => segments
      end subroutine produceResult
#endif
  end function parseURI

  function isAbsoluteURI(u) result(p)
    type(URI), intent(in) :: u
    logical :: p

#ifdef DUMMYLIB
    p = .false.
#else
    p = associated(u%scheme).or.associated(u%authority)
    if (.not.p.and.size(u%segments(1)%s)>0) then
      p = u%segments(1)%s(1)=="/"
    endif
#endif
  end function isAbsoluteURI

  function rebaseURI(u1, u2) result(u3)
    type(URI), pointer :: u1, u2
    type(URI), pointer :: u3

    u3 => null()
#ifndef DUMMYLIB

    if (associated(u2%scheme).or.associated(u2%authority)) then
      u3 => copyURI(u2)
      return
    endif

    allocate(u3)
    if (associated(u1%scheme)) u3%scheme => vs_vs_alloc(u1%scheme)
    if (associated(u1%authority)) u3%authority => vs_vs_alloc(u1%authority)

    u3%segments => appendPaths(u1%segments, u2%segments)
    u3%path => expressSegments(u3%segments)

    if (associated(u2%query)) u3%query => vs_vs_alloc(u2%query)
    if (associated(u2%fragment)) u3%fragment => vs_vs_alloc(u2%fragment)
#endif
  end function rebaseURI

#ifndef DUMMYLIB
  function appendPaths(seg1, seg2) result(seg3)
    type(path_segment), pointer :: seg1(:), seg2(:)
    type(path_segment), pointer :: seg3(:)

    type(path_segment), pointer :: temp(:)

    integer :: i, n, n2

    if (size(seg2(1)%s)==0) then
      seg3 => normalizePath(seg1)
      return
    elseif (seg2(1)%s(1)=="/") then
      seg3 => normalizePath(seg2)
      return
    endif

    n = size(seg1) + size(seg2)
    i = size(seg1)
    if (seg1(i)%s(size(seg1(i)%s))/="/") &
      n = n - 1

    allocate(temp(n))
    n2 = 1
    do i = 1, size(seg1)
      if (i==size(seg1).and.seg1(i)%s(size(seg1(i)%s))/="/") exit ! it's a file
      temp(n2)%s => vs_vs_alloc(seg1(i)%s)
      n2 = n2 + 1
    enddo
      
    do i = 1, size(seg2)
      temp(n2)%s => vs_vs_alloc(seg2(i)%s)
      n2 = n2 + 1
    enddo

    seg3 => normalizePath(temp)
    do i = 1, size(temp)
      deallocate(temp(i)%s)
    enddo
    deallocate(temp)

  end function appendPaths

  function normalizepath(seg1) result(seg2)
    type(path_segment), pointer :: seg1(:)
    type(path_segment), pointer :: seg2(:)

    integer :: i, n, n2, parents
    character, pointer :: tmp(:) 
    

    ! If the last of the input segments are
    ! equal to '.' or '..', append a slash
    ! so the rest of the subroutine works.

    if ((str_vs(seg1(size(seg1))%s) == '.').or. &
        (str_vs(seg1(size(seg1))%s) == '..')) then
        tmp => vs_vs_alloc(seg1(size(seg1))%s)
        deallocate(seg1(size(seg1))%s)
        seg1(size(seg1))%s => vs_str_alloc(str_vs(tmp)//"/")
        deallocate(tmp)
    endif

    n = 0
    parents = 0
    do i = 1, size(seg1)
      if (str_vs(seg1(i)%s)//"x"=="./x") then
        continue
      elseif (str_vs(seg1(i)%s)//"x"=="../x") then
        if (n>0) then
          n = n - 1
        else
          parents = parents + 1
        endif
      else
        n = n + 1
      endif
    enddo

    n = n + parents
    allocate(seg2(n))

    n2 = parents
    do i = 1, parents
      seg2(i)%s => vs_str_alloc("../")
    enddo
    do i = 1, size(seg1)
      if (str_vs(seg1(i)%s)//"x"=="./x") then
        continue
      elseif (str_vs(seg1(i)%s)//"x"=="../x") then
        if (n2>parents) then
          if (n2<=n) deallocate(seg2(n2)%s)
          n2 = n2 - 1
        endif
      else
        n2 = n2 + 1
        if (n2>0.and.n2<=n) &
          seg2(n2)%s => vs_vs_alloc(seg1(i)%s)
      endif
    enddo

  end function normalizepath

  function expressSegments(seg1) result(s)
    type(path_segment), pointer :: seg1(:)
    character, pointer :: s(:)

    integer :: i, n

    n = 0

    do i = 1, size(seg1)
      n = n + size(seg1(i)%s)
    enddo
    allocate(s(n))
    n = 1
    do i = 1, size(seg1)
      s(n:n+size(seg1(i)%s)-1) = seg1(i)%s
      n = n + size(seg1(i)%s)
    enddo
  end function expressSegments

  pure function expressURI_len(u) result(n)
    type(URI), intent(in) :: u
    integer :: n

    n = 0
    if (associated(u%scheme)) &
      n = size(u%scheme) + 1
    if (associated(u%authority)) &
      n = n + pctEncode_len(str_vs(u%authority), unreserved//sub_delims//"@:") + 2
    !FIXME - I suspect that ';' as the first character of a segment should be escaped
    n = n + pctEncode_len(str_vs(u%path), pchar//";"//"/")
    if (associated(u%query)) &
      n = n + pctEncode_len(str_vs(u%query), uric) + 1
    if (associated(u%fragment)) &
      n = n + pctEncode_len(str_vs(u%fragment), uric) + 1

  end function expressURI_len

  function expressURI(u) result(URIstring)
    type(URI), intent(in) :: u
    character(len=expressURI_len(u)) :: URIstring

    integer :: i, j
    URIstring=""
    i = 1
    if (associated(u%scheme)) then
      URIstring(:size(u%scheme)+1) = str_vs(u%scheme)//":"
      i = i + size(u%scheme) + 1
    endif
    if (associated(u%authority)) then
      j = pctEncode_len(str_vs(u%authority), unreserved//sub_delims//"@:")
      URIstring(i:i+j+1) = &
        "//"//pctEncode(str_vs(u%authority), unreserved//sub_delims//"@:")
      i = i + j + 2
    endif
    if (size(u%path)>0) then
      !FIXME - I suspect that ';' as the first character of a segment should be escaped
      j = pctEncode_len(str_vs(u%path), pchar//";"//"/")
      URIstring(i:i+j-1) = pctEncode(str_vs(u%path), pchar//";"//"/")
      i = i + j
    endif
    if (associated(u%query)) then
      j = pctEncode_len(str_vs(u%query), uric)
      URIstring(i:i+j) = "?"//pctEncode(str_vs(u%query), uric)
      i = i + j + 1
    endif
    if (associated(u%fragment)) then
      j = pctEncode_len(str_vs(u%fragment), uric)
      URIstring(i:i+j) = "#"//pctEncode(str_vs(u%fragment), uric)
    endif

  end function expressURI

  subroutine dumpURI(u)
    type(URI), intent(in) :: u
    integer :: i
    if (associated(u%scheme)) then
      write(*,*) "scheme: ", str_vs(u%scheme)
    else
      write(*,*) "scheme UNDEFINED"
    endif
    if (associated(u%authority)) then
      write(*,*) "authority: ", str_vs(u%authority)
    else
      write(*,*) "authority UNDEFINED"
    endif
    if (associated(u%userinfo)) then
      write(*,*) "userinfo: ", str_vs(u%userinfo)
    else
      write(*,*) "userinfo UNDEFINED"
    endif
    if (associated(u%host)) then
      write(*,*) "host: ", str_vs(u%host)
    else
      write(*,*) "host UNDEFINED"
    endif
    if (u%port>0) then
      write(*,*) "port: ", str(u%port)
    else
      write(*,*) "port UNDEFINED"
    endif
    if (associated(u%path)) then
      write(*,*) "path: ", str_vs(u%path)
    else
      write(*,*) "path UNDEFINED"
    endif
    if (associated(u%segments)) then
      do i = 1, size(u%segments)
        write(*,*) "    segment: ", str_vs(u%segments(i)%s)
      enddo
    endif
    if (associated(u%query)) then
      write(*,*) "query: ", str_vs(u%query)
    else
      write(*,*) "query UNDEFINED"
    endif
    if (associated(u%fragment)) then
      write(*,*) "fragment: ", str_vs(u%fragment)
    else
      write(*,*) "fragment UNDEFINED"
    endif
  end subroutine dumpURI
#endif

  function copyURI(u1) result(u2)
    type(URI), pointer :: u1
    type(URI), pointer :: u2
#ifndef DUMMYLIB
    integer :: i

    if (.not.associated(u1)) then
#endif
      u2 => null()
#ifndef DUMMYLIB
      return
    endif
    allocate(u2)
    u2%scheme => vs_vs_alloc(u1%scheme)
    u2%authority => vs_vs_alloc(u1%authority)
    u2%userinfo => vs_vs_alloc(u1%userinfo)
    u2%host => vs_vs_alloc(u1%host)
    u2%port = u1%port
    u2%path => vs_vs_alloc(u1%path)
    allocate(u2%segments(size(u1%segments)))
    do i = 1, size(u1%segments)
      u2%segments(i)%s => vs_vs_alloc(u1%segments(i)%s)
    enddo
    u2%query => vs_vs_alloc(u1%query)
    u2%fragment => vs_vs_alloc(u1%fragment)
#endif
  end function copyURI


  subroutine destroyURI(u)
    type(URI), pointer :: u
#ifndef DUMMYLIB
    integer :: i
    if (associated(u%scheme)) deallocate(u%scheme)
    if (associated(u%authority)) deallocate(u%authority)
    if (associated(u%userinfo)) deallocate(u%userinfo)
    if (associated(u%host)) deallocate(u%host)
    if (associated(u%path)) deallocate(u%path)
    if (associated(u%segments)) then
      do i = 1, size(u%segments)
        deallocate(u%segments(i)%s)
      enddo
      deallocate(u%segments)
    endif
    if (associated(u%query)) deallocate(u%query)
    if (associated(u%fragment)) deallocate(u%fragment)

    deallocate(u)
#endif
  end subroutine destroyURI

  function hasScheme(u) result(p)
    type(URI), pointer :: u
    logical :: p

    p = .false.
#ifndef DUMMYLIB
    if (.not.associated(u)) return
    p = associated(u%scheme)
#endif
  end function hasScheme

  function getScheme(u) result(s)
    type(URI), pointer :: u

#ifndef DUMMYLIB
    character(len=size(u%scheme)) :: s
    s = str_vs(u%scheme)
#else
    character(len=1) :: s
    s = ""
#endif
  end function getScheme

  function hasAuthority(u) result(p)
    type(URI), pointer :: u
    logical :: p

    p = .false.
#ifndef DUMMYLIB
    if (.not.associated(u)) return
    p = associated(u%authority)
#endif
  end function hasAuthority

  function getAuthority(u) result(s)
    type(URI), pointer :: u
#ifndef DUMMYLIB
    character(len=size(u%authority)) :: s
    s = str_vs(u%authority)
#else
    character(len=1) :: s
    s = ""
#endif
  end function getAuthority

  function hasUserinfo(u) result(p)
    type(URI), pointer :: u
    logical :: p

    p = .false.
#ifndef DUMMYLIB
    if (.not.associated(u)) return
    p = associated(u%userinfo)
#endif
  end function hasUserinfo

  function getUserinfo(u) result(s)
    type(URI), pointer :: u
#ifndef DUMMYLIB
    character(len=size(u%userinfo)) :: s
    s = str_vs(u%userinfo)
#else
    character(len=1) :: s
    s = ""
#endif
  end function getUserinfo

  function hasHost(u) result(p)
    type(URI), pointer :: u
    logical :: p

    p = .false.
#ifndef DUMMYLIB
    if (.not.associated(u)) return
    p = associated(u%host)
#endif
  end function hasHost

  function getHost(u) result(s)
    type(URI), pointer :: u
#ifndef DUMMYLIB
    character(len=size(u%host)) :: s
    s = str_vs(u%host)
#else
    character(len=1) :: s
    s = ""
#endif
  end function getHost

  function hasPort(u) result(p)
    type(URI), pointer :: u
    logical :: p

    p = .false.
#ifndef DUMMYLIB
    if (.not.associated(u)) return
    p = u%port > 0
#endif
  end function hasPort

  function getPort(u) result(n)
    type(URI), pointer :: u
    integer :: n
#ifndef DUMMYLIB
    n = u%port
#else
    n = 0
#endif
  end function getPort

  function getPath(u) result(s)
    type(URI), pointer :: u
#ifndef DUMMYLIB
    character(len=size(u%path)) :: s
    s = str_vs(u%path)
#else
    character(len=1) :: s
    s = ""
#endif
  end function getPath

  function hasQuery(u) result(p)
    type(URI), pointer :: u
    logical :: p

    p = .false.
#ifndef DUMMYLIB
    if (.not.associated(u)) return
    p = associated(u%query)
#endif
  end function hasQuery

  function getQuery(u) result(s)
    type(URI), pointer :: u
#ifndef DUMMYLIB
    character(len=size(u%query)) :: s
    s = str_vs(u%query)
#else
    character(len=1) :: s
    s = ""
#endif
  end function getQuery

  function hasFragment(u) result(p)
    type(URI), pointer :: u
    logical :: p

    p = .false.
#ifndef DUMMYLIB
    if (.not.associated(u)) return
    p = associated(u%fragment)
#endif
  end function hasFragment

  function getFragment(u) result(s)
    type(URI), pointer :: u
#ifndef DUMMYLIB
    character(len=size(u%fragment)) :: s
    s = str_vs(u%fragment)
#else
    character(len=1) :: s
    s = ""
#endif
  end function getFragment

end module fox_m_utils_uri
