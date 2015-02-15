program short_test

  use fox_m_fsys_realtypes, only: sp, dp
  use FoX_common, only : rts, countrts

  print*, "test.0.1.1"
  call stringdatascalar("abcd", "abcd", 1, 0)
  print*, "test.0.1.2"
  call stringdatascalar(" abcd ", "abcd", 1, 0)
  print*, "test.0.1.3"
  call stringdatascalar(" abc d ", "abc d", 1, 0)
  print*, "test.0.1.4"
  call stringdatascalar(" abc  d ", "abc d", 1, 0)
  print*, "test.0.1.5"
  call stringdatascalar(" abc "//achar(13)//"  d ", "abc d", 1, 0)

  print*, "test.0.2.1"
  call logicaldatascalar("true", .true., 1, 0)
  print*, "test.0.2.2"
  call logicaldatascalar("false", .false., 1, 0)
  print*, "test.0.2.3"
  call logicaldatascalar("1", .true., 1, 0)
  print*, "test.0.2.4"
  call logicaldatascalar("0", .false., 1, 0)
  print*, "test.0.2.5"
  call logicaldatascalar("  true  ", .true., 1, 0)
  print*, "test.0.2.6"
  call logicaldatascalar("  false  ", .false., 1, 0)
  print*, "test.0.2.7"
  call logicaldatascalar("  1 ", .true., 1, 0)
  print*, "test.0.2.8"
  call logicaldatascalar("  0 ", .false., 1, 0)
  print*, "test.0.2.9"
  call logicaldatascalar("  .tru ", .false., 0, 2)
  print*, "test.0.2.10"
  call logicaldatascalar("  false. ", .false., 0, 2)

  print*, "test.0.3.1"
  call integerdatascalar(" 1 ", 1, 1, 0)
  print*, "test.0.3.2"
  call integerdatascalar(" -10 ", -10, 1, 0)
  print*, "test.0.3.3"
  call integerdatascalar(" -1.3 ", 0, 0, 2)

  print*, "test.0.4.1"
  call realspdatascalar(" 1 ", 1.0, 1, 0)
  print*, "test.0.4.2"
  call realspdatascalar(" -1.35D2", -1.35e2, 1, 0)
  print*, "test.0.4.3"
  call realspdatascalar(" -1.35FG", 0.0, 0, 2)

  print*, "test.0.5.1"
  call realdpdatascalar(" 1 ", 1.0_dp, 1, 0)
  print*, "test.0.5.2"
  call realdpdatascalar(" -1.35D2", -1.35e2_dp, 1, 0)
  print*, "test.0.5.3"
  call realdpdatascalar(" -1.35FG", 0.0_dp, 0, 2)

  print*, "test.0.6.1"
  call complexspdatascalar(" 1,2 ", (1.0,2.0), 1, 0)
  print*, "test.0.6.2"
  call complexspdatascalar(" (1)+i(2)", (1.0,2.0), 1, 0)
  print*, "test.0.6.3"
  call complexspdatascalar(" arse", (0.0,0.0), 0, 2)

  print*, "test.0.7.1"
  call complexdpdatascalar(" 1,2 ", (1.0_dp,2.0_dp), 1, 0)
  print*, "test.0.7.2"
  call complexdpdatascalar(" (1)+i(2)", (1.0_dp,2.0_dp), 1, 0)
  print*, "test.0.7.3"
  call complexdpdatascalar(" arse", (0.0_dp,0.0_dp), 0, 2)

  print*,"test.1.1.1"
  call stringdataarray("a b c d e", (/"a", "b", "c", "d", "e"/), 5, 0)
  print*,"test.1.1.2"
  call stringdataarray("a b c d", (/"a", "b", "c", "d", " "/), 4, -1)
  print*,"test.1.1.3"
  call stringdataarray("a b c d e f", (/"a", "b", "c", "d", "e"/), 5, 1)

  print*,"test.1.2.1"
  call stringdataarray("a, b, c, d, e", (/"a ", " b", " c", " d", " e"/), &
    5, 0, ",")
  print*,"test.1.2.2"
  call stringdataarray("a, b, c, d", (/"a ", " b", " c", " d", "  "/), &
    4, -1, ",")
  print*,"test.1.2.3"
  call stringdataarray("a, b, c, d, e, f", (/"a ", " b", " c", " d", " e"/), &
    5, 1, ",")
  print*,"test.1.2.4"
  call stringdataarray("a, b, c, d, e, ", (/"a ", " b", " c", " d", " e"/), &
    5, 1, ",")
  print*,"test.1.2.5"
  call stringdataarray("a, b, c, d, e,", (/"a ", " b", " c", " d", " e"/), &
    5, 0, ",")

  print*,"test.1.2.6"
  call stringdataarray("a, b, c, d, e", (/"a ", " b", " c", " d", " e"/), &
    5, 0, csv=.true.)
  print*,"test.1.2.7"
  call stringdataarray("a, b, c, d", (/"a ", " b", " c", " d", "  "/), &
    4, -1, csv=.true.)
  print*,"test.1.2.8"
  call stringdataarray("a, b, c, d, e, f", (/"a ", " b", " c", " d", " e"/), &
    5, 1, csv=.true.)
  print*,"test.1.2.9" 
  call stringdataarray("a, b, c, d, e, ", (/"a ", " b", " c", " d", " e"/), &
    5, 1, csv=.true.)
  print*,"test.1.2.10"
  call stringdataarray("a, b, c, d, e,", (/"a ", " b", " c", " d", " e"/), &
    5, 1, csv=.true.)
  print*,"test.1.2.11"
  call stringdataarray("a, b, c,, e", (/"a ", " b", " c", "  ", " e"/), & 
    5, 0, csv=.true.)
  print*,"test.1.2.12"
  call stringdataarray("a, b, c, d,", (/"a ", " b", " c", " d", "  "/), &  
    5, 0, csv=.true.)
  print*,"test.1.2.13"
  call stringdataarray(", b, c, d, e", (/"  ", " b", " c", " d", " e"/), & 
    5, 0, csv=.true.)

  print*,"test.1.3.1"
  call stringdatamatrix("a b c d e f", &
    reshape((/"a", "b", "c", "d", "e", "f"/), (/2,3/)), &
    6, 0)
  print*,"test.1.3.2"
  call stringdatamatrix("a b c d", &
    reshape((/"a", "b", "c", "d", " ", " "/), (/2,3/)), &
    4, -1)
  print*,"test.1.3.3"
  call stringdatamatrix("a b c d e f g", &
    reshape((/"a", "b", "c", "d", "e", "f"/), (/2,3/)), &
    6, 1)

  print*,"test.1.4.1"
  call stringdatamatrix("a, b, c, d, e, f", &
    reshape((/"a ", " b", " c", " d", " e", " f"/), (/2,3/)), &
    6, 0, ",")
  print*,"test.1.4.2"
  call stringdatamatrix("a, b, c, d", &
    reshape((/"a ", " b", " c", " d", "  ", "  "/), (/2,3/)), &
    4, -1, ",")
  print*,"test.1.4.3"
  call stringdatamatrix("a, b, c, d, e, f, g", &
    reshape((/"a ", " b", " c", " d", " e", " f"/), (/2,3/)), &
    6, 1, ",")
  print*,"test.1.4.4"
  call stringdatamatrix("a, b, c, d, e, f, ", &
    reshape((/"a ", " b", " c", " d", " e", " f"/), (/2,3/)), &
    6, 1, ",")
  print*,"test.1.4.5"
  call stringdatamatrix("a, b, c, d, e, f,", &
    reshape((/"a ", " b", " c", " d", " e", " f"/), (/2,3/)), &
    6, 0, ",")

  print*, "test.2.1.1"
  call logicaldataarray("true false false false true", &
    (/.true., .false., .false., .false., .true./), 5, 0)
  print*, "test.2.1.2"
  call logicaldataarray("true false false false", &
    (/.true., .false., .false., .false., .false./), 4, -1)
  print*, "test.2.1.3"
  call logicaldataarray("true false false false true true", &
    (/.true., .false., .false., .false., .true./), 5, 1)
  print*, "test.2.1.4"
  call logicaldataarray("true, false, false, false true", &
    (/.true., .false., .false., .false., .true./), 5, 0)
  print*, "test.2.1.5"
  call logicaldataarray("true, false, false, ,false true", &
    (/.true., .false., .false., .false., .false./), 3, 2)

  print*, "test.2.2.1"
  call logicaldatamatrix("true false false false true true", &
    reshape((/.true., .false., .false., .false., .true., .true./), (/2,3/)), &
    6, 0)
  print*, "test.2.2.2"
  call logicaldatamatrix("true false false false", &
    reshape((/.true., .false., .false., .false., .false., .false./), (/2,3/)), &
    4, -1)
  print*, "test.2.2.3"
  call logicaldatamatrix("true false false false true true false", &
    reshape((/.true., .false., .false., .false., .true., .true./), (/2,3/)), &
    6, 1)

  print*, "test.3.1.1"
  call integerdataarray("1 2 3 4 5", &
    (/1, 2, 3, 4, 5/), 5, 0)
  print*, "test.3.1.2"
  call integerdataarray("1 2 3 4", &
    (/1, 2, 3, 4, 0/), 4, -1)
  print*, "test.3.1.3"
  call integerdataarray("1 2 3 4 5 6", &
    (/1, 2, 3, 4, 5/), 5, 1)

  print*, "test.3.1.4"
  call integerdataarray("1, 2, 3, 4, 5", & 
    (/1, 2, 3, 4, 5/), 5, 0)
  print*, "test.3.1.5"
  call integerdataarray("1, 2, 3, 4,", & 
    (/1, 2, 3, 4, 0/), 4, -1)
  print*, "test.3.1.6"
  call integerdataarray("1, 2, 3, 4, 5, 6,", & 
    (/1, 2, 3, 4, 5/), 5, 1)

  print*, "test.4.1.1"
  call realspdataarray("1.0 -2.e5 +3.44e-3 4. 5.09090909", & 
    (/1.0e0, -2.e5, 3.44e-3, 4.e0, 5.09090909e0/), 5, 0)
  print*, "test.4.1.2"
  call realspdataarray("1.0 -2.e5 +3.44e-3 4.", & 
    (/1.0e0, -2.e5, 3.44e-3, 4.e0, 0.0e0/), 4, -1)
  print*, "test.4.1.3"
  call realspdataarray("1.0 -2.e5 +3.44e-3 4. 5.09090909 -0.3", & 
    (/1.0e0, -2.e5, 3.44e-3, 4.e0, 5.09090909e0/), 5, 1)

  print*, "test.4.1.4"
  call realspdataarray("1.0, -2.e5 +3.44e-3, 4. 5.09090909", & 
    (/1.0e0, -2.e5, 3.44e-3, 4.e0, 5.09090909e0/), 5, 0)
  print*, "test.4.1.5"
  call realspdataarray("1.0, -2.e5 +3.44e-3, 4.", &  
    (/1.0e0, -2.e5, 3.44e-3, 4.e0, 0.0e0/), 4, -1)
  print*, "test.4.1.6"
  call realspdataarray("1.0, -2.e5 +3.44e-3, 4. 5.09090909, -0.3", & 
    (/1.0e0, -2.e5, 3.44e-3, 4.e0, 5.09090909e0/), 5, 1)

  print*, "test.5.1.1"
  call realdpdataarray("1.0 -2.e5 +3.44e-3 4. 5.09090909", & 
    (/1.0d0, -2.d5, 3.44d-3, 4.d0, 5.09090909d0/), 5, 0)
  print*, "test.5.1.2"
  call realdpdataarray("1.0 -2.e5 +3.44e-3 4.", &  
    (/1.0d0, -2.d5, 3.44d-3, 4.d0, 0.0d0/), 4, -1)
  print*, "test.5.1.3"
  call realdpdataarray("1.0 -2.e5 +3.44e-3 4. 5.09090909 -0.3", & 
    (/1.0d0, -2.d5, 3.44d-3, 4.d0, 5.09090909d0/), 5, 1)

  print*, "test.5.1.4"
  call realdpdataarray("1.0, -2.e5 +3.44e-3, 4. 5.09090909", & 
    (/1.0d0, -2.d5, 3.44d-3, 4.d0, 5.09090909d0/), 5, 0)
  print*, "test.5.1.5"
  call realdpdataarray("1.0, -2.e5 +3.44e-3, 4.", &  
    (/1.0d0, -2.d5, 3.44d-3, 4.d0, 0.0d0/), 4, -1)
  print*, "test.5.1.6"
  call realdpdataarray("1.0, -2.e5 +3.44e-3, 4. 5.09090909, -0.3", & 
    (/1.0d0, -2.d5, 3.44d-3, 4.d0, 5.09090909d0/), 5, 1)

  print*, "test.6.1.1"
  call cmplxspdataarray("(1.0)+i(-2.e5) (+3.44e-3)+i(4.) (5.09090909)+i(1.0) (-2e5)+i(+3.44e-3) (4.)+i(5.09090909)", & 
    (/(1.0e0,-2.e5), (3.44e-3,4.e0), (5.09090909e0,1.0e0), (-2e5,3.44e-3), (4.,5.09090909)/), 5, 0)
  print*, "test.6.1.2"
  call cmplxspdataarray("(1.0)+i(-2.e5) (+3.44e-3)+i(4.) (5.09090909)+i(1.0) (-2e5)+i(+3.44e-3)", & 
    (/(1.0e0,-2.e5), (3.44e-3,4.e0), (5.09090909e0,1.0e0), (-2e5,3.44e-3), (0.,0.)/), 4, -1)
  print*, "test.6.1.3"
  call cmplxspdataarray("(1.0)+i(-2.e5) (+3.44e-3)+i(4.) (5.09090909)+i(1.0) (-2e5)+i(+3.44e-3) (4.)+i(5.09090909) (1)+i(1)", & 
    (/(1.0e0,-2.e5), (3.44e-3,4.e0), (5.09090909e0,1.0e0), (-2e5,3.44e-3), (4.,5.09090909)/), 5, 1)

  print*, "test.6.1.4"
  call cmplxspdataarray("1.0,-2.e5 +3.44e-3,4. 5.09090909,1.0 -2e5,+3.44e-3 4.,5.09090909", & 
    (/(1.0e0,-2.e5), (3.44e-3,4.e0), (5.09090909e0,1.0e0), (-2e5,3.44e-3), (4.,5.09090909)/), 5, 0)
  print*, "test.6.1.5"
  call cmplxspdataarray("1.0,-2.e5 +3.44e-3,4. 5.09090909,1.0 -2e5,+3.44e-3", & 
    (/(1.0e0,-2.e5), (3.44e-3,4.e0), (5.09090909e0,1.0e0), (-2e5,3.44e-3), (0.,0.)/), 4, -1)
  print*, "test.6.1.6"
  call cmplxspdataarray("1.0,-2.e5 +3.44e-3,4. 5.09090909,1.0 -2e5,+3.44e-3 4.,5.09090909 1,1", & 
    (/(1.0e0,-2.e5), (3.44e-3,4.e0), (5.09090909e0,1.0e0), (-2e5,3.44e-3), (4.,5.09090909)/), 5, 1)

contains


    subroutine stringdatascalar(string, array, num, iostat)
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: array
      integer, intent(in) :: num, iostat

      character(len=len(array)) :: temp
      integer :: n, i 

      call rts(string, temp, num=n, iostat=i)

      if (temp/=array) then
        print*, "Different array"
        stop
      endif
      if (i/=iostat) &
        print*, "Wrong iostat" 
      if (n/=num) &
        print*, "Wrong num"
      ! Don't test countrts on string scalar rts output,
      ! as in this case rts does not split on whitespace.
    end subroutine stringdatascalar

    subroutine logicaldatascalar(string, array, num, iostat) 
      character(len=*), intent(in) :: string
      logical, intent(in) :: array
      integer, intent(in) :: num, iostat

      logical :: temp 
      integer :: n, i 

      call rts(string, temp, num=n, iostat=i)

      if (temp.neqv.array) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat" 
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,.false.)/=num)).or. & 
          ((i>1).and.(countrts(string,.false.)/=0)))      &
        print*, "Countrts wrong", countrts(string,.false.)
    end subroutine logicaldatascalar

    subroutine integerdatascalar(string, array, num, iostat) 
      character(len=*), intent(in) :: string
      integer, intent(in) :: array
      integer, intent(in) :: num, iostat

      integer :: temp 
      integer :: n, i 

      call rts(string, temp, num=n, iostat=i)

      if (temp/=array) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat" 
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1)/=num)).or. & 
          ((i>1).and.(countrts(string,1)/=0)))      &
        print*, "Countrts wrong", countrts(string,1)
    end subroutine integerdatascalar

    subroutine realspdatascalar(string, array, num, iostat) 
      character(len=*), intent(in) :: string
      real(sp), intent(in) :: array
      integer, intent(in) :: num, iostat

      real(sp) :: temp 
      integer :: n, i 

      call rts(string, temp, num=n, iostat=i)

      if (array==0.0) then
        if (temp>1e-5) then
          print*, "Different array"
        endif   
      elseif (abs(temp-array)/array>1e-5) then
        print*, "Different array"
      endif  
      if (i/=iostat) &
        print*, "Wrong iostat" 
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1.0)/=num)).or. & 
          ((i>1).and.(countrts(string,1.0)/=0)))      &
        print*, "Countrts wrong", countrts(string,1.0)
    end subroutine realspdatascalar

    subroutine realdpdatascalar(string, array, num, iostat) 
      character(len=*), intent(in) :: string
      real(dp), intent(in) :: array
      integer, intent(in) :: num, iostat

      real(dp) :: temp 
      integer :: n, i 

      call rts(string, temp, num=n, iostat=i)

      if (array==0.0) then
        if (temp>1e-5) then
          print*, "Different array"
        endif   
      elseif (abs(temp-array)/array>1e-5) then
        print*, "Different array"
      endif  
      if (i/=iostat) &
        print*, "Wrong iostat" 
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1.0_dp)/=num)).or. & 
          ((i>1).and.(countrts(string,1.0_dp)/=0)))      &
        print*, "Countrts wrong", countrts(string,1.0_dp)
    end subroutine realdpdatascalar

    subroutine complexspdatascalar(string, array, num, iostat) 
      character(len=*), intent(in) :: string
      complex(sp), intent(in) :: array
      integer, intent(in) :: num, iostat

      complex(sp) :: temp 
      integer :: n, i 

      call rts(string, temp, num=n, iostat=i)

      if (real(array)==0.0) then
        if (real(temp)>1e-5) then
          print*, "Different array"
        endif   
      elseif (abs(real(temp)-real(array))/real(array)>1e-5) then
        print*, "Different array"
      endif
      if (aimag(array)==0.0) then
        if (aimag(temp)>1e-5) then
          print*, "Different array"
        endif   
      elseif (abs(aimag(temp)-aimag(array))/aimag(array)>1e-5) then
        print*, "Different array"
      endif
      if (i/=iostat) &
        print*, "Wrong iostat" 
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,cmplx(1.0,0.0))/=num)).or. & 
          ((i>1).and.(countrts(string,cmplx(1.0,0.0))/=0)))      &
        print*, "Countrts wrong", countrts(string,cmplx(1.0,0.0))
    end subroutine complexspdatascalar

    subroutine complexdpdatascalar(string, array, num, iostat) 
      character(len=*), intent(in) :: string
      complex(dp), intent(in) :: array
      integer, intent(in) :: num, iostat

      complex(dp) :: temp 
      integer :: n, i 

      call rts(string, temp, num=n, iostat=i)

      if (real(array)==0.0) then
        if (real(temp)>1e-5) then
          print*, "Different array"
        endif   
      elseif (abs(real(temp)-real(array))/real(array)>1e-5) then
        print*, "Different array"
      endif  
      if (aimag(array)==0.0) then
        if (aimag(temp)>1e-5) then
          print*, "Different array"
        endif   
      elseif (abs(aimag(temp)-aimag(array))/aimag(array)>1e-5) then
        print*, "Different array"
      endif  
      if (i/=iostat) &
        print*, "Wrong iostat" 
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,cmplx(1.0_dp,0.0_dp))/=num)).or. & 
          ((i>1).and.(countrts(string,cmplx(1.0_dp,0.0_dp))/=0)))      &
        print*, "Countrts wrong", countrts(string,cmplx(1.0_dp,0.0_dp))
    end subroutine complexdpdatascalar

    subroutine stringdataarray(string, array, num, iostat, sep, csv)
      character(len=*), intent(in) :: string
      character(len=*), dimension(:), intent(in) :: array
      integer, intent(in) :: num, iostat
      character, intent(in), optional :: sep
      logical, intent(in), optional :: csv

      character(len=len(array)) :: temp(size(array))
      integer :: n, i

      call rts(string, temp, separator=sep, csv=csv, num=n, iostat=i)

      if (any(temp/=array)) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string," ",sep,csv)/=num)).or. & 
          ((i>1).and.(countrts(string," ",sep,csv)/=0)))      &
        print*, "Countrts wrong", countrts(string," ",sep,csv)
    end subroutine stringdataarray

    subroutine stringdatamatrix(string, array, num, iostat, sep, csv)
      character(len=*), intent(in) :: string
      character(len=*), dimension(:, :), intent(in) :: array
      integer, intent(in) :: num, iostat
      character, intent(in), optional :: sep
      logical, intent(in), optional :: csv

      character(len=len(array)) :: temp(size(array,1),size(array,2))
      integer :: n, i

      call rts(string, temp, separator=sep, csv=csv, num=n, iostat=i)

      if (any(temp/=array)) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string," ",sep,csv)/=num)).or. & 
          ((i>1).and.(countrts(string," ",sep,csv)/=0)))      &
        print*, "Countrts wrong", countrts(string," ",sep,csv)
    end subroutine stringdatamatrix

    subroutine logicaldataarray(string, array, num, iostat)
      character(len=*), intent(in) :: string
      logical, dimension(:), intent(in) :: array
      integer, intent(in) :: num, iostat

      logical :: temp(size(array))
      integer :: n, i

      call rts(string, temp, n, i)

      if (any(temp.neqv.array)) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,.false.)/=num)).or. & 
          ((i>1).and.(countrts(string,.false.)/=0)))      &
        print*, "Countrts wrong", countrts(string,.false.)
    end subroutine logicaldataarray

    subroutine logicaldatamatrix(string, array, num, iostat)
      character(len=*), intent(in) :: string
      logical, dimension(:, :), intent(in) :: array
      integer, intent(in) :: num, iostat

      logical :: temp(size(array,1),size(array,2))
      integer :: n, i

      call rts(string, temp, n, i)

      if (any(temp.neqv.array)) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,.false.)/=num)).or. & 
          ((i>1).and.(countrts(string,.false.)/=0)))      &
        print*, "Countrts wrong", countrts(string,.false.)
    end subroutine logicaldatamatrix

    subroutine integerdataarray(string, array, num, iostat)
      character(len=*), intent(in) :: string
      integer, dimension(:), intent(in) :: array
      integer, intent(in) :: num, iostat

      integer :: temp(size(array))
      integer :: n, i

      call rts(string, temp, n, i)

      if (any(temp/=array)) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1)/=num)).or. & 
          ((i>1).and.(countrts(string,1)/=0)))      &
        print*, "Countrts wrong", countrts(string,1)
    end subroutine integerdataarray

    subroutine integerdatamatrix(string, array, num, iostat)
      character(len=*), intent(in) :: string
      integer, dimension(:, :), intent(in) :: array
      integer, intent(in) :: num, iostat

      integer :: temp(size(array,1),size(array,2))
      integer :: n, i

      call rts(string, temp, n, i)

      if (any(temp/=array)) &
        print*, "Different array"
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1.0)/=num)).or. & 
          ((i>1).and.(countrts(string,1.0)/=0)))      &
        print*, "Countrts wrong", countrts(string,1.0)
    end subroutine integerdatamatrix

    subroutine realspdataarray(string, array, num, iostat)
      character(len=*), intent(in) :: string
      real(sp), dimension(:), intent(in) :: array
      integer, intent(in) :: num, iostat

      real(sp) :: temp(size(array))
      integer :: n, i, j

      call rts(string, temp, n, i)

      do j = 1, size(array)
        if (array(j)==0.0) then
          if (temp(j)>1e-5) then
            print*, "Different array"
            exit
          endif
        elseif (abs(temp(j)-array(j))/array(j)>1e-5) then
          print*, "Different array"
          exit
        endif
      enddo 
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1.0)/=num)).or. & 
          ((i>1).and.(countrts(string,1.0)/=0)))      &
        print*, "Countrts wrong", countrts(string,1.0)
    end subroutine realspdataarray

    subroutine realspdatamatrix(string, array, num, iostat)
      character(len=*), intent(in) :: string
      real(sp), dimension(:, :), intent(in) :: array
      integer, intent(in) :: num, iostat

      real(sp) :: temp(size(array,1),size(array,2))
      integer :: n, i, j, k

      call rts(string, temp, n, i)

      loop: do j = 1, size(array,1)
        do k = 1, size(array,2)
          if (array(j,k)==0.0) then
            if (temp(j,k)>1e-5) then
              exit loop
            endif
           elseif (abs(temp(j,k)-array(j,k))/array(j,k)>1e-5) then
            print*, "Different array"
            exit loop
          endif
        enddo
      enddo loop
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1.0_dp)/=num)).or. & 
          ((i>1).and.(countrts(string,1.0_dp)/=0)))      &
        print*, "Countrts wrong", countrts(string,1.0_dp)
    end subroutine realspdatamatrix

    subroutine realdpdataarray(string, array, num, iostat)
      character(len=*), intent(in) :: string
      real(dp), dimension(:), intent(in) :: array
      integer, intent(in) :: num, iostat

      real(dp) :: temp(size(array))
      integer :: n, i, j

      call rts(string, temp, n, i)

      do j = 1, size(array)
        if (array(j)==0.0) then
          if (temp(j)>1e-5) then
            print*, "Different array"
            exit
          endif
        elseif (abs(temp(j)-array(j))/array(j)>1e-5) then
          print*, "Different array"
          exit
        endif
      enddo 
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1.0_dp)/=num)).or. & 
          ((i>1).and.(countrts(string,1.0_dp)/=0)))      &
        print*, "Countrts wrong", countrts(string,1.0_dp)
    end subroutine realdpdataarray

    subroutine realdpdatamatrix(string, array, num, iostat)
      character(len=*), intent(in) :: string
      real(dp), dimension(:, :), intent(in) :: array
      integer, intent(in) :: num, iostat

      real(dp) :: temp(size(array,1),size(array,2))
      integer :: n, i, j, k

      call rts(string, temp, n, i)

      loop: do j = 1, size(array,1)
        do k = 1, size(array,2)
          if (array(j,k)==0.0) then
            if (temp(j,k)>1e-5) then
              print*, "Different array"
              exit loop
            endif
           elseif (abs(temp(j,k)-array(j,k))/array(j,k)>1e-5) then
            print*, "Different array"
            exit loop
          endif
        enddo
      enddo loop
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,1.0_dp)/=num)).or. & 
          ((i>1).and.(countrts(string,1.0_dp)/=0)))      &
        print*, "Countrts wrong", countrts(string,1.0_dp)
    end subroutine realdpdatamatrix

    subroutine cmplxspdataarray(string, array, num, iostat)
      character(len=*), intent(in) :: string
      complex(sp), dimension(:), intent(in) :: array
      integer, intent(in) :: num, iostat

      complex(sp) :: temp(size(array))
      real(sp) :: a(2), t(2)
      integer :: n, i, j, m

      call rts(string, temp, n, i)

      do j = 1, size(array)
        a = transfer(array(j), a)
        t = transfer(temp(j), t)
        do m = 1, 2
          if (a(m)==0.0) then
            if (t(m)>1e-5) then
              print*, "Different array"
              exit
            endif
          elseif (abs(t(m)-a(m))/a(m)>1e-5) then
            print*, "Different array"
            exit
          endif
        enddo
      enddo 
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,cmplx(1.0,0.0))/=num)).or. & 
          ((i>1).and.(countrts(string,cmplx(1.0,0.0))/=0)))      &
        print*, "Countrts wrong", countrts(string,cmplx(1.0,0.0))
    end subroutine cmplxspdataarray

    subroutine cmplxspdatamatrix(string, array, num, iostat)
      character(len=*), intent(in) :: string
      complex(sp), dimension(:, :), intent(in) :: array
      integer, intent(in) :: num, iostat

      complex(sp) :: temp(size(array,1),size(array,2))
      real(sp) :: a(2), t(2)
      integer :: n, i, j, k, m

      call rts(string, temp, n, i)

      loop: do j = 1, size(array,1)
        do k = 1, size(array,2)
          a = transfer(array(j,k), a)
          t = transfer(temp(j,k), t)
          do m = 1, 2
            if (a(m)==0.0) then
              if (t(m)>1e-5) then
                print*, "Different array"
                exit loop
              endif
            elseif (abs(t(m)-a(m))/a(m)>1e-5) then
              print*, "Different array"
              exit loop
            endif
          enddo
        enddo
      enddo loop
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,cmplx(1.0,0.0))/=num)).or. & 
          ((i>1).and.(countrts(string,cmplx(1.0,0.0))/=0)))      &
        print*, "Countrts wrong", countrts(string,cmplx(1.0,0.0))
    end subroutine cmplxspdatamatrix

    subroutine cmplxdpdataarray(string, array, num, iostat)
      character(len=*), intent(in) :: string
      complex(dp), dimension(:), intent(in) :: array
      integer, intent(in) :: num, iostat

      complex(dp) :: temp(size(array))
      real(dp) :: a(2), t(2)
      integer :: n, i, j, m

      call rts(string, temp, n, i)

      do j = 1, size(array)
        a = transfer(array(j), a)
        t = transfer(temp(j), t)
        do m = 1, 2
          if (a(m)==0.0) then
            if (t(m)>1e-5) then
              print*, "Different array"
              exit
            endif
          elseif (abs(t(m)-a(m))/a(m)>1e-5) then
            print*, "Different array"
            exit
          endif
        enddo
      enddo 
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,cmplx(1.0_dp,0.0_dp))/=num)).or. & 
          ((i>1).and.(countrts(string,cmplx(1.0_dp,0.0_dp))/=0)))      &
        print*, "Countrts wrong", countrts(string,cmplx(1.0_dp,0.0_dp))
    end subroutine cmplxdpdataarray

    subroutine cmplxdpdatamatrix(string, array, num, iostat)
      character(len=*), intent(in) :: string
      complex(dp), dimension(:, :), intent(in) :: array
      integer, intent(in) :: num, iostat

      complex(dp) :: temp(size(array,1),size(array,2))
      real(dp) :: a(2), t(2)
      integer :: n, i, j, k, m

      call rts(string, temp, n, i)

      loop: do j = 1, size(array,1)
        do k = 1, size(array,2)
          a = transfer(array(j,k), a)
          t = transfer(temp(j,k), t)
          do m = 1, 2
            if (a(m)==0.0) then
              if (t(m)>1e-5) then
                print*, "Different array"
                exit loop
              endif
            elseif (abs(t(m)-a(m))/a(m)>1e-5) then
              print*, "Different array"
              exit loop
            endif
          enddo
        enddo
      enddo loop
      if (i/=iostat) &
        print*, "Wrong iostat"
      if (n/=num) &
        print*, "Wrong num"
      if (((i<=0).and.(countrts(string,cmplx(1.0_dp,0.0_dp))/=num)).or. & 
          ((i>1).and.(countrts(string,cmplx(1.0_dp,0.0_dp))/=0)))      &
        print*, "Countrts wrong", countrts(string,cmplx(1.0_dp,0.0_dp))
    end subroutine cmplxdpdatamatrix

end program short_test
