module endf

  use global, only: int_to_str

contains

!=====================================================================
! REACTION_NAME gives the name of the reaction for a given MT value
!=====================================================================

  function reaction_name(MT) result(string)

    integer, intent(in) :: MT
    character(20)       :: string

    select case (MT)
    case (1)
       string = '(n,total)'
    case (2)
       string = '(elastic)'
    case (4)
       string = '(n,level)'
    case (11)
       string = '(n,2nd)'
    case (16)
       string = '(n,2n)'
    case (17)
       string = '(n,3n)'
    case (18)
       string = '(n,fission)'
    case (19)
       string = '(n,f)'
    case (20)
       string = '(n,2nf)'
    case (21)
       string = '(n,2nf)'
    case (22)
       string = '(n,na)'
    case (23)
       string = '(n,n3a)'
    case (24)
       string = '(n,2na)'
    case (25)
       string = '(n,3na)'
    case (28)
       string = '(n,np)'
    case (29)
       string = '(n,n2a)'
    case (30)
       string = '(n,2n2a)'
    case (32)
       string = '(n,nd)'
    case (33)
       string = '(n,nt)'
    case (34)
       string = '(n,nHe-3)'
    case (35)
       string = '(n,nd3a)'
    case (36)
       string = '(n,nt2a)'
    case (37)
       string = '(n,4n)'
    case (38)
       string = '(n,3nf)'
    case (41)
       string = '(n,2np)'
    case (42)
       string = '(n,3np)'
    case (44)
       string = '(n,2np)'
    case (45)
       string = '(n,npa)'
    case (51 : 90)
       string = '(n,n' // trim(int_to_str(MT-50)) // ')'
    case (91)
       string = '(n,nc)'
    case (102)
       string = '(n,gamma)'
    case (103)
       string = '(n,p)'
    case (104)
       string = '(n,d)'
    case (105)
       string = '(n,t)'
    case (106)
       string = '(n,3He)'
    case (107)
       string = '(n,a)'
    case (108)
       string = '(n,2a)'
    case (109)
       string = '(n,3a)'
    case (111)
       string = '(n,2p)'
    case (112)
       string = '(n,pa)'
    case (113)
       string = '(n,t2a)'
    case (114)
       string = '(n,d2a)'
    case (115)
       string = '(n,pd)'
    case (116)
       string = '(n,pt)'
    case (117)
       string = '(n,da)'
    case (201)
       string = '(n,Xn)'
    case (202)
       string = '(n,Xgamma)'
    case (203)
       string = '(n,Xp)'
    case (204)
       string = '(n,Xd)'
    case (205)
       string = '(n,Xt)'
    case (206)
       string = '(n,X3He)'
    case (207)
       string = '(n,Xa)'
    case (444)
       string = '(damage)'
    case (600 : 648)
       string = '(n,p' // trim(int_to_str(MT-600)) // ')'
    case (649)
       string = '(n,pc)'
    case (650 : 698)
       string = '(n,d' // trim(int_to_str(MT-650)) // ')'
    case (699)
       string = '(n,dc)'
    case (700 : 748)
       string = '(n,t' // trim(int_to_str(MT-600)) // ')'
    case (749)
       string = '(n,tc)'
    case (750 : 798)
       string = '(n,3He' // trim(int_to_str(MT-650)) // ')'
    case (799)
       string = '(n,3Hec)'
    case (800 : 848)
       string = '(n,a' // trim(int_to_str(MT-800)) // ')'
    case (849)
       string = '(n,tc)'
    case default
       string = 'MT=' // trim(int_to_str(MT))
    end select

  end function reaction_name
       
end module endf
