module avg_urr_xs_values

  use ace_header, only: Nuclide
  use global

  implicit none

contains

  subroutine set_avg_urr_xs(i_nuc)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    integer :: i_nuc ! nuclide index
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

    select case(nuc % MAT)
    case(9228)
      nuc % avg_urr_n =&
        (/&
        1.2138065E+01_8,&
        1.2117005E+01_8,&
        1.2055394E+01_8,&
        1.1997653E+01_8,&
        1.1944853E+01_8,&
        1.1899382E+01_8,&
        1.1855636E+01_8,&
        1.1815674E+01_8,&
        1.1778704E+01_8,&
        1.1651874E+01_8,&
        1.1590295E+01_8,&
        1.1448052E+01_8,&
        1.1348528E+01_8,&
        1.1324268E+01_8&
        /)

      nuc % avg_urr_f =&
        (/&
        5.5282348E+00_8,&
        5.2379606E+00_8,&
        4.4374975E+00_8,&
        3.9295148E+00_8,&
        3.5764610E+00_8,&
        3.3060439E+00_8,&
        3.0981449E+00_8,&
        2.9350292E+00_8,&
        2.7951273E+00_8,&
        2.4437594E+00_8,&
        2.3198605E+00_8,&
        2.0856048E+00_8,&
        1.9555880E+00_8,&
        1.9289955E+00_8&
        /)

      nuc % avg_urr_g =&
        (/&
        1.9943743E+00_8,&
        1.8883931E+00_8,&
        1.5884363E+00_8,&
        1.4014364E+00_8,&
        1.2686110E+00_8,&
        1.1706243E+00_8,&
        1.0934998E+00_8,&
        1.0324131E+00_8,&
        9.8174210E-01_8,&
        8.5472689E-01_8,&
        8.0842188E-01_8,&
        7.2382245E-01_8,&
        6.7598517E-01_8,&
        6.6768202E-01_8&
        /)

      nuc % avg_urr_x =&
        (/&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO&
        /)

    case(9237)
      nuc % avg_urr_n =&
        (/&
        1.385E+01_8,&
        1.369E+01_8,&
        1.357E+01_8,&
        1.342E+01_8,&
        1.328E+01_8,&
        1.315E+01_8,&
        1.306E+01_8,&
        1.305E+01_8,&
        1.293E+01_8,&
        1.279E+01_8,&
        1.266E+01_8,&
        1.244E+01_8,&
        1.224E+01_8,&
        1.207E+01_8,&
        1.191E+01_8,&
        1.163E+01_8,&
        1.139E+01_8,&
        1.129E+01_8&
        /)

      nuc % avg_urr_f =&
        (/&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO&
        /)

      nuc % avg_urr_g =&
        (/&
        5.296E-01_8,&
        4.957E-01_8,&
        4.681E-01_8,&
        4.373E-01_8,&
        4.065E-01_8,&
        3.826E-01_8,&
        3.631E-01_8,&
        3.623E-01_8,&
        3.194E-01_8,&
        2.888E-01_8,&
        2.641E-01_8,&
        2.286E-01_8,&
        2.043E-01_8,&
        1.871E-01_8,&
        1.739E-01_8,&
        1.566E-01_8,&
        1.459E-01_8,&
        1.423E-01_8&
        /)

      nuc % avg_urr_x =&
        (/&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        ZERO,&
        6.315E-02_8,&
        1.304E-01_8,&
        1.936E-01_8,&
        3.016E-01_8,&
        3.890E-01_8,&
        4.615E-01_8,&
        5.229E-01_8,&
        6.207E-01_8,&
        6.959E-01_8,&
        7.244E-01_8&
        /)

    case default
      if (nuc % LSSF == 1) then
        write(*, '(A40,I4,A66)') &
          & 'Averaged URR cross sections for MAT ', nuc % MAT, &
          & ' are not pre-computed and must be calculated at initialization'
      end if

    end select

  end subroutine set_avg_urr_xs

end module avg_urr_xs_values
