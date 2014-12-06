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
        0.000E+00_8,&
        0.000E+00_8,&
        0.000E+00_8,&
        0.000E+00_8,&
        0.000E+00_8,&
        0.000E+00_8,&
        0.000E+00_8,&
        0.000E+00_8,&
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
      write(*, '(A40,I4,A66)') 'Averaged URR cross sections for MAT ', nuc % MAT,&
        & ' are not pre-computed and must be calculated at initialization'

    end select

  end subroutine set_avg_urr_xs

end module avg_urr_xs_values
