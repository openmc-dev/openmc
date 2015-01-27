module avg_urr_xs_values

  use global
  use unresolved, only: Isotope, isotopes

  implicit none

contains

  subroutine set_avg_urr_xs(i)

    type(Isotope), pointer :: tope => null() ! nuclide object pointer
    integer :: i ! isotope index
!$omp threadprivate(nuc)

    tope => isotopes(i)

    select case(tope % MAT)
    case(9228)
      tope % avg_urr_n =&
        (/&
        1.208415E+01_8, &
        1.206545E+01_8, &
        1.200110E+01_8, &
        1.194555E+01_8, &
        1.189861E+01_8, &
        1.184984E+01_8, &
        1.180893E+01_8, &
        1.176896E+01_8, &
        1.173066E+01_8, &
        1.160858E+01_8, &
        1.155669E+01_8, &
        1.142788E+01_8, &
        1.133661E+01_8, &
        1.131365E+01_8  &
        /)

      tope % avg_urr_f =&
        (/&
        5.636060E+00_8, &
        5.354305E+00_8, &
        4.569279E+00_8, &
        4.071116E+00_8, &
        3.728682E+00_8, &
        3.466185E+00_8, &
        3.267486E+00_8, &
        3.111314E+00_8, &
        2.975452E+00_8, &
        2.654352E+00_8, &
        2.529428E+00_8, &
        2.322383E+00_8, &
        2.204007E+00_8, &
        2.180881E+00_8  &
        /)

      tope % avg_urr_g =&
        (/&
        2.034389E+00_8, &
        1.926581E+00_8, &
        1.636855E+00_8, &
        1.453476E+00_8, &
        1.325482E+00_8, &
        1.231341E+00_8, &
        1.158980E+00_8, &
        1.101734E+00_8, &
        1.051473E+00_8, &
        9.328347E-01_8, &
        8.891014E-01_8, &
        8.109868E-01_8, &
        7.656194E-01_8, &
        7.579173E-01_8  &
        /)

      tope % avg_urr_x =&
        (/&
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO  &
        /)

    case(9237)
      tope % avg_urr_n =&
        (/&
        1.382182E+01_8, &
        1.366665E+01_8, &
        1.354169E+01_8, &
        1.339594E+01_8, &
        1.324607E+01_8, &
        1.311889E+01_8, &
        1.300442E+01_8, &
        1.299837E+01_8, &
        1.291292E+01_8, &
        1.274591E+01_8, &
        1.261738E+01_8, &
        1.239658E+01_8, &
        1.219202E+01_8, &
        1.201576E+01_8, &
        1.185615E+01_8, &
        1.158227E+01_8, &
        1.133988E+01_8, &
        1.123894E+01_8  &
        /)

      tope % avg_urr_f =&
        (/&
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO  &
        /)

      tope % avg_urr_g =&
        (/&
        5.310637E-01_8, &
        4.961408E-01_8, &
        4.686683E-01_8, &
        4.378293E-01_8, &
        4.069400E-01_8, &
        3.827348E-01_8, &
        3.637040E-01_8, &
        3.625148E-01_8, &
        3.456107E-01_8, &
        2.888114E-01_8, &
        2.641575E-01_8, &
        2.286368E-01_8, &
        2.044441E-01_8, &
        1.870973E-01_8, &
        1.741106E-01_8, &
        1.566557E-01_8, &
        1.460217E-01_8, &
        1.424347E-01_8  &
        /)

      tope % avg_urr_x =&
        (/&
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        ZERO, &
        2.279359E-03_8, &
        1.302187E-01_8, &
        1.931377E-01_8, &
        3.014849E-01_8, &
        3.890843E-01_8, &
        4.609592E-01_8, &
        5.229405E-01_8, &
        6.199341E-01_8, &
        6.972089E-01_8, &
        7.249919E-01_8  &
        /)
    case default
      if (tope % LSSF == 1) then
        write(*, '(A40,I4,A66)') &
          & 'Averaged URR cross sections for MAT ', tope % MAT, &
          & ' are not pre-computed and must be calculated at initialization'
      end if

    end select

  end subroutine set_avg_urr_xs

end module avg_urr_xs_values
