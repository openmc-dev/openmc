module avg_urr_xs_values

  use global
  use unresolved, only: Isotope, isotopes, represent_params

  implicit none

contains

  subroutine set_avg_urr_xs(i)

    type(Isotope), pointer :: tope => null() ! nuclide object pointer
    integer :: i ! isotope index
!$omp threadprivate(nuc)

    tope => isotopes(i)

    if (represent_params == CONTINUOUS) then
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

      case(9437)
        tope % avg_urr_n =&
          (/&
0          /)

        tope % avg_urr_f =&
          (/&
 0         /)

        tope % avg_urr_g =&
          (/&
  0        /)

        tope % avg_urr_x =&
          (/&
   0       /)

      case default
        if (tope % LSSF == 1) then
          write(*, '(A40,I4,A66)') &
            & 'Averaged URR cross sections for MAT ', tope % MAT, &
            & ' are not pre-computed and must be calculated at initialization'
        end if
      end select

    else if (represent_params == DISCRETE) then
      select case(tope % MAT)
      case(9228)
        tope % avg_urr_n =&
          (/&
          12.08406_8, &
          12.06536_8, &
          12.00102_8, &
          11.94548_8, &
          11.89855_8, &
          11.84985_8, &
          11.80883_8, &
          11.76891_8, &
          11.73063_8, &
          11.60849_8, &
          11.55670_8, &
          11.42774_8, &
          11.33675_8, &
          11.31362_8  &
          /)

        tope % avg_urr_f =&
          (/&
          5.636044_8, &
          5.354288_8, &
          4.569263_8, &
          4.071107_8, &
          3.728670_8, &
          3.466121_8, &
          3.267571_8, &
          3.111236_8, &
          2.975452_8, &
          2.654328_8, &
          2.529419_8, &
          2.322576_8, &
          2.204009_8, &
          2.180565_8  &
          /)

        tope % avg_urr_g =&
          (/&
          2.034381_8, &
          1.926577_8, &
          1.636851_8, &
          1.453473_8, &
          1.325479_8, &
          1.231350_8, &
          1.158947_8, &
          1.101774_8, &
          1.051463_8, &
          0.9328223_8, &
          0.8891086_8, &
          0.8109948_8, &
          0.7656376_8, &
          0.7578801_8  &
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
          13.82207_8, &
          13.66642_8, &
          13.54062_8, &
          13.39396_8, &
          13.24554_8, &
          13.11851_8, &
          13.00414_8, &
          12.99775_8, &
          12.87910_8, &
          12.74689_8, &
          12.61919_8, &
          12.38779_8, &
          12.18756_8, &
          12.01593_8, &
          11.85803_8, &
          11.58140_8, &
          11.33963_8, &
          11.24521_8  &
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
          0.5311968_8, &
          0.4960280_8, &
          0.4686543_8, &
          0.4377895_8, &
          0.4069119_8, &
          0.3827407_8, &
          0.3637049_8, &
          0.3625604_8, &
          0.3191216_8, &
          0.2890935_8, &
          0.2641849_8, &
          0.2285740_8, &
          0.2041515_8, &
          0.1870316_8, &
          0.1740800_8, &
          0.1565132_8, &
          0.1458820_8, &
          0.1424402_8  &
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
          0.06297613_8, &
          0.1303260_8, &
          0.1934372_8, &
          0.3015513_8, &
          0.3886899_8, &
          0.4616134_8, &
          0.5228022_8, &
          0.6208314_8, &
          0.6960096_8, &
          0.7254137_8  &
          /)

      case(9437)
        tope % avg_urr_n =&
          (/&
0          /)

        tope % avg_urr_f =&
          (/&
0          /)

        tope % avg_urr_g =&
          (/&
0          /)

        tope % avg_urr_x =&
          (/&
0          /)

      case default
        if (tope % LSSF == 1) then
          write(*, '(A40,I4,A66)') &
            & 'Averaged URR cross sections for MAT ', tope % MAT, &
            & ' are not pre-computed and must be calculated at initialization'
        end if
      end select
    end if

  end subroutine set_avg_urr_xs

end module avg_urr_xs_values
