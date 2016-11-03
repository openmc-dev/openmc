module cross_section

  use algorithm,        only: binary_search
  use constants
  use energy_grid,      only: log_spacing
  use error,            only: fatal_error
  use global
  use list_header,      only: ListElemInt
  use material_header,  only: Material
  use math,             only: faddeeva, broaden_wmp_polynomials
  use multipole_header, only: FORM_RM, FORM_MLBW, MP_EA, RM_RT, RM_RA, RM_RF, &
                              MLBW_RT, MLBW_RX, MLBW_RA, MLBW_RF, FIT_T, FIT_A,&
                              FIT_F, MultipoleArray
  use nuclide_header
  use particle_header,  only: Particle
  use random_lcg,       only: prn, future_prn, prn_set_stream
  use sab_header,       only: SAlphaBeta

  implicit none

contains

!===============================================================================
! CALCULATE_XS determines the macroscopic cross sections for the material the
! particle is currently traveling through.
!===============================================================================

  subroutine calculate_xs(p)

    type(Particle), intent(inout) :: p

    integer :: i             ! loop index over nuclides
    integer :: i_nuclide     ! index into nuclides array
    integer :: i_sab         ! index into sab_tables array
    integer :: j             ! index in mat % i_sab_nuclides
    integer :: i_grid        ! index into logarithmic mapping array or material
                             ! union grid
    real(8) :: atom_density  ! atom density of a nuclide
    logical :: check_sab     ! should we check for S(a,b) table?

    ! Set all material macroscopic cross sections to zero
    material_xs % total          = ZERO
    material_xs % elastic        = ZERO
    material_xs % absorption     = ZERO
    material_xs % fission        = ZERO
    material_xs % nu_fission     = ZERO

    ! Exit subroutine if material is void
    if (p % material == MATERIAL_VOID) return

    associate (mat => materials(p % material))
      ! Find energy index on energy grid
      i_grid = int(log(p % E/energy_min_neutron)/log_spacing)

      ! Determine if this material has S(a,b) tables
      check_sab = (mat % n_sab > 0)

      ! Initialize position in i_sab_nuclides
      j = 1

      ! Add contribution from each nuclide in material
      do i = 1, mat % n_nuclides
        ! ========================================================================
        ! CHECK FOR S(A,B) TABLE

        i_sab = 0

        ! Check if this nuclide matches one of the S(a,b) tables specified -- this
        ! relies on i_sab_nuclides being in sorted order
        if (check_sab) then
          if (i == mat % i_sab_nuclides(j)) then
            ! Get index in sab_tables
            i_sab = mat % i_sab_tables(j)

            ! If particle energy is greater than the highest energy for the S(a,b)
            ! table, don't use the S(a,b) table
            if (p % E > sab_tables(i_sab) % data(1) % threshold_inelastic) i_sab = 0

            ! Increment position in i_sab_nuclides
            j = j + 1

            ! Don't check for S(a,b) tables if there are no more left
            if (j > mat % n_sab) check_sab = .false.
          end if
        end if

        ! ========================================================================
        ! CALCULATE MICROSCOPIC CROSS SECTION

        ! Determine microscopic cross sections for this nuclide
        i_nuclide = mat % nuclide(i)

        ! Calculate microscopic cross section for this nuclide
        if (p % E /= micro_xs(i_nuclide) % last_E &
             .or. p % sqrtkT /= micro_xs(i_nuclide) % last_sqrtkT) then
          call calculate_nuclide_xs(i_nuclide, i_sab, p % E, i_grid, p % sqrtkT)
        else if (i_sab /= micro_xs(i_nuclide) % last_index_sab) then
          call calculate_nuclide_xs(i_nuclide, i_sab, p % E, i_grid, p % sqrtkT)
        end if

        ! ========================================================================
        ! ADD TO MACROSCOPIC CROSS SECTION

        ! Copy atom density of nuclide in material
        atom_density = mat % atom_density(i)

        ! Add contributions to material macroscopic total cross section
        material_xs % total = material_xs % total + &
             atom_density * micro_xs(i_nuclide) % total

        ! Add contributions to material macroscopic scattering cross section
        material_xs % elastic = material_xs % elastic + &
             atom_density * micro_xs(i_nuclide) % elastic

        ! Add contributions to material macroscopic absorption cross section
        material_xs % absorption = material_xs % absorption + &
             atom_density * micro_xs(i_nuclide) % absorption

        ! Add contributions to material macroscopic fission cross section
        material_xs % fission = material_xs % fission + &
             atom_density * micro_xs(i_nuclide) % fission

        ! Add contributions to material macroscopic nu-fission cross section
        material_xs % nu_fission = material_xs % nu_fission + &
             atom_density * micro_xs(i_nuclide) % nu_fission
      end do
    end associate

  end subroutine calculate_xs

!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array at the energy of the given particle
!===============================================================================

  subroutine calculate_nuclide_xs(i_nuclide, i_sab, E, i_log_union, sqrtkT)
    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy
    integer, intent(in) :: i_log_union ! index into logarithmic mapping array or
                                       ! material union energy grid
    real(8), intent(in) :: sqrtkT    ! Square root of kT, material dependent

    logical :: use_mp ! true if XS can be calculated with windowed multipole
    integer :: i_temp ! index for temperature
    integer :: i_grid ! index on nuclide energy grid
    integer :: i_low  ! lower logarithmic mapping index
    integer :: i_high ! upper logarithmic mapping index
    real(8) :: f      ! interp factor on nuclide energy grid
    real(8) :: kT     ! temperature in eV
    real(8) :: sigT, sigA, sigF ! Intermediate multipole variables

    associate (nuc => nuclides(i_nuclide))
      ! Check to see if there is multipole data present at this energy
      use_mp = .false.
      if (nuc % mp_present) then
        if (E >= nuc % multipole % start_E .and. &
             E <= nuc % multipole % end_E) then
          use_mp = .true.
        end if
      end if

      ! Evaluate multipole or interpolate
      if (use_mp) then
        ! Call multipole kernel
        call multipole_eval(nuc % multipole, E, sqrtkT, sigT, sigA, sigF)

        micro_xs(i_nuclide) % total = sigT
        micro_xs(i_nuclide) % absorption = sigA
        micro_xs(i_nuclide) % elastic = sigT - sigA

        if (nuc % fissionable) then
          micro_xs(i_nuclide) % fission = sigF
          micro_xs(i_nuclide) % nu_fission = sigF * nuc % nu(E, EMISSION_TOTAL)
        else
          micro_xs(i_nuclide) % fission    = ZERO
          micro_xs(i_nuclide) % nu_fission = ZERO
        end if

        ! Ensure these values are set
        ! Note, the only time either is used is in one of 4 places:
        ! 1. physics.F90 - scatter - For inelastic scatter.
        ! 2. physics.F90 - sample_fission - For partial fissions.
        ! 3. tally.F90 - score_general - For tallying on MTxxx reactions.
        ! 4. cross_section.F90 - calculate_urr_xs - For unresolved purposes.
        ! It is worth noting that none of these occur in the resolved
        ! resonance range, so the value here does not matter.  index_temp is
        ! set to -1 to force a segfault in case a developer messes up and tries
        ! to use it with multipole.
        micro_xs(i_nuclide) % index_temp    = -1
        micro_xs(i_nuclide) % index_grid    = 0
        micro_xs(i_nuclide) % interp_factor = ZERO

      else
        ! Find the appropriate temperature index.
        kT = sqrtkT**2
        select case (temperature_method)
        case (TEMPERATURE_NEAREST)
          i_temp = minloc(abs(nuclides(i_nuclide) % kTs - kT), dim=1)

        case (TEMPERATURE_INTERPOLATION)
          ! Find temperatures that bound the actual temperature
          do i_temp = 1, size(nuc % kTs) - 1
            if (nuc % kTs(i_temp) <= kT .and. kT < nuc % kTs(i_temp + 1)) exit
          end do

          ! Randomly sample between temperature i and i+1
          f = (kT - nuc % kTs(i_temp)) / &
               (nuc % kTs(i_temp + 1) - nuc % kTs(i_temp))
          if (f > prn()) i_temp = i_temp + 1
        end select

        associate (grid => nuc % grid(i_temp), xs => nuc % sum_xs(i_temp))
          ! Determine the energy grid index using a logarithmic mapping to
          ! reduce the energy range over which a binary search needs to be
          ! performed

          if (E < grid % energy(1)) then
            i_grid = 1
          elseif (E > grid % energy(size(grid % energy))) then
            i_grid = size(grid % energy) - 1
          else
            ! Determine bounding indices based on which equal log-spaced
            ! interval the energy is in
            i_low  = grid % grid_index(i_log_union)
            i_high = grid % grid_index(i_log_union + 1) + 1

            ! Perform binary search over reduced range
            i_grid = binary_search(grid % energy(i_low:i_high), &
                 i_high - i_low + 1, E) + i_low - 1
          end if

          ! check for rare case where two energy points are the same
          if (grid % energy(i_grid) == grid % energy(i_grid + 1)) &
               i_grid = i_grid + 1

          ! calculate interpolation factor
          f = (E - grid % energy(i_grid)) / &
               (grid % energy(i_grid + 1) - grid % energy(i_grid))

          micro_xs(i_nuclide) % index_temp    = i_temp
          micro_xs(i_nuclide) % index_grid    = i_grid
          micro_xs(i_nuclide) % interp_factor = f

          ! Initialize nuclide cross-sections to zero
          micro_xs(i_nuclide) % fission    = ZERO
          micro_xs(i_nuclide) % nu_fission = ZERO

          ! Calculate microscopic nuclide total cross section
          micro_xs(i_nuclide) % total = (ONE - f) * xs % total(i_grid) &
               + f * xs % total(i_grid + 1)

          ! Calculate microscopic nuclide elastic cross section
          micro_xs(i_nuclide) % elastic = (ONE - f) * xs % elastic(i_grid) &
               + f * xs % elastic(i_grid + 1)

          ! Calculate microscopic nuclide absorption cross section
          micro_xs(i_nuclide) % absorption = (ONE - f) * xs % absorption( &
               i_grid) + f * xs % absorption(i_grid + 1)

          if (nuc % fissionable) then
            ! Calculate microscopic nuclide total cross section
            micro_xs(i_nuclide) % fission = (ONE - f) * xs % fission(i_grid) &
                 + f * xs % fission(i_grid + 1)

            ! Calculate microscopic nuclide nu-fission cross section
            micro_xs(i_nuclide) % nu_fission = (ONE - f) * xs % nu_fission( &
                 i_grid) + f * xs % nu_fission(i_grid + 1)
          end if
        end associate
      end if

      ! Initialize sab treatment to false
      micro_xs(i_nuclide) % index_sab   = NONE
      micro_xs(i_nuclide) % elastic_sab = ZERO

      ! Initialize URR probability table treatment to false
      micro_xs(i_nuclide) % use_ptable  = .false.

      ! If there is S(a,b) data for this nuclide, we need to do a few
      ! things. Since the total cross section was based on non-S(a,b) data, we
      ! need to correct it by subtracting the non-S(a,b) elastic cross section and
      ! then add back in the calculated S(a,b) elastic+inelastic cross section.

      if (i_sab > 0) call calculate_sab_xs(i_nuclide, i_sab, E, sqrtkT)

      ! if the particle is in the unresolved resonance range and there are
      ! probability tables, we need to determine cross sections from the table

      if (urr_ptables_on .and. nuc % urr_present .and. .not. use_mp) then
        if (E > nuc % urr_data(i_temp) % energy(1) .and. E < nuc % &
             urr_data(i_temp) % energy(nuc % urr_data(i_temp) % n_energy)) then
          call calculate_urr_xs(i_nuclide, i_temp, E)
        end if
      end if

      micro_xs(i_nuclide) % last_E = E
      micro_xs(i_nuclide) % last_index_sab = i_sab
      micro_xs(i_nuclide) % last_sqrtkT = sqrtkT
    end associate

  end subroutine calculate_nuclide_xs

!===============================================================================
! CALCULATE_SAB_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range. These cross sections replace
! whatever data were taken from the normal Nuclide table.
!===============================================================================

  subroutine calculate_sab_xs(i_nuclide, i_sab, E, sqrtkT)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy
    real(8), intent(in) :: sqrtkT    ! temperature

    integer :: i_grid    ! index on S(a,b) energy grid
    integer :: i_temp    ! temperature index
    real(8) :: f         ! interp factor on S(a,b) energy grid
    real(8) :: inelastic ! S(a,b) inelastic cross section
    real(8) :: elastic   ! S(a,b) elastic cross section
    real(8) :: kT

    ! Set flag that S(a,b) treatment should be used for scattering
    micro_xs(i_nuclide) % index_sab = i_sab

    ! Determine temperature for S(a,b) table
    kT = sqrtkT**2
    if (temperature_method == TEMPERATURE_NEAREST) then
      ! If using nearest temperature, do linear search on temperature
      do i_temp = 1, size(sab_tables(i_sab) % kTs)
        if (abs(sab_tables(i_sab) % kTs(i_temp) - kT) < &
             K_BOLTZMANN*temperature_tolerance) exit
      end do
    else
      ! Find temperatures that bound the actual temperature
      do i_temp = 1, size(sab_tables(i_sab) % kTs) - 1
        if (sab_tables(i_sab) % kTs(i_temp) <= kT .and. &
             kT < sab_tables(i_sab) % kTs(i_temp + 1)) exit
      end do

      ! Randomly sample between temperature i and i+1
      f = (kT - sab_tables(i_sab) % kTs(i_temp)) / &
           (sab_tables(i_sab) % kTs(i_temp + 1) - sab_tables(i_sab) % kTs(i_temp))
      if (f > prn()) i_temp = i_temp + 1
    end if


    ! Get pointer to S(a,b) table
    associate (sab => sab_tables(i_sab) % data(i_temp))

      ! Get index and interpolation factor for inelastic grid
      if (E < sab % inelastic_e_in(1)) then
        i_grid = 1
        f = ZERO
      else
        i_grid = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, E)
        f = (E - sab%inelastic_e_in(i_grid)) / &
             (sab%inelastic_e_in(i_grid+1) - sab%inelastic_e_in(i_grid))
      end if

      ! Calculate S(a,b) inelastic scattering cross section
      inelastic = (ONE - f) * sab % inelastic_sigma(i_grid) + &
           f * sab % inelastic_sigma(i_grid + 1)

      ! Check for elastic data
      if (E < sab % threshold_elastic) then
        ! Determine whether elastic scattering is given in the coherent or
        ! incoherent approximation. For coherent, the cross section is
        ! represented as P/E whereas for incoherent, it is simply P

        if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
          if (E < sab % elastic_e_in(1)) then
            ! If energy is below that of the lowest Bragg peak, the elastic
            ! cross section will be zero
            elastic = ZERO
          else
            i_grid = binary_search(sab % elastic_e_in, &
                 sab % n_elastic_e_in, E)
            elastic = sab % elastic_P(i_grid) / E
          end if
        else
          ! Determine index on elastic energy grid
          if (E < sab % elastic_e_in(1)) then
            i_grid = 1
          else
            i_grid = binary_search(sab % elastic_e_in, &
                 sab % n_elastic_e_in, E)
          end if

          ! Get interpolation factor for elastic grid
          f = (E - sab%elastic_e_in(i_grid))/(sab%elastic_e_in(i_grid+1) - &
               sab%elastic_e_in(i_grid))

          ! Calculate S(a,b) elastic scattering cross section
          elastic = (ONE - f) * sab % elastic_P(i_grid) + &
               f * sab % elastic_P(i_grid + 1)
        end if
      else
        ! No elastic data
        elastic = ZERO
      end if
    end associate

    ! Correct total and elastic cross sections
    micro_xs(i_nuclide) % total = micro_xs(i_nuclide) % total - &
         micro_xs(i_nuclide) % elastic + inelastic + elastic
    micro_xs(i_nuclide) % elastic = inelastic + elastic

    ! Store S(a,b) elastic cross section for sampling later
    micro_xs(i_nuclide) % elastic_sab = elastic

    ! Save temperature index
    micro_xs(i_nuclide) % index_temp_sab = i_temp

  end subroutine calculate_sab_xs

!===============================================================================
! CALCULATE_URR_XS determines cross sections in the unresolved resonance range
! from probability tables
!===============================================================================

  subroutine calculate_urr_xs(i_nuclide, i_temp, E)
    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_temp    ! temperature index
    real(8), intent(in) :: E         ! energy

    integer :: i_energy     ! index for energy
    integer :: i_low        ! band index at lower bounding energy
    integer :: i_up         ! band index at upper bounding energy
    real(8) :: f            ! interpolation factor
    real(8) :: r            ! pseudo-random number
    real(8) :: elastic      ! elastic cross section
    real(8) :: capture      ! (n,gamma) cross section
    real(8) :: fission      ! fission cross section
    real(8) :: inelastic    ! inelastic cross section

    micro_xs(i_nuclide) % use_ptable = .true.

    associate (nuc => nuclides(i_nuclide), urr => nuclides(i_nuclide) % urr_data(i_temp))
      ! determine energy table
      i_energy = 1
      do
        if (E < urr % energy(i_energy + 1)) exit
        i_energy = i_energy + 1
      end do

      ! determine interpolation factor on table
      f = (E - urr % energy(i_energy)) / &
           (urr % energy(i_energy + 1) - urr % energy(i_energy))

      ! sample probability table using the cumulative distribution

      ! Random numbers for xs calculation are sampled from a separated stream.
      ! This guarantees the randomness and, at the same time, makes sure we reuse
      ! random number for the same nuclide at different temperatures, therefore
      ! preserving correlation of temperature in probability tables.
      call prn_set_stream(STREAM_URR_PTABLE)
      r = future_prn(int(i_nuclide, 8))
      call prn_set_stream(STREAM_TRACKING)

      i_low = 1
      do
        if (urr % prob(i_energy, URR_CUM_PROB, i_low) > r) exit
        i_low = i_low + 1
      end do
      i_up = 1
      do
        if (urr % prob(i_energy + 1, URR_CUM_PROB, i_up) > r) exit
        i_up = i_up + 1
      end do

      ! determine elastic, fission, and capture cross sections from probability
      ! table
      if (urr % interp == LINEAR_LINEAR) then
        elastic = (ONE - f) * urr % prob(i_energy, URR_ELASTIC, i_low) + &
             f * urr % prob(i_energy + 1, URR_ELASTIC, i_up)
        fission = (ONE - f) * urr % prob(i_energy, URR_FISSION, i_low) + &
             f * urr % prob(i_energy + 1, URR_FISSION, i_up)
        capture = (ONE - f) * urr % prob(i_energy, URR_N_GAMMA, i_low) + &
             f * urr % prob(i_energy + 1, URR_N_GAMMA, i_up)
      elseif (urr % interp == LOG_LOG) then
        ! Get logarithmic interpolation factor
        f = log(E / urr % energy(i_energy)) / &
             log(urr % energy(i_energy + 1) / urr % energy(i_energy))

        ! Calculate elastic cross section/factor
        elastic = ZERO
        if (urr % prob(i_energy, URR_ELASTIC, i_low) > ZERO .and. &
             urr % prob(i_energy + 1, URR_ELASTIC, i_up) > ZERO) then
          elastic = exp((ONE - f) * log(urr % prob(i_energy, URR_ELASTIC, &
               i_low)) + f * log(urr % prob(i_energy + 1, URR_ELASTIC, &
               i_up)))
        end if

        ! Calculate fission cross section/factor
        fission = ZERO
        if (urr % prob(i_energy, URR_FISSION, i_low) > ZERO .and. &
             urr % prob(i_energy + 1, URR_FISSION, i_up) > ZERO) then
          fission = exp((ONE - f) * log(urr % prob(i_energy, URR_FISSION, &
               i_low)) + f * log(urr % prob(i_energy + 1, URR_FISSION, &
               i_up)))
        end if

        ! Calculate capture cross section/factor
        capture = ZERO
        if (urr % prob(i_energy, URR_N_GAMMA, i_low) > ZERO .and. &
             urr % prob(i_energy + 1, URR_N_GAMMA, i_up) > ZERO) then
          capture = exp((ONE - f) * log(urr % prob(i_energy, URR_N_GAMMA, &
               i_low)) + f * log(urr % prob(i_energy + 1, URR_N_GAMMA, &
               i_up)))
        end if
      end if

      ! Determine treatment of inelastic scattering
      inelastic = ZERO
      if (urr % inelastic_flag > 0) then
        ! Get index on energy grid and interpolation factor
        i_energy = micro_xs(i_nuclide) % index_grid
        f = micro_xs(i_nuclide) % interp_factor

        ! Determine inelastic scattering cross section
        associate (xs => nuc % reactions(nuc % urr_inelastic) % xs(i_temp))
          if (i_energy >= xs % threshold) then
            inelastic = (ONE - f) * xs % value(i_energy - xs % threshold + 1) + &
                 f * xs % value(i_energy - xs % threshold + 2)
          end if
        end associate
      end if

      ! Multiply by smooth cross-section if needed
      if (urr % multiply_smooth) then
        elastic = elastic * micro_xs(i_nuclide) % elastic
        capture = capture * (micro_xs(i_nuclide) % absorption - &
             micro_xs(i_nuclide) % fission)
        fission = fission * micro_xs(i_nuclide) % fission
      end if

      ! Check for negative values
      if (elastic < ZERO) elastic = ZERO
      if (fission < ZERO) fission = ZERO
      if (capture < ZERO) capture = ZERO

      ! Set elastic, absorption, fission, and total cross sections. Note that the
      ! total cross section is calculated as sum of partials rather than using the
      ! table-provided value
      micro_xs(i_nuclide) % elastic = elastic
      micro_xs(i_nuclide) % absorption = capture + fission
      micro_xs(i_nuclide) % fission = fission
      micro_xs(i_nuclide) % total = elastic + inelastic + capture + fission

      ! Determine nu-fission cross section
      if (nuc % fissionable) then
        micro_xs(i_nuclide) % nu_fission = nuc % nu(E, EMISSION_TOTAL) * &
             micro_xs(i_nuclide) % fission
      end if
    end associate

  end subroutine calculate_urr_xs

!===============================================================================
! FIND_ENERGY_INDEX determines the index on the union energy grid at a certain
! energy
!===============================================================================

  pure function find_energy_index(mat, E) result(i)
    type(Material), intent(in) :: mat ! pointer to current material
    real(8),        intent(in) :: E   ! energy of particle
    integer                    :: i   ! energy grid index

    ! if the energy is outside of energy grid range, set to first or last
    ! index. Otherwise, do a binary search through the union energy grid.
    if (E <= mat % e_grid(1)) then
      i = 1
    elseif (E > mat % e_grid(mat % n_grid)) then
      i = mat % n_grid - 1
    else
      i = binary_search(mat % e_grid, mat % n_grid, E)
    end if

  end function find_energy_index

!===============================================================================
! MULTIPOLE_EVAL evaluates the windowed multipole equations for cross
! sections in the resolved resonance regions
!===============================================================================

  subroutine multipole_eval(multipole, E, sqrtkT, sigT, sigA, sigF)
    type(MultipoleArray), intent(in) :: multipole ! The windowed multipole
                                                  !  object to process.
    real(8), intent(in)              :: E         ! The energy at which to
                                                  !  evaluate the cross section
    real(8), intent(in)              :: sqrtkT    ! The temperature in the form
                                                  !  sqrt(kT), at which
                                                  !  to evaluate the XS.
    real(8), intent(out)             :: sigT      ! Total cross section
    real(8), intent(out)             :: sigA      ! Absorption cross section
    real(8), intent(out)             :: sigF      ! Fission cross section
    complex(8) :: psi_chi  ! The value of the psi-chi function for the
                           !  asymptotic form
    complex(8) :: c_temp   ! complex temporary variable
    complex(8) :: w_val    ! The faddeeva function evaluated at Z
    complex(8) :: Z        ! sqrt(atomic weight ratio / kT) * (sqrt(E) - pole)
    complex(8) :: sigT_factor(multipole % num_l)
    real(8) :: broadened_polynomials(multipole % fit_order + 1)
    real(8) :: sqrtE       ! sqrt(E), eV
    real(8) :: invE        ! 1/E, eV
    real(8) :: dopp        ! sqrt(atomic weight ratio / kT) = 1 / (2 sqrt(xi))
    real(8) :: temp        ! real temporary value
    integer :: i_pole      ! index of pole
    integer :: i_poly      ! index of curvefit
    integer :: i_window    ! index of window
    integer :: startw      ! window start pointer (for poles)
    integer :: endw        ! window end pointer

    ! ==========================================================================
    ! Bookkeeping

    ! Define some frequently used variables.
    sqrtE = sqrt(E)
    invE = ONE / E
    dopp = multipole % sqrtAWR / sqrtkT

    ! Locate us.
    i_window = floor((sqrtE - sqrt(multipole % start_E)) / multipole % spacing &
         + ONE)
    startw = multipole % w_start(i_window)
    endw = multipole % w_end(i_window)

    ! Fill in factors.
    if (startw <= endw) then
      call compute_sigT_factor(multipole, sqrtE, sigT_factor)
    end if

    ! Initialize the ouptut cross sections.
    sigT = ZERO
    sigA = ZERO
    sigF = ZERO

    ! ==========================================================================
    ! Add the contribution from the curvefit polynomial.

    if (sqrtkT /= ZERO .and. multipole % broaden_poly(i_window) == 1) then
      ! Broaden the curvefit.
      call broaden_wmp_polynomials(E, dopp, multipole % fit_order + 1, &
           broadened_polynomials)
      do i_poly = 1, multipole % fit_order+1
        sigT = sigT + multipole % curvefit(FIT_T, i_poly, i_window) &
             * broadened_polynomials(i_poly)
        sigA = sigA + multipole % curvefit(FIT_A, i_poly, i_window) &
             * broadened_polynomials(i_poly)
        sigF = sigF + multipole % curvefit(FIT_F, i_poly, i_window) &
             * broadened_polynomials(i_poly)
      end do
    else ! Evaluate as if it were a polynomial
      temp = invE
      do i_poly = 1, multipole % fit_order+1
        sigT = sigT + multipole % curvefit(FIT_T, i_poly, i_window) * temp
        sigA = sigA + multipole % curvefit(FIT_A, i_poly, i_window) * temp
        sigF = sigF + multipole % curvefit(FIT_F, i_poly, i_window) * temp
        temp = temp * sqrtE
      end do
    end if

    ! ==========================================================================
    ! Add the contribution from the poles in this window.

    if (sqrtkT == ZERO) then
      ! If at 0K, use asymptotic form.
      do i_pole = startw, endw
        psi_chi = -ONEI / (multipole % data(MP_EA, i_pole) - sqrtE)
        c_temp = psi_chi / E
        if (multipole % formalism == FORM_MLBW) then
          sigT = sigT + real(multipole % data(MLBW_RT, i_pole) * c_temp * &
                             sigT_factor(multipole % l_value(i_pole))) &
                      + real(multipole % data(MLBW_RX, i_pole) * c_temp)
          sigA = sigA + real(multipole % data(MLBW_RA, i_pole) * c_temp)
          sigF = sigF + real(multipole % data(MLBW_RF, i_pole) * c_temp)
        else if (multipole % formalism == FORM_RM) then
          sigT = sigT + real(multipole % data(RM_RT, i_pole) * c_temp * &
                             sigT_factor(multipole % l_value(i_pole)))
          sigA = sigA + real(multipole % data(RM_RA, i_pole) * c_temp)
          sigF = sigF + real(multipole % data(RM_RF, i_pole) * c_temp)
        end if
      end do
    else
      ! At temperature, use Faddeeva function-based form.
      if (endw >= startw) then
        do i_pole = startw, endw
          Z = (sqrtE - multipole % data(MP_EA, i_pole)) * dopp
          w_val = faddeeva(Z) * dopp * invE * SQRT_PI
          if (multipole % formalism == FORM_MLBW) then
            sigT = sigT + real((multipole % data(MLBW_RT, i_pole) * &
                          sigT_factor(multipole % l_value(i_pole)) + &
                          multipole % data(MLBW_RX, i_pole)) * w_val)
            sigA = sigA + real(multipole % data(MLBW_RA, i_pole) * w_val)
            sigF = sigF + real(multipole % data(MLBW_RF, i_pole) * w_val)
          else if (multipole % formalism == FORM_RM) then
            sigT = sigT + real(multipole % data(RM_RT, i_pole) * w_val * &
                               sigT_factor(multipole % l_value(i_pole)))
            sigA = sigA + real(multipole % data(RM_RA, i_pole) * w_val)
            sigF = sigF + real(multipole % data(RM_RF, i_pole) * w_val)
          end if
        end do
      end if
    end if
  end subroutine multipole_eval

!===============================================================================
! COMPUTE_SIGT_FACTOR calculates the sigT_factor, a factor inside of the sigT
! equation not present in the sigA and sigF equations.
!===============================================================================

  subroutine compute_sigT_factor(multipole, sqrtE, sigT_factor)
    type(MultipoleArray), intent(in)  :: multipole
    real(8),              intent(in)  :: sqrtE
    complex(8),           intent(out) :: sigT_factor(multipole % num_l)

    integer :: iL
    real(8) :: twophi(multipole % num_l)
    real(8) :: arg

    do iL = 1, multipole % num_l
      twophi(iL) = multipole % pseudo_k0RS(iL) * sqrtE
      if (iL == 2) then
        twophi(iL) = twophi(iL) - atan(twophi(iL))
      else if (iL == 3) then
        arg = 3.0_8 * twophi(iL) / (3.0_8 - twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      else if (iL == 4) then
        arg = twophi(iL) * (15.0_8 - twophi(iL)**2) &
             / (15.0_8 - 6.0_8 * twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      end if
    end do

    twophi = 2.0_8 * twophi
    sigT_factor = cmplx(cos(twophi), -sin(twophi), KIND=8)
  end subroutine compute_sigT_factor

!===============================================================================
! 0K_ELASTIC_XS determines the microscopic 0K elastic cross section
! for a given nuclide at the trial relative energy used in resonance scattering
!===============================================================================

  pure function elastic_xs_0K(E, nuc) result(xs_out)
    real(8),       intent(in) :: E      ! trial energy
    type(Nuclide), intent(in) :: nuc    ! target nuclide at temperature
    real(8)                   :: xs_out ! 0K xs at trial energy

    integer :: i_grid ! index on nuclide energy grid
    real(8) :: f      ! interp factor on nuclide energy grid

    ! Determine index on nuclide energy grid
    if (E < nuc % energy_0K(1)) then
      i_grid = 1
    elseif (E > nuc % energy_0K(nuc % n_grid_0K)) then
      i_grid = nuc % n_grid_0K - 1
    else
      i_grid = binary_search(nuc % energy_0K, nuc % n_grid_0K, E)
    end if

    ! check for rare case where two energy points are the same
    if (nuc % energy_0K(i_grid) == nuc % energy_0K(i_grid+1)) then
      i_grid = i_grid + 1
    end if

    ! calculate interpolation factor
    f = (E - nuc % energy_0K(i_grid)) &
         & / (nuc % energy_0K(i_grid + 1) - nuc % energy_0K(i_grid))

    ! Calculate microscopic nuclide elastic cross section
    xs_out = (ONE - f) * nuc % elastic_0K(i_grid) &
         & + f * nuc % elastic_0K(i_grid + 1)

  end function elastic_xs_0K

end module cross_section
