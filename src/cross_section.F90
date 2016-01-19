module cross_section

  use ace_header,       only: Nuclide, SAlphaBeta, Reaction, UrrData
  use constants
  use energy_grid,      only: grid_method, log_spacing
  use error,            only: fatal_error
  use fission,          only: nu_total
  use global
  use list_header,      only: ListElemInt
  use material_header,  only: Material
  use math,             only: w, broaden_n_polynomials
  use multipole_header, only: FORM_RM, FORM_MLBW, MP_EA, RM_RT, RM_RA, RM_RF, &
                              MLBW_RT, MLBW_RX, MLBW_RA, MLBW_RF, FIT_T, FIT_A, FIT_F, &
                              MultipoleArray, max_poly, max_L, max_poles
  use particle_header,  only: Particle
  use random_lcg,       only: prn
  use search,           only: binary_search

  implicit none

  ! Allocatable arrays for multipole that are allocated once for speed purposes
  complex(8), allocatable :: sigT_factor(:)
  real(8), allocatable :: twophi(:)
  real(8), allocatable :: broadened_polynomials(:)
  logical :: mp_already_alloc = .false.

!$omp threadprivate(sigT_factor, twophi, broadened_polynomials, mp_already_alloc)

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
    type(Material), pointer :: mat ! current material

    ! Set all material macroscopic cross sections to zero
    material_xs % total          = ZERO
    material_xs % elastic        = ZERO
    material_xs % absorption     = ZERO
    material_xs % fission        = ZERO
    material_xs % nu_fission     = ZERO

    ! Exit subroutine if material is void
    if (p % material == MATERIAL_VOID) return

    mat => materials(p % material)

    ! Find energy index on energy grid
    if (grid_method == GRID_MAT_UNION) then
      i_grid = find_energy_index(mat, p % E)
    else if (grid_method == GRID_LOGARITHM) then
      i_grid = int(log(p % E/energy_min_neutron)/log_spacing)
    end if

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
          if (p % E > sab_tables(i_sab) % threshold_inelastic) i_sab = 0

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
      if (p % E /= micro_xs(i_nuclide) % last_E .or. mat % sqrtkT /= micro_xs(i_nuclide) % last_sqrtkT) then
        call calculate_nuclide_xs(i_nuclide, i_sab, p % E, p % material, i, i_grid, mat % sqrtkT)
      else if (i_sab /= micro_xs(i_nuclide) % last_index_sab) then
        call calculate_nuclide_xs(i_nuclide, i_sab, p % E, p % material, i, i_grid, mat % sqrtkT)
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

  end subroutine calculate_xs

!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array at the energy of the given particle
!===============================================================================

  subroutine calculate_nuclide_xs(i_nuclide, i_sab, E, i_mat, i_nuc_mat, i_log_union, sqrtkT)
    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy
    integer, intent(in) :: i_mat     ! index into materials array
    integer, intent(in) :: i_nuc_mat ! index into nuclides array for a material
    integer, intent(in) :: i_log_union ! index into logarithmic mapping array or
                                       ! material union energy grid
    real(8), intent(in) :: sqrtkT    ! Square root of kT, material dependent

    integer :: i_grid ! index on nuclide energy grid
    integer :: i_low  ! lower logarithmic mapping index
    integer :: i_high ! upper logarithmic mapping index
    real(8) :: f      ! interp factor on nuclide energy grid
    real(8) :: sigT, sigA, sigF ! Intermediate multipole variables
    type(Nuclide),  pointer :: nuc
    type(Material), pointer :: mat

    ! Set pointer to nuclide and material
    nuc => nuclides(i_nuclide)
    mat => materials(i_mat)

    ! If MP, don't interpolate, it's all already baked in.
    if (nuc % mp_present .and. &
         (E >= nuc % multipole % start_E/1.0e6_8 .and.&
          E <= nuc % multipole % end_E/1.0e6_8)) then

      ! Call multipole kernel
      call multipole_eval(nuc % multipole, E, sqrtkT, sigT, sigA, sigF)

      micro_xs(i_nuclide) % total = sigT
      micro_xs(i_nuclide) % absorption = sigA
      micro_xs(i_nuclide) % elastic = sigT - sigA

      if (nuc % fissionable) then
        micro_xs(i_nuclide) % fission = sigF
        micro_xs(i_nuclide) % nu_fission = sigF * nu_total(nuc, E)
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
      ! resonance range, so the value here does not matter.
      micro_xs(i_nuclide) % index_grid    = 0
      micro_xs(i_nuclide) % interp_factor = ZERO
    else
      ! Determine index on nuclide energy grid
      select case (grid_method)
      case (GRID_MAT_UNION)

        i_grid = mat % nuclide_grid_index(i_nuc_mat, i_log_union)

      case (GRID_LOGARITHM)
        ! Determine the energy grid index using a logarithmic mapping to reduce
        ! the energy range over which a binary search needs to be performed

        if (E < nuc % energy(1)) then
          i_grid = 1
        elseif (E > nuc % energy(nuc % n_grid)) then
          i_grid = nuc % n_grid - 1
        else
          ! Determine bounding indices based on which equal log-spaced interval
          ! the energy is in
          i_low  = nuc % grid_index(i_log_union)
          i_high = nuc % grid_index(i_log_union + 1) + 1

          ! Perform binary search over reduced range
          i_grid = binary_search(nuc % energy(i_low:i_high), &
               i_high - i_low + 1, E) + i_low - 1
        end if

      case (GRID_NUCLIDE)
        ! Perform binary search on the nuclide energy grid in order to determine
        ! which points to interpolate between

        if (E <= nuc % energy(1)) then
          i_grid = 1
        elseif (E > nuc % energy(nuc % n_grid)) then
          i_grid = nuc % n_grid - 1
        else
          i_grid = binary_search(nuc % energy, nuc % n_grid, E)
        end if

      end select

      ! check for rare case where two energy points are the same
      if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) i_grid = i_grid + 1

      ! calculate interpolation factor
      f = (E - nuc%energy(i_grid))/(nuc%energy(i_grid+1) - nuc%energy(i_grid))

      micro_xs(i_nuclide) % index_grid    = i_grid
      micro_xs(i_nuclide) % interp_factor = f

      ! Initialize nuclide cross-sections to zero
      micro_xs(i_nuclide) % fission    = ZERO
      micro_xs(i_nuclide) % nu_fission = ZERO

      ! Calculate microscopic nuclide total cross section
      micro_xs(i_nuclide) % total = (ONE - f) * nuc % total(i_grid) &
           + f * nuc % total(i_grid+1)

      ! Calculate microscopic nuclide elastic cross section
      micro_xs(i_nuclide) % elastic = (ONE - f) * nuc % elastic(i_grid) &
           + f * nuc % elastic(i_grid+1)

      ! Calculate microscopic nuclide absorption cross section
      micro_xs(i_nuclide) % absorption = (ONE - f) * nuc % absorption( &
           i_grid) + f * nuc % absorption(i_grid+1)

      if (nuc % fissionable) then
        ! Calculate microscopic nuclide total cross section
        micro_xs(i_nuclide) % fission = (ONE - f) * nuc % fission(i_grid) &
             + f * nuc % fission(i_grid+1)

        ! Calculate microscopic nuclide nu-fission cross section
        micro_xs(i_nuclide) % nu_fission = (ONE - f) * nuc % nu_fission( &
             i_grid) + f * nuc % nu_fission(i_grid+1)
      end if
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

    if (i_sab > 0) call calculate_sab_xs(i_nuclide, i_sab, E)

    ! if the particle is in the unresolved resonance range and there are
    ! probability tables, we need to determine cross sections from the table

    if (urr_ptables_on .and. nuc % urr_present) then
      if (E > nuc % urr_data % energy(1) .and. &
           E < nuc % urr_data % energy(nuc % urr_data % n_energy)) then
        call calculate_urr_xs(i_nuclide, E)
      end if
    end if

    micro_xs(i_nuclide) % last_E = E
    micro_xs(i_nuclide) % last_index_sab = i_sab

  end subroutine calculate_nuclide_xs

!===============================================================================
! CALCULATE_SAB_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range. These cross sections replace
! whatever data were taken from the normal Nuclide table.
!===============================================================================

  subroutine calculate_sab_xs(i_nuclide, i_sab, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy

    integer :: i_grid    ! index on S(a,b) energy grid
    real(8) :: f         ! interp factor on S(a,b) energy grid
    real(8) :: inelastic ! S(a,b) inelastic cross section
    real(8) :: elastic   ! S(a,b) elastic cross section
    type(SAlphaBeta), pointer :: sab

    ! Set flag that S(a,b) treatment should be used for scattering
    micro_xs(i_nuclide) % index_sab = i_sab

    ! Get pointer to S(a,b) table
    sab => sab_tables(i_sab)

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

    ! Correct total and elastic cross sections
    micro_xs(i_nuclide) % total = micro_xs(i_nuclide) % total - &
         micro_xs(i_nuclide) % elastic + inelastic + elastic
    micro_xs(i_nuclide) % elastic = inelastic + elastic

    ! Store S(a,b) elastic cross section for sampling later
    micro_xs(i_nuclide) % elastic_sab = elastic

  end subroutine calculate_sab_xs

!===============================================================================
! CALCULATE_URR_XS determines cross sections in the unresolved resonance range
! from probability tables
!===============================================================================

  subroutine calculate_urr_xs(i_nuclide, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    real(8), intent(in) :: E         ! energy

    integer :: i            ! loop index
    integer :: i_energy     ! index for energy
    integer :: i_low        ! band index at lower bounding energy
    integer :: i_up         ! band index at upper bounding energy
    integer :: same_nuc_idx ! index of same nuclide
    real(8) :: f            ! interpolation factor
    real(8) :: r            ! pseudo-random number
    real(8) :: elastic      ! elastic cross section
    real(8) :: capture      ! (n,gamma) cross section
    real(8) :: fission      ! fission cross section
    real(8) :: inelastic    ! inelastic cross section
    logical :: same_nuc     ! do we know the xs for this nuclide at this energy?
    type(UrrData),  pointer :: urr
    type(Nuclide),  pointer :: nuc

    micro_xs(i_nuclide) % use_ptable = .true.

    ! get pointer to probability table
    nuc => nuclides(i_nuclide)
    urr => nuc % urr_data

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

    ! if we're dealing with a nuclide that we've previously encountered at
    ! this energy but a different temperature, use the original random number to
    ! preserve correlation of temperature in probability tables
    same_nuc = .false.
    do i = 1, nuc % nuc_list % size()
      if (E /= ZERO .and. E == micro_xs(nuc % nuc_list % data(i)) % last_E) then
        same_nuc = .true.
        same_nuc_idx = i
        exit
      end if
    end do

    if (same_nuc) then
      r = micro_xs(nuc % nuc_list % data(same_nuc_idx)) % last_prn
    else
      r = prn()
      micro_xs(i_nuclide) % last_prn = r
    end if

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
      associate (rxn => nuc % reactions(nuc % urr_inelastic))
        if (i_energy >= rxn % threshold) then
          inelastic = (ONE - f) * rxn % sigma(i_energy - rxn%threshold + 1) + &
               f * rxn % sigma(i_energy - rxn%threshold + 2)
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
      micro_xs(i_nuclide) % nu_fission = nu_total(nuc, E) * &
           micro_xs(i_nuclide) % fission
    end if

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
! MULTIPOLE_EVAL_ALLOCATE allocates fixed-length arrays that vary based on
! what nuclides are loaded into the problem
!===============================================================================

  subroutine multipole_eval_allocate()
    allocate(sigT_factor(max_L))
    allocate(twophi(max_L))
    allocate(broadened_polynomials(max_poly))

    mp_already_alloc = .true.
  end subroutine

!===============================================================================
! MULTIPOLE_EVAL evaluates the windowed multipole equations for cross
! sections in the resolved resonance regions
!===============================================================================

  subroutine multipole_eval(multipole, Emev, sqrtkT, sigT, sigA, sigF)
    type(MultipoleArray), intent(in) :: multipole                    ! The windowed multipole object to process.
    real(8), intent(in)              :: Emev                         ! The energy at which to evaluate the cross section in MeV
    real(8), intent(in)              :: sqrtkT                       ! The temperature in the form sqrt(kT (in eV)), at which to evaluate the cross section.
    real(8), intent(out)             :: sigT                         ! Total cross section
    real(8), intent(out)             :: sigA                         ! Absorption cross section
    real(8), intent(out)             :: sigF                         ! Fission cross section
    complex(8) :: psi_ki   ! The value of the psi-ki function for the asymptotic form
    complex(8) :: c_temp   ! complex temporary variable 
    complex(8) :: w_val    ! The faddeeva function evaluated at Z
    complex(8) :: Z        ! sqrt(atomic weight ratio / kT) * (sqrt(E) - pole)
    real(8) :: sqrtE       ! sqrt(E), eV
    real(8) :: invE        ! 1/E, eV
    real(8) :: dopp        ! sqrt(atomic weight ratio / kT)
    real(8) :: dopp_ecoef  ! sqrt(atomic weight ratio * pi / kT) / E 
    real(8) :: temp        ! real temporary value
    real(8) :: E           ! energy, eV
    integer :: iP          ! index of pole
    integer :: iC          ! index of curvefit
    integer :: iW          ! index of window
    integer :: startw      ! window start pointer (for poles)
    integer :: startw_1    ! window start pointer - 1
    integer :: startw_endw ! window start pointer - window end pointer
    integer :: endw        ! window end pointer

    ! Convert to eV
    E = Emev * 1.0e6_8

    sqrtE = sqrt(E)
    invE = ONE/E

    if(.not. mp_already_alloc) then
      call multipole_eval_allocate()
    end if

    ! Locate us
    iW = floor((sqrtE - sqrt(multipole % start_E))/multipole % spacing + ONE)

    startw = multipole % w_start(iW)
    startw_1 = startw - 1 ! This is an index shift parameter.
    endw = multipole % w_end(iW)
    startw_endw = endw - startw + 1

    ! Fill in factors
    if (startw <= endw) then
      call fill_factors(multipole, sqrtE, sigT_factor, twophi, multipole % num_l)
    end if

    ! Generate some doppler broadening parameters

    ! dopp_ecoef is inverse of dopp, divided by E, multiplied by sqrt(pi).
    dopp = multipole % sqrtAWR / sqrtKT
    dopp_ecoef = dopp * invE * SQRT_PI

    sigT = ZERO
    sigA = ZERO
    sigF = ZERO
    ! Evaluate linefit first

    if(sqrtkT /= 0 .and. multipole % broaden_poly(iW) == 1) then ! Broaden the curvefit.
      call broaden_n_polynomials(E, dopp, multipole % fit_order + 1, broadened_polynomials)

      do iC = 1, multipole % fit_order+1
        sigT = sigT + multipole % curvefit(FIT_T, iC, iW)*broadened_polynomials(iC)
        sigA = sigA + multipole % curvefit(FIT_A, iC, iW)*broadened_polynomials(iC)
        if (multipole % fissionable) then
          sigF = sigF + multipole % curvefit(FIT_F, iC, iW)*broadened_polynomials(iC)
        end if
      end do
    else ! Evaluate as if it were a polynomial
      temp = invE
      do iC = 1, multipole % fit_order+1

        sigT = sigT + multipole % curvefit(FIT_T, iC, iW)*temp
        sigA = sigA + multipole % curvefit(FIT_A, iC, iW)*temp
        if (multipole % fissionable) then
          sigF = sigF + multipole % curvefit(FIT_F, iC, iW)*temp
        end if

        temp = temp * sqrtE
      end do
    end if

    ! Then get the poles we want and broaden them.

    if (sqrtkT == ZERO) then
      ! If at 0K, use asymptotic form.
      do iP = startw, endw
        psi_ki = -ONEI/(multipole % data(MP_EA, iP) - sqrtE)
        c_temp = psi_ki/E
        if (multipole % formalism == FORM_MLBW) then
          sigT = sigT + real(multipole % data(MLBW_RT, iP) * c_temp * &
                             sigT_factor(multipole % l_value(iP))) &
                      + real(multipole % data(MLBW_RX, iP) * c_temp)
          sigA = sigA + real(multipole % data(MLBW_RA, iP) * c_temp)
          sigF = sigF + real(multipole % data(MLBW_RF, iP) * c_temp)
        else if (multipole % formalism == FORM_RM) then
          sigT = sigT + real(multipole % data(RM_RT, iP) * c_temp* &
                             sigT_factor(multipole % l_value(iP)))
          sigA = sigA + real(multipole % data(RM_RA, iP) * c_temp)
          sigF = sigF + real(multipole % data(RM_RF, iP) * c_temp)
        end if
      end do
    else
      ! At temperature, use Faddeeva function-based form.
      if(endw >= startw) then
        do iP = startw, endw
          Z = (sqrtE - multipole % data(MP_EA, iP)) * dopp
          w_val = w(Z) * dopp_ecoef

          if (multipole % formalism == FORM_MLBW) then
            sigT = sigT + real((multipole % data(MLBW_RT, iP) * &
                          sigT_factor(multipole%l_value(iP)) + &
                          multipole % data(MLBW_RX, iP)) * w_val)
            sigA = sigA + real(multipole % data(MLBW_RA, iP) * w_val)
            sigF = sigF + real(multipole % data(MLBW_RF, iP) * w_val)
          else if (multipole % formalism == FORM_RM) then
            sigT = sigT + real(multipole % data(RM_RT, iP) * w_val * &
                               sigT_factor(multipole % l_value(iP)))
            sigA = sigA + real(multipole % data(RM_RA, iP) * w_val)
            sigF = sigF + real(multipole % data(RM_RF, iP) * w_val)
          end if
        end do
      end if
    end if
  end subroutine

!===============================================================================
! FILL_FACTORS calculates the value of phi, the hardsphere phase shift factor,
! and sigT_factor, a factor inside of the sigT equation not present in the
! sigA and sigF equations.
!===============================================================================

  subroutine fill_factors(multipole, sqrtE, sigT_factor, twophi, max_L)
    type(MultipoleArray), intent(in)   :: multipole
    real(8), intent(in)                  :: sqrtE
    integer, intent(in)                  :: max_L
    complex(8), intent(out)              :: sigT_factor(max_L)
    real(8), intent(out)                 :: twophi(max_L)

    integer :: iL
    real(8) :: arg

    do iL = 1, max_L
      twophi(iL) = multipole % pseudo_k0RS(iL) * sqrtE
      if (iL == 2) then
        twophi(iL) = twophi(iL) - atan(twophi(iL))
      else if (iL == 3) then
        arg = 3.0_8 * twophi(iL) / (3.0_8 - twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      else if (iL == 4) then
        arg = twophi(iL) * (15.0_8 - twophi(iL)**2) / (15.0_8 - 6.0_8 * twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      end if
    end do

    twophi = 2.0_8 * twophi
    sigT_factor = cmplx(cos(twophi),-sin(twophi), KIND=8)
  end subroutine

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
