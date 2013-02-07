module cross_section

  use ace_header,      only: Nuclide, SAlphaBeta, Reaction, UrrData
  use constants
  use error,           only: fatal_error
  use fission,         only: nu_total
  use global
  use material_header, only: Material
  use random_lcg,      only: prn
  use search,          only: binary_search

  implicit none

contains

!===============================================================================
! CALCULATE_XS determines the macroscopic cross sections for the material the
! particle is currently traveling through.
!===============================================================================

  subroutine calculate_xs()

    integer :: i             ! loop index over nuclides
    integer :: i_nuclide     ! index into nuclides array
    integer :: i_sab         ! index into sab_tables array
    integer :: j             ! index in mat % i_sab_nuclides
    real(8) :: atom_density  ! atom density of a nuclide
    logical :: check_sab     ! should we check for S(a,b) table?
    type(Material),  pointer :: mat => null() ! current material

    ! Set all material macroscopic cross sections to zero
    material_xs % total      = ZERO
    material_xs % elastic    = ZERO
    material_xs % absorption = ZERO
    material_xs % fission    = ZERO
    material_xs % nu_fission = ZERO
    material_xs % kappa_fission  = ZERO

    ! Exit subroutine if material is void
    if (p % material == MATERIAL_VOID) return

    mat => materials(p % material)

    ! Find energy index on unionized grid
    if (grid_method == GRID_UNION) call find_energy_index()

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
      if (p % E /= micro_xs(i_nuclide) % last_E) then
        call calculate_nuclide_xs(i_nuclide, i_sab)
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
           
      ! Add contributions to material macroscopic energy release from fission
      material_xs % kappa_fission = material_xs % kappa_fission + &
           atom_density * micro_xs(i_nuclide) % kappa_fission
    end do

  end subroutine calculate_xs

!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array at the energy of the given particle
!===============================================================================

  subroutine calculate_nuclide_xs(i_nuclide, i_sab)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array

    integer :: i_grid ! index on nuclide energy grid
    real(8) :: f      ! interp factor on nuclide energy grid
    type(Nuclide),   pointer :: nuc => null()

    ! Set pointer to nuclide
    nuc => nuclides(i_nuclide)

    ! Determine index on nuclide energy grid
    select case (grid_method)
    case (GRID_UNION)
      ! If we're using the unionized grid with pointers, finding the index on
      ! the nuclide energy grid is as simple as looking up the pointer

      i_grid = nuc % grid_index(p % index_grid)

    case (GRID_NUCLIDE)
      ! If we're not using the unionized grid, we have to do a binary search on
      ! the nuclide energy grid in order to determine which points to
      ! interpolate between

      if (p % E < nuc % energy(1)) then
        i_grid = 1
      elseif (p % E > nuc % energy(nuc % n_grid)) then
        i_grid = nuc % n_grid - 1
      else
        i_grid = binary_search(nuc % energy, nuc % n_grid, p % E)
      end if

    end select

    ! check for rare case where two energy points are the same
    if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) i_grid = i_grid + 1

    ! calculate interpolation factor
    f = (p%E - nuc%energy(i_grid))/(nuc%energy(i_grid+1) - nuc%energy(i_grid))

    micro_xs(i_nuclide) % index_grid    = i_grid
    micro_xs(i_nuclide) % interp_factor = f

    ! Initialize sab treatment to false
    micro_xs(i_nuclide) % index_sab   = NONE
    micro_xs(i_nuclide) % elastic_sab = ZERO
    micro_xs(i_nuclide) % use_ptable  = .false.

    ! Initialize nuclide cross-sections to zero
    micro_xs(i_nuclide) % fission    = ZERO
    micro_xs(i_nuclide) % nu_fission = ZERO
    micro_xs(i_nuclide) % kappa_fission  = ZERO

    ! Calculate microscopic nuclide total cross section
    micro_xs(i_nuclide) % total = (ONE - f) * nuc % total(i_grid) &
         + f * nuc % total(i_grid+1)

    ! Calculate microscopic nuclide total cross section
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
           
      ! Calculate microscopic nuclide kappa-fission cross section
      ! The ENDF standard (ENDF-102) states that MT 18 stores
      ! the fission energy as the Q_value (fission(1))
      micro_xs(i_nuclide) % kappa_fission = &
           nuc % reactions(nuc % index_fission(1)) % Q_value * &
           micro_xs(i_nuclide) % fission
    end if

    ! If there is S(a,b) data for this nuclide, we need to do a few
    ! things. Since the total cross section was based on non-S(a,b) data, we
    ! need to correct it by subtracting the non-S(a,b) elastic cross section and
    ! then add back in the calculated S(a,b) elastic+inelastic cross section.

    if (i_sab > 0) call calculate_sab_xs(i_nuclide, i_sab)

    ! if the particle is in the unresolved resonance range and there are
    ! probability tables, we need to determine cross sections from the table

    if (urr_ptables_on .and. nuc % urr_present) then
      if (p % E > nuc % urr_data % energy(1) .and. &
           p % E < nuc % urr_data % energy(nuc % urr_data % n_energy)) then
        call calculate_urr_xs(i_nuclide)
      end if
    end if

    ! Set last evaluated energy -- if we're in S(a,b) region, force
    ! re-calculation of cross-section
    if (i_sab == 0) then
      micro_xs(i_nuclide) % last_E = p % E
    else
      micro_xs(i_nuclide) % last_E = ZERO
    end if

  end subroutine calculate_nuclide_xs

!===============================================================================
! CALCULATE_SAB_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range. These cross sections replace
! whatever data were taken from the normal Nuclide table.
!===============================================================================

  subroutine calculate_sab_xs(i_nuclide, i_sab)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array

    integer :: i_grid    ! index on S(a,b) energy grid
    real(8) :: f         ! interp factor on S(a,b) energy grid
    real(8) :: inelastic ! S(a,b) inelastic cross section
    real(8) :: elastic   ! S(a,b) elastic cross section
    type(SAlphaBeta), pointer :: sab => null()

    ! Set flag that S(a,b) treatment should be used for scattering
    micro_xs(i_nuclide) % index_sab = i_sab

    ! Get pointer to S(a,b) table
    sab => sab_tables(i_sab)

    ! Get index and interpolation factor for inelastic grid
    if (p%E < sab % inelastic_e_in(1)) then
      i_grid = 1
      f = ZERO
    else
      i_grid = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, p%E)
      f = (p%E - sab%inelastic_e_in(i_grid)) / & 
           (sab%inelastic_e_in(i_grid+1) - sab%inelastic_e_in(i_grid))
    end if

    ! Calculate S(a,b) inelastic scattering cross section
    inelastic = (ONE - f) * sab % inelastic_sigma(i_grid) + &
         f * sab % inelastic_sigma(i_grid + 1)

    ! Check for elastic data
    if (p % E < sab % threshold_elastic) then
      ! Determine whether elastic scattering is given in the coherent or
      ! incoherent approximation. For coherent, the cross section is
      ! represented as P/E whereas for incoherent, it is simply P

      if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
        if (p % E < sab % elastic_e_in(1)) then
          ! If energy is below that of the lowest Bragg peak, the elastic
          ! cross section will be zero
          elastic = ZERO
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, p % E)
          elastic = sab % elastic_P(i_grid) / p % E
        end if
      else
        ! Determine index on elastic energy grid
        if (p % E < sab % elastic_e_in(1)) then
          i_grid = 1
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, p%E)
        end if

        ! Get interpolation factor for elastic grid
        f = (p%E - sab%elastic_e_in(i_grid))/(sab%elastic_e_in(i_grid+1) - &
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

  subroutine calculate_urr_xs(i_nuclide)

    integer, intent(in) :: i_nuclide ! index into nuclides array

    integer :: i_energy   ! index for energy
    integer :: i_table    ! index for table
    real(8) :: f          ! interpolation factor
    real(8) :: r          ! pseudo-random number
    real(8) :: elastic    ! elastic cross section
    real(8) :: capture    ! (n,gamma) cross section
    real(8) :: fission    ! fission cross section
    real(8) :: inelastic  ! inelastic cross section
    type(UrrData),  pointer :: urr => null()
    type(Nuclide),  pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()

    micro_xs(i_nuclide) % use_ptable = .true.

    ! get pointer to probability table
    nuc => nuclides(i_nuclide)
    urr => nuc % urr_data

    ! determine energy table
    i_energy = 1
    do
      if (p % E < urr % energy(i_energy + 1)) exit
      i_energy = i_energy + 1
    end do

    ! determine interpolation factor on table
    f = (p % E - urr % energy(i_energy)) / &
         (urr % energy(i_energy + 1) - urr % energy(i_energy))

    ! sample probability table using the cumulative distribution
    r = prn()
    i_table = 1
    do
      if (urr % prob(i_energy, URR_CUM_PROB, i_table) > r) exit
      i_table = i_table + 1
    end do

    ! determine elastic, fission, and capture cross sections from probability
    ! table
    if (urr % interp == LINEAR_LINEAR) then
      elastic = (ONE - f) * urr % prob(i_energy, URR_ELASTIC, i_table) + &
           f * urr % prob(i_energy + 1, URR_ELASTIC, i_table)
      fission = (ONE - f) * urr % prob(i_energy, URR_FISSION, i_table) + &
           f * urr % prob(i_energy + 1, URR_FISSION, i_table)
      capture = (ONE - f) * urr % prob(i_energy, URR_N_GAMMA, i_table) + &
           f * urr % prob(i_energy + 1, URR_N_GAMMA, i_table)
    elseif (urr % interp == LOG_LOG) then
      ! Get logarithmic interpolation factor
      f = log(p % E / urr % energy(i_energy)) / &
           log(urr % energy(i_energy + 1) / urr % energy(i_energy))

      ! Calculate elastic cross section/factor
      elastic = ZERO
      if (urr % prob(i_energy, URR_ELASTIC, i_table) > ZERO) then
        elastic = exp((ONE - f) * log(urr % prob(i_energy, URR_ELASTIC, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_ELASTIC, &
             i_table)))
      end if

      ! Calculate fission cross section/factor
      fission = ZERO
      if (urr % prob(i_energy, URR_FISSION, i_table) > ZERO) then
        fission = exp((ONE - f) * log(urr % prob(i_energy, URR_FISSION, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_FISSION, &
             i_table)))
      end if

      ! Calculate capture cross section/factor
      capture = ZERO
      if (urr % prob(i_energy, URR_N_GAMMA, i_table) > ZERO) then
        capture = exp((ONE - f) * log(urr % prob(i_energy, URR_N_GAMMA, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_N_GAMMA, &
             i_table)))
      end if
    end if

    ! Determine treatment of inelastic scattering
    inelastic = ZERO
    if (urr % inelastic_flag > 0) then
      ! Get pointer to inelastic scattering reaction
      rxn => nuc % reactions(nuc % urr_inelastic)

      ! Get index on energy grid and interpolation factor
      i_energy = micro_xs(i_nuclide) % index_grid
      f = micro_xs(i_nuclide) % interp_factor

      ! Determine inelastic scattering cross section
      if (i_energy >= rxn % threshold) then
        inelastic = (ONE - f) * rxn % sigma(i_energy - rxn%threshold + 1) + &
             f * rxn % sigma(i_energy - rxn%threshold + 2)
      end if
    end if

    ! Multiply by smooth cross-section if needed
    if (urr % multiply_smooth) then
      elastic = elastic * micro_xs(i_nuclide) % elastic
      capture = capture * (micro_xs(i_nuclide) % absorption - &
           micro_xs(i_nuclide) % fission)
      fission = fission * micro_xs(i_nuclide) % fission
    end if

    ! Set elastic, absorption, fission, and total cross sections. Note that the
    ! total cross section is calculated as sum of partials rather than using the
    ! table-provided value
    micro_xs(i_nuclide) % elastic = elastic
    micro_xs(i_nuclide) % absorption = capture + fission
    micro_xs(i_nuclide) % fission = fission
    micro_xs(i_nuclide) % total = elastic + inelastic + capture + fission

    ! Determine nu-fission cross section
    if (nuc % fissionable) then
      micro_xs(i_nuclide) % nu_fission = nu_total(nuc, p % E) * &
           micro_xs(i_nuclide) % fission
    end if

  end subroutine calculate_urr_xs

!===============================================================================
! FIND_ENERGY_INDEX determines the index on the union energy grid and the
! interpolation factor for a particle at a certain energy
!===============================================================================

  subroutine find_energy_index()

    integer :: i_grid ! index on union energy grid
    real(8) :: E      ! energy of particle
    real(8) :: interp ! interpolation factor

    ! copy particle's energy
    E = p % E

    ! if particle's energy is outside of energy grid range, set to first or last
    ! index. Otherwise, do a binary search through the union energy grid.
    if (E < e_grid(1)) then
      i_grid = 1
    elseif (E > e_grid(n_grid)) then
      i_grid = n_grid - 1
    else
      i_grid = binary_search(e_grid, n_grid, E)
    end if

    ! calculate the interpolation factor -- note this will be outside of [0,1)
    ! for a particle outside the energy range of the union grid
    interp = (E - e_grid(i_grid))/(e_grid(i_grid+1) - e_grid(i_grid))

    ! set particle attributes
    p % index_grid = i_grid
    p % interp     = interp

  end subroutine find_energy_index

end module cross_section
