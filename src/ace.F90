module ace

  use angleenergy_header, only: AngleEnergy
  use constants
  use distribution_univariate, only: Uniform, Equiprobable, Tabular
  use endf, only: is_fission, is_disappearance
  use endf_header, only: Constant1D, Tabulated1D, Polynomial
  use energy_distribution, only: TabularEquiprobable, LevelInelastic, &
       ContinuousTabular, MaxwellEnergy, Evaporation, WattEnergy
  use error, only: fatal_error, warning
  use global
  use list_header, only: ListInt
  use material_header, only: Material
  use nuclide_header
  use output, only: write_message
  use product_header, only: ReactionProduct
  use sab_header
  use set_header, only: SetChar
  use secondary_correlated, only: CorrelatedAngleEnergy
  use secondary_kalbach, only: KalbachMann
  use secondary_nbody, only: NBodyPhaseSpace
  use secondary_uncorrelated, only: UncorrelatedAngleEnergy
  use string, only: to_str, to_lower

  implicit none

  integer                   :: JXS(32)   ! Pointers into ACE XSS tables
  integer                   :: NXS(16)   ! Descriptors for ACE XSS tables
  real(8), allocatable      :: XSS(:)    ! Cross section data
  integer                   :: XSS_index ! Current index in XSS data

  private :: JXS
  private :: NXS
  private :: XSS

contains

!===============================================================================
! READ_ACE_XS reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_ace_xs()

    integer :: i            ! index in materials array
    integer :: j            ! index over nuclides in material
    integer :: k            ! index over S(a,b) tables in material
    integer :: n            ! index over resonant scatterers
    integer :: i_listing    ! index in xs_listings array
    integer :: i_nuclide    ! index in nuclides
    integer :: i_sab        ! index in sab_tables
    integer :: m            ! position for sorting
    integer :: temp_nuclide ! temporary value for sorting
    integer :: temp_table   ! temporary value for sorting
    character(12)  :: name  ! name of isotope, e.g. 92235.03c
    character(12)  :: alias ! alias of nuclide, e.g. U-235.03c
    type(Material),   pointer :: mat
    type(Nuclide), pointer :: nuc
    type(SAlphaBeta), pointer :: sab
    type(SetChar) :: already_read

    ! allocate arrays for ACE table storage and cross section cache
    allocate(nuclides(n_nuclides_total))
    allocate(sab_tables(n_sab_tables))
!$omp parallel
    allocate(micro_xs(n_nuclides_total))
!$omp end parallel

    ! ==========================================================================
    ! READ ALL ACE CROSS SECTION TABLES

    ! Loop over all files
    MATERIAL_LOOP: do i = 1, n_materials
      mat => materials(i)

      NUCLIDE_LOOP: do j = 1, mat % n_nuclides
        name = mat % names(j)

        if (.not. already_read % contains(name)) then
          i_listing = xs_listing_dict % get_key(to_lower(name))
          i_nuclide = nuclide_dict % get_key(to_lower(name))
          name  = xs_listings(i_listing) % name
          alias = xs_listings(i_listing) % alias

          ! Keep track of what listing is associated with this nuclide
          nuc => nuclides(i_nuclide)
          nuc % listing = i_listing

          ! Read the ACE table into the appropriate entry on the nuclides
          ! array
          call read_ace_table(i_nuclide, i_listing)

          ! 0K resonant scatterer information, if treating resonance scattering
          if (treat_res_scat) then
            do n = 1, n_res_scatterers_total
              if (name == nuclides_0K(n) % name) then
                nuclides(i_nuclide) % resonant = .true.
                nuclides(i_nuclide) % name_0K = nuclides_0K(n) % name_0K
                nuclides(i_nuclide) % name_0K = trim(nuclides(i_nuclide) % &
                     & name_0K)
                nuclides(i_nuclide) % scheme = nuclides_0K(n) % scheme
                nuclides(i_nuclide) % scheme = trim(nuclides(i_nuclide) % &
                     & scheme)
                nuclides(i_nuclide) % E_min = nuclides_0K(n) % E_min
                nuclides(i_nuclide) % E_max = nuclides_0K(n) % E_max
                if (.not. already_read % contains(nuclides(i_nuclide) % &
                     & name_0K)) then
                  i_listing = xs_listing_dict % get_key(nuclides(i_nuclide) % &
                       & name_0K)
                  call read_ace_table(i_nuclide, i_listing)
                end if
                exit
              end if
            end do
          end if

          ! Add name and alias to dictionary
          call already_read % add(name)
          call already_read % add(alias)
        end if
      end do NUCLIDE_LOOP

      SAB_LOOP: do k = 1, mat % n_sab
        ! Get name of S(a,b) table
        name = mat % sab_names(k)

        if (.not. already_read % contains(name)) then
          i_listing = xs_listing_dict % get_key(to_lower(name))
          i_sab  = sab_dict % get_key(to_lower(name))

          ! Read the ACE table into the appropriate entry on the sab_tables
          ! array
          call read_ace_table(i_sab, i_listing)

          ! Add name to dictionary
          call already_read % add(name)
        end if
      end do SAB_LOOP
    end do MATERIAL_LOOP

    ! ==========================================================================
    ! ASSIGN S(A,B) TABLES TO SPECIFIC NUCLIDES WITHIN MATERIALS

    MATERIAL_LOOP2: do i = 1, n_materials
      ! Get pointer to material
      mat => materials(i)

      ASSIGN_SAB: do k = 1, mat % n_sab
        ! In order to know which nuclide the S(a,b) table applies to, we need to
        ! search through the list of nuclides for one which has a matching zaid
        sab => sab_tables(mat % i_sab_tables(k))

        ! Loop through nuclides and find match
        FIND_NUCLIDE: do j = 1, mat % n_nuclides
          if (any(sab % zaid == nuclides(mat % nuclide(j)) % zaid)) then
            mat % i_sab_nuclides(k) = j
            exit FIND_NUCLIDE
          end if
        end do FIND_NUCLIDE

        ! Check to make sure S(a,b) table matched a nuclide
        if (mat % i_sab_nuclides(k) == NONE) then
          call fatal_error("S(a,b) table " // trim(mat % sab_names(k)) &
               &// " did not match any nuclide on material " &
               &// trim(to_str(mat % id)))
        end if
      end do ASSIGN_SAB

      ! If there are multiple S(a,b) tables, we need to make sure that the
      ! entries in i_sab_nuclides are sorted or else they won't be applied
      ! correctly in the cross_section module. The algorithm here is a simple
      ! insertion sort -- don't need anything fancy!

      if (mat % n_sab > 1) then
        SORT_SAB: do k = 2, mat % n_sab
          ! Save value to move
          m = k
          temp_nuclide = mat % i_sab_nuclides(k)
          temp_table   = mat % i_sab_tables(k)

          MOVE_OVER: do
            ! Check if insertion value is greater than (m-1)th value
            if (temp_nuclide >= mat % i_sab_nuclides(m-1)) exit

            ! Move values over until hitting one that's not larger
            mat % i_sab_nuclides(m) = mat % i_sab_nuclides(m-1)
            mat % i_sab_tables(m)   = mat % i_sab_tables(m-1)
            m = m - 1

            ! Exit if we've reached the beginning of the list
            if (m == 1) exit
          end do MOVE_OVER

          ! Put the original value into its new position
          mat % i_sab_nuclides(m) = temp_nuclide
          mat % i_sab_tables(m)   = temp_table
        end do SORT_SAB
      end if

      ! Deallocate temporary arrays for names of nuclides and S(a,b) tables
      if (allocated(mat % names)) deallocate(mat % names)

    end do MATERIAL_LOOP2

    ! Avoid some valgrind leak errors
    call already_read % clear()

    ! Loop around material
    MATERIAL_LOOP3: do i = 1, n_materials

      ! Get material
      mat => materials(i)

      ! Loop around nuclides in material
      NUCLIDE_LOOP2: do j = 1, mat % n_nuclides

        ! Check for fission in nuclide
        if (nuclides(mat % nuclide(j)) % fissionable) then
          mat % fissionable = .true.
          exit NUCLIDE_LOOP2
        end if

      end do NUCLIDE_LOOP2

    end do MATERIAL_LOOP3

    ! Show which nuclide results in lowest energy for neutron transport
    do i = 1, n_nuclides_total
      if (nuclides(i) % energy(nuclides(i) % n_grid) == energy_max_neutron) then
        call write_message("Maximum neutron transport energy: " // &
             trim(to_str(energy_max_neutron)) // " MeV for " // &
             trim(adjustl(nuclides(i) % name)), 6)
        exit
      end if
    end do

  end subroutine read_ace_xs

!===============================================================================
! READ_ACE_TABLE reads a single cross section table in either ASCII or binary
! format. This routine reads the header data for each table and then calls
! appropriate subroutines to parse the actual data.
!===============================================================================

  subroutine read_ace_table(i_table, i_listing)
    integer, intent(in) :: i_table   ! index in nuclides/sab_tables
    integer, intent(in) :: i_listing ! index in xs_listings

    integer       :: i             ! loop index for XSS records
    integer       :: j, j1, j2     ! indices in XSS
    integer       :: record_length ! Fortran record length
    integer       :: location      ! location of ACE table
    integer       :: entries       ! number of entries on each record
    integer       :: length        ! length of ACE table
    integer       :: unit_ace      ! file unit
    integer       :: zaids(16)     ! list of ZAIDs (only used for S(a,b))
    integer       :: filetype      ! filetype (ASCII or BINARY)
    real(8)       :: kT            ! temperature of table
    real(8)       :: awrs(16)      ! list of atomic weight ratios (not used)
    real(8)       :: awr           ! atomic weight ratio for table
    logical       :: file_exists   ! does ACE library exist?
    logical       :: data_0K       ! are we reading 0K data?
    character(7)  :: readable      ! is ACE library readable?
    character(10) :: name          ! name of ACE table
    character(10) :: date_         ! date ACE library was processed
    character(10) :: mat           ! material identifier
    character(70) :: comment       ! comment for ACE table
    character(MAX_FILE_LEN) :: filename ! path to ACE cross section library
    type(Nuclide), pointer :: nuc
    type(SAlphaBeta), pointer :: sab
    type(XsListing),  pointer :: listing

    ! determine path, record length, and location of table
    listing => xs_listings(i_listing)
    filename      = listing % path
    record_length = listing % recl
    location      = listing % location
    entries       = listing % entries
    filetype      = listing % filetype

    ! Check if ACE library exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
      call fatal_error("ACE library '" // trim(filename) // "' does not exist!")
    elseif (readable(1:3) == 'NO') then
      call fatal_error("ACE library '" // trim(filename) // "' is not readable!&
           & Change file permissions with chmod command.")
    end if

    ! display message
    call write_message("Loading ACE cross section table: " // listing % name, 6)

    if (filetype == ASCII) then
      ! =======================================================================
      ! READ ACE TABLE IN ASCII FORMAT

      ! Find location of table
      open(NEWUNIT=unit_ace, FILE=filename, STATUS='old', ACTION='read')
      rewind(UNIT=unit_ace)
      do i = 1, location - 1
        read(UNIT=unit_ace, FMT=*)
      end do

      ! Read first line of header
      read(UNIT=unit_ace, FMT='(A10,2G12.0,1X,A10)') name, awr, kT, date_

      ! Check that correct xs was found -- if cross_sections.xml is broken, the
      ! location of the table may be wrong
      if(adjustl(name) /= adjustl(listing % name)) then
        call fatal_error("XS listing entry " // trim(listing % name) // " did &
             &not match ACE data, " // trim(name) // " found instead.")
      end if

      ! Read more header and NXS and JXS
      read(UNIT=unit_ace, FMT=100) comment, mat, &
           (zaids(i), awrs(i), i=1,16), NXS, JXS
100   format(A70,A10/4(I7,F11.0)/4(I7,F11.0)/4(I7,F11.0)/4(I7,F11.0)/&
           ,8I9/8I9/8I9/8I9/8I9/8I9)

      ! determine table length
      length = NXS(1)
      allocate(XSS(length))

      ! Read XSS array
      read(UNIT=unit_ace, FMT='(4G20.0)') XSS

      ! Close ACE file
      close(UNIT=unit_ace)

    elseif (filetype == BINARY) then
      ! =======================================================================
      ! READ ACE TABLE IN BINARY FORMAT

      ! Open ACE file
      open(NEWUNIT=unit_ace, FILE=filename, STATUS='old', ACTION='read', &
           ACCESS='direct', RECL=record_length)

      ! Read all header information
      read(UNIT=unit_ace, REC=location) name, awr, kT, date_, &
           comment, mat, (zaids(i), awrs(i), i=1,16), NXS, JXS

      ! determine table length
      length = NXS(1)
      allocate(XSS(length))

      ! Read remaining records with XSS
      do i = 1, (length + entries - 1)/entries
        j1 = 1 + (i-1)*entries
        j2 = min(length, j1 + entries - 1)
        read(UNIT=UNIT_ACE, REC=location + i) (XSS(j), j=j1,j2)
      end do

      ! Close ACE file
      close(UNIT=unit_ace)
    end if

    ! ==========================================================================
    ! PARSE DATA BASED ON NXS, JXS, AND XSS ARRAYS

    select case(listing % type)
    case (ACE_NEUTRON)

      ! only read in a resonant scatterers info once
      nuc => nuclides(i_table)
      data_0K = .false.
      if (trim(adjustl(name)) == nuc % name_0K) then
        data_0K = .true.
      else
        nuc % name = name
        nuc % awr  = awr
        nuc % kT   = kT
        nuc % zaid = listing % zaid
      end if

      ! read all blocks
      call read_esz(nuc, data_0K)

      ! don't read unnecessary 0K data for resonant scatterers
      if (data_0K) then
        continue
      else
        call read_reactions(nuc)
        call read_nu_data(nuc)
        call read_energy_dist(nuc)
        call read_angular_dist(nuc)
        call read_unr_res(nuc)
      end if

      ! for fissionable nuclides, precalculate microscopic nu-fission cross
      ! sections so that we don't need to call the nu_total function during
      ! cross section lookups (except if we're dealing w/ 0K data for resonant
      ! scatterers)

      if (nuc % fissionable .and. .not. data_0K) then
        call generate_nu_fission(nuc)
      end if

    case (ACE_THERMAL)
      sab => sab_tables(i_table)
      sab % name = name
      sab % awr  = awr
      sab % kT   = kT
      ! Find sab % n_zaid
      do i = 1, 16
        if (zaids(i) == 0) then
          sab % n_zaid = i - 1
          exit
        end if
      end do
      allocate(sab % zaid(sab % n_zaid))
      sab % zaid = zaids(1: sab % n_zaid)

      call read_thermal_data(sab)
    end select

    deallocate(XSS)

  end subroutine read_ace_table

!===============================================================================
! READ_ESZ - reads through the ESZ block. This block contains the energy grid,
! total xs, absorption xs, elastic scattering xs, and heating numbers.
!===============================================================================

  subroutine read_esz(nuc, data_0K)
    type(Nuclide), intent(inout) :: nuc
    logical,          intent(in)    :: data_0K ! are we reading 0K data?

    integer :: NE ! number of energy points for total and elastic cross sections
    integer :: i  ! index in 0K elastic xs array for this nuclide

    real(8) :: xs_cdf_sum = ZERO ! xs cdf value

    ! determine number of energy points
    NE = NXS(3)

    ! allocate storage for energy grid and cross section arrays

    ! read in 0K data if we've already read in non-0K data
    if (data_0K) then
      nuc % n_grid_0K = NE
      allocate(nuc % energy_0K(NE))
      allocate(nuc % elastic_0K(NE))
      allocate(nuc % xs_cdf(NE))
      nuc % elastic_0K = ZERO
      nuc % xs_cdf = ZERO
      XSS_index = 1
      nuc % energy_0K = get_real(NE)

      ! Skip total and absorption
      XSS_index = XSS_index + 2*NE

      ! Continue reading elastic scattering and heating
      nuc % elastic_0K = get_real(NE)

      do i = 1, nuc % n_grid_0K - 1

        ! Negative cross sections result in a CDF that is not monotonically
        ! increasing. Set all negative xs values to ZERO.
        if (nuc % elastic_0K(i) < ZERO) nuc % elastic_0K(i) = ZERO

        ! build xs cdf
        xs_cdf_sum = xs_cdf_sum &
             + (sqrt(nuc % energy_0K(i)) * nuc % elastic_0K(i) &
             + sqrt(nuc % energy_0K(i+1)) * nuc % elastic_0K(i+1)) / TWO &
             * (nuc % energy_0K(i+1) - nuc % energy_0K(i))
        nuc % xs_cdf(i) = xs_cdf_sum
      end do

    else ! read in non-0K data
      nuc % n_grid = NE
      allocate(nuc % energy(NE))
      allocate(nuc % total(NE))
      allocate(nuc % elastic(NE))
      allocate(nuc % fission(NE))
      allocate(nuc % nu_fission(NE))
      allocate(nuc % absorption(NE))

      ! initialize cross sections
      nuc % total      = ZERO
      nuc % elastic    = ZERO
      nuc % fission    = ZERO
      nuc % nu_fission = ZERO
      nuc % absorption = ZERO

      ! Read data from XSS -- only the energy grid, elastic scattering and heating
      ! cross section values are actually read from here. The total and absorption
      ! cross sections are reconstructed from the partial reaction data.

      XSS_index = 1
      nuc % energy = get_real(NE)

      ! Skip total and absorption
      XSS_index = XSS_index + 2*NE

      ! Continue reading elastic scattering and heating
      nuc % elastic = get_real(NE)

      ! Determine if minimum/maximum energy for this nuclide is greater/less
      ! than the previous
      energy_min_neutron = max(energy_min_neutron, nuc%energy(1))
      energy_max_neutron = min(energy_max_neutron, nuc%energy(NE))
    end if

  end subroutine read_esz

!===============================================================================
! READ_NU_DATA reads data given on the number of neutrons emitted from fission
! as a function of the incoming energy of a neutron. This data may be broken
! down into prompt and delayed neutrons emitted as well.
!===============================================================================

  subroutine read_nu_data(nuc)
    type(Nuclide), intent(inout) :: nuc

    integer :: i, j   ! loop index
    integer :: idx    ! index in XSS
    integer :: KNU    ! location for nu data
    integer :: LNU    ! type of nu data (polynomial or tabular)
    integer :: NR     ! number of interpolation regions
    integer :: NE     ! number of energies
    integer :: NPCR   ! number of delayed neutron precursor groups
    integer :: LOCC   ! location of energy distributions for given MT
    integer :: LAW
    integer :: IDAT
    real(8) :: total_group_probability
    type(Tabulated1D) :: yield_delayed
    type(Tabulated1D) :: group_probability

    if (JXS(2) == 0) then
      ! Nuclide is not fissionable
      return
    end if

    ! Determine number of delayed neutron precursors
    if (JXS(24) > 0) then
      NPCR = NXS(8)
    else
      NPCR = 0
    end if
    nuc % n_precursor = NPCR

    ! Check to make sure nuclide does not have more than the maximum number
    ! of delayed groups
    if (NPCR > MAX_DELAYED_GROUPS) then
      call fatal_error("Encountered nuclide with " // trim(to_str(NPCR)) &
           // " delayed groups while the maximum number of delayed groups is " &
           // trim(to_str(MAX_DELAYED_GROUPS)))
    end if

    associate (rx => nuc % reactions(nuc % index_fission(1)))
      ! Allocate space for prompt/delayed neutron products
      allocate(rx % products(1 + NPCR))
      rx % products(:) % particle = NEUTRON

      if (XSS(JXS(2)) > 0) then
        ! =======================================================================
        ! PROMPT OR TOTAL NU DATA

        ! If delayed data is present, then prompt data must be present. Otherwise
        ! the product represents 'total' neutron emission
        if (JXS(24) > 0) then
          rx % products(1) % emission_mode = EMISSION_PROMPT
        else
          rx % products(1) % emission_mode = EMISSION_TOTAL
        end if

        KNU = JXS(2)
        LNU = nint(XSS(KNU))
        if (LNU == 1) then
          ! Polynomial data
          allocate(Polynomial :: rx % products(1) % yield)

          ! determine order of polynomial and read coefficients
          select type (yield => rx % products(1) % yield)
          type is (Polynomial)
            call yield % from_ace(XSS, KNU + 1)
          end select

        elseif (LNU == 2) then
          ! Tabulated data
          allocate(Tabulated1D :: rx % products(1) % yield)

          select type(yield => rx % products(1) % yield)
          type is (Tabulated1D)
            call yield % from_ace(XSS, KNU + 1)
          end select

        end if

      elseif (XSS(JXS(2)) < 0) then
        ! =======================================================================
        ! PROMPT AND TOTAL NU DATA

        rx % products(1) % emission_mode = EMISSION_PROMPT

        KNU = JXS(2) + 1
        LNU = nint(XSS(KNU))
        if (LNU == 1) then
          ! Polynomial data
          allocate(Polynomial :: rx % products(1) % yield)

          ! determine order of polynomial and read coefficients
          select type (yield => rx % products(1) % yield)
          type is (Polynomial)
            call yield % from_ace(XSS, KNU + 1)
          end select

        elseif (LNU == 2) then
          ! Tabulated data
          allocate(Tabulated1D :: rx % products(1) % yield)

          select type(yield => rx % products(1) % yield)
          type is (Tabulated1D)
            call yield % from_ace(XSS, KNU + 1)
          end select
        end if

        KNU = JXS(2) + nint(abs(XSS(JXS(2)))) + 1
        LNU = nint(XSS(KNU))
        if (LNU == 1) then
          ! Polynomial data
          allocate(Polynomial :: nuc % total_nu)

          ! determine order of polynomial and read coefficients
          select type (yield => nuc % total_nu)
          type is (Polynomial)
            call yield % from_ace(XSS, KNU + 1)
          end select

        elseif (LNU == 2) then
          ! Tabulated data
          allocate(Tabulated1D :: nuc % total_nu)

          select type(yield => nuc % total_nu)
          type is (Tabulated1D)
            call yield % from_ace(XSS, KNU + 1)
          end select
        end if
      end if

      if (JXS(24) > 0) then
        ! =======================================================================
        ! DELAYED NU DATA

        ! Read total yield of delayed neutrons
        call yield_delayed % from_ace(XSS, JXS(24) + 1)

        idx = JXS(25)
        total_group_probability = ZERO
        do i = 1, NPCR
          ! Set emission mode and decay rate
          rx % products(1 + i) % emission_mode = EMISSION_DELAYED
          rx % products(1 + i) % decay_rate = XSS(idx)

          ! Read probability for this precursor group
          call group_probability % from_ace(XSS, idx + 1)

          ! Set yield based on product of group probability and delayed yield
          if (all(group_probability % y == group_probability % y(1))) then
            allocate(Tabulated1D :: rx % products(1 + i) % yield)
            select type (yield => rx % products(1 + i) % yield)
            type is (Tabulated1D)
              yield = yield_delayed
              yield % y(:) = yield % y(:) * group_probability % y(1)
              total_group_probability = total_group_probability + group_probability % y(1)
            end select
          else
            call fatal_error("Delayed neutron with energy-dependent group &
                 &probability not implemented")
          end if

          ! Advance position
          NR = nint(XSS(idx + 1))
          NE = nint(XSS(idx + 2 + 2*NR))
          idx = idx + 3 + 2*(NR + NE)

          ! =======================================================================
          ! DELAYED NEUTRON ENERGY DISTRIBUTION

          ! Read energy distribution
          LOCC = nint(XSS(JXS(26) + i - 1))

          ! Determine law and location of data
          LAW = nint(XSS(JXS(27) + LOCC))
          IDAT = nint(XSS(JXS(27) + LOCC + 1))

          ! read energy distribution data
          associate(p => rx % products(1 + i))
            allocate(p % applicability(1))
            allocate(p % distribution(1))
            call get_energy_dist(p % distribution(1) % obj, LAW, JXS(27), IDAT, &
                 ZERO, ZERO)

            select type (aedist => p % distribution(1) % obj)
            type is (UncorrelatedAngleEnergy)
              aedist % fission = .true.
            end select
          end associate
        end do

        ! Renormalize delayed neutron yields to reflect fact that in ACE file, the
        ! sum of the group probabilities is not exactly one
        do i = 1, NPCR
          select type (yield => rx % products(1 + i) % yield)
          type is (Tabulated1D)
            yield % y(:) = yield % y(:) / total_group_probability
          end select
        end do
      end if

      ! Assign products to other fission reactions
      do i = 2, nuc % n_fission
        j = nuc % index_fission(i)
        allocate(nuc % reactions(j) % products(1 + NPCR))
        nuc % reactions(j) % products(:) = rx % products(:)
      end do
    end associate

  end subroutine read_nu_data

!===============================================================================
! READ_REACTIONS - Get the list of reaction MTs for this cross-section
! table. The MT values are somewhat arbitrary. Also read in Q-values, neutron
! multiplicities, and cross-sections.
!===============================================================================

  subroutine read_reactions(nuc)
    type(Nuclide), intent(inout) :: nuc

    integer :: i         ! loop indices
    integer :: i_fission ! index in nuc % index_fission
    integer :: LMT       ! index of MT list in XSS
    integer :: NMT       ! Number of reactions
    integer :: JXS4      ! index of Q values in XSS
    integer :: JXS5      ! index of neutron multiplicities in XSS
    integer :: JXS7      ! index of reactions cross-sections in XSS
    integer :: LXS       ! location of cross-section locators
    integer :: LOCA      ! location of cross-section for given MT
    integer :: IE        ! reaction's starting index on energy grid
    integer :: NE        ! number of energies
    real(8) :: y
    type(ListInt) :: MTs

    LMT  = JXS(3)
    JXS4 = JXS(4)
    JXS5 = JXS(5)
    LXS  = JXS(6)
    JXS7 = JXS(7)
    NMT  = NXS(4)

    ! allocate array of reactions. Add one since we need to include an elastic
    ! scattering channel
    nuc % n_reaction = NMT + 1
    allocate(nuc % reactions(NMT+1))

    ! Store elastic scattering cross-section on reaction one -- note that the
    ! sigma array is not allocated or stored for elastic scattering since it is
    ! already stored in nuc % elastic
    associate (rxn => nuc % reactions(1))
      rxn % MT = 2
      rxn % Q_value = ZERO
      allocate(rxn % products(1))
      rxn % products(1) % particle = NEUTRON
      allocate(Constant1D :: rxn % products(1) % yield)
      select type(yield => rxn % products(1) % yield)
      type is (Constant1D)
        yield % y = 1
      end select
      rxn % threshold = 1
      rxn % scatter_in_cm = .true.
      allocate(rxn % products(1) % distribution(1))
      allocate(UncorrelatedAngleEnergy :: rxn % products(1) % distribution(1) % obj)
    end associate

    ! Add contribution of elastic scattering to total cross section
    nuc % total = nuc % total + nuc % elastic

    ! By default, set nuclide to not fissionable and then change if fission
    ! reactions are encountered
    nuc % fissionable = .false.
    nuc % has_partial_fission = .false.
    nuc % n_fission = 0
    i_fission = 0

    do i = 1, NMT
      associate (rxn => nuc % reactions(i+1))
        ! read MT number, Q-value, and neutrons produced
        rxn % MT = int(XSS(LMT + i - 1))
        rxn % Q_value = XSS(JXS4 + i - 1)
        rxn % scatter_in_cm = (nint(XSS(JXS5 + i - 1)) < 0)

        if (.not. is_fission(rxn % MT)) then
          allocate(rxn % products(1))
          rxn % products(1) % particle = NEUTRON

          y = abs(nint(XSS(JXS5 + i - 1)))
          if (y > 100) then
            ! Read energy-dependent multiplicities

            ! Set flag and allocate space for Tabulated1D to store yield
            allocate(Tabulated1D :: rxn % products(1) % yield)

            ! Read yield function
            select type (yield => rxn % products(1) % yield)
            type is (Tabulated1D)
              XSS_index = JXS(11) + int(y) - 101
              call yield % from_ace(XSS, XSS_index)
            end select
          else
            ! Integral yield
            allocate(Constant1D :: rxn % products(1) % yield)
            select type (yield => rxn % products(1) % yield)
            type is (Constant1D)
              yield % y = y
            end select
          end if
        end if

        ! read starting energy index
        LOCA = int(XSS(LXS + i - 1))
        IE   = int(XSS(JXS7 + LOCA - 1))
        rxn % threshold = IE

        ! read number of energies cross section values
        NE = int(XSS(JXS7 + LOCA))
        allocate(rxn % sigma(NE))
        XSS_index = JXS7 + LOCA + 1
        rxn % sigma = get_real(NE)
      end associate
    end do

    ! Create set of MT values
    do i = 1, size(nuc % reactions)
      call MTs % append(nuc % reactions(i) % MT)
      call nuc%reaction_index%add_key(nuc%reactions(i)%MT, i)
    end do

    ! Create total, absorption, and fission cross sections
    do i = 2, size(nuc % reactions)
      associate (rxn => nuc % reactions(i))
        IE = rxn % threshold
        NE = size(rxn % sigma)

        ! Skip total inelastic level scattering, gas production cross sections
        ! (MT=200+), etc.
        if (rxn % MT == N_LEVEL) cycle
        if (rxn % MT > N_5N2P .and. rxn % MT < N_P0) cycle

        ! Skip level cross sections if total is available
        if (rxn % MT >= N_P0 .and. rxn % MT <= N_PC .and. MTs % contains(N_P)) cycle
        if (rxn % MT >= N_D0 .and. rxn % MT <= N_DC .and. MTs % contains(N_D)) cycle
        if (rxn % MT >= N_T0 .and. rxn % MT <= N_TC .and. MTs % contains(N_T)) cycle
        if (rxn % MT >= N_3HE0 .and. rxn % MT <= N_3HEC .and. MTs % contains(N_3HE)) cycle
        if (rxn % MT >= N_A0 .and. rxn % MT <= N_AC .and. MTs % contains(N_A)) cycle
        if (rxn % MT >= N_2N0 .and. rxn % MT <= N_2NC .and. MTs % contains(N_2N)) cycle

        ! Add contribution to total cross section
        nuc % total(IE:IE+NE-1) = nuc % total(IE:IE+NE-1) + rxn % sigma

        ! Add contribution to absorption cross section
        if (is_disappearance(rxn % MT)) then
          nuc % absorption(IE:IE+NE-1) = nuc % absorption(IE:IE+NE-1) + rxn % sigma
        end if

        ! Information about fission reactions
        if (rxn % MT == N_FISSION) then
          allocate(nuc % index_fission(1))
        elseif (rxn % MT == N_F) then
          allocate(nuc % index_fission(PARTIAL_FISSION_MAX))
          nuc % has_partial_fission = .true.
        end if

        ! Add contribution to fission cross section
        if (is_fission(rxn % MT)) then
          nuc % fissionable = .true.
          nuc % fission(IE:IE+NE-1) = nuc % fission(IE:IE+NE-1) + rxn % sigma

          ! Also need to add fission cross sections to absorption
          nuc % absorption(IE:IE+NE-1) = nuc % absorption(IE:IE+NE-1) + rxn % sigma

          ! If total fission reaction is present, there's no need to store the
          ! reaction cross-section since it was copied to nuc % fission
          if (rxn % MT == N_FISSION) deallocate(rxn % sigma)

          ! Keep track of this reaction for easy searching later
          i_fission = i_fission + 1
          nuc % index_fission(i_fission) = i
          nuc % n_fission = nuc % n_fission + 1
        end if
      end associate
    end do

    ! Clear MTs set
    call MTs % clear()

  end subroutine read_reactions

!===============================================================================
! READ_ANGULAR_DIST parses the angular distribution for each reaction with
! secondary neutrons
!===============================================================================

  subroutine read_angular_dist(nuc)
    type(Nuclide), intent(inout) :: nuc

    integer :: LOCB   ! location of angular distribution for given MT
    integer :: NE     ! number of incoming energies
    integer :: NP     ! number of points for cosine distribution
    integer :: i      ! index in reactions array
    integer :: j      ! index over incoming energies
    integer :: k      ! index over energy distributions
    integer :: interp
    integer, allocatable :: LC(:) ! locator

    ! loop over all reactions with secondary neutrons -- NXS(5) does not include
    ! elastic scattering
    do i = 1, NXS(5) + 1
      associate (rxn => nuc%reactions(i))
        ! find location of angular distribution
        LOCB = int(XSS(JXS(8) + i - 1))

        ! Angular distribution given as part of a correlated angle-energy distribution
        if (LOCB == -1) cycle

        ! No angular distribution data are given for this reaction, isotropic
        ! scattering is assumed (in CM if TY < 0 and in LAB if TY > 0)
        if (LOCB == 0) cycle

        ! Loop over each separate energy distribution. Even though there is only
        ! "one" angular distribution, it is repeated as many times as there are
        ! energy distributions for this reaction since the
        ! UncorrelatedAngleEnergy type holds one angle and energy distribution.
        do k = 1, size(rxn % products(1) % distribution)
          select type (aedist => rxn % products(1) % distribution(k) % obj)
          type is (UncorrelatedAngleEnergy)
            ! allocate space for incoming energies and locations
            NE = int(XSS(JXS(9) + LOCB - 1))
            allocate(aedist % angle % energy(NE))
            allocate(aedist % angle % distribution(NE))
            allocate(LC(NE))

            ! read incoming energy grid and location of nucs
            XSS_index = JXS(9) + LOCB
            aedist % angle % energy(:) = get_real(NE)
            LC(:) = get_int(NE)

            ! determine dize of data block
            do j = 1, NE
              if (LC(j) == 0) then
                ! isotropic
                allocate(Uniform :: aedist % angle % distribution(j) % obj)
                select type (adist => aedist % angle % distribution(j) % obj)
                type is (Uniform)
                  adist % a = -ONE
                  adist % b = ONE
                end select

              elseif (LC(j) > 0) then
                ! 32 equiprobable bins
                allocate(Equiprobable :: aedist % angle % distribution(j) % obj)
                select type (adist => aedist % angle % distribution(j) % obj)
                type is (Equiprobable)
                  allocate(adist % x(33))
                end select

              elseif (LC(j) < 0) then
                ! tabular distribution
                allocate(Tabular :: aedist % angle % distribution(j) % obj)
              end if
            end do

            ! read angular distribution -- currently this does not actually parse the
            ! angular distribution tables for each incoming energy, that must be done
            ! on-the-fly
            do j = 1, NE
              XSS_index = JXS(9) + abs(LC(j)) - 1
              select type(adist => aedist % angle % distribution(j) % obj)
              type is (Equiprobable)
                adist % x(:) = get_real(33)
              type is (Tabular)
                ! determine interpolation and number of points
                interp = nint(XSS(XSS_index))
                NP = nint(XSS(XSS_index + 1))

                ! Get probability density data
                XSS_index = XSS_index + 2
                allocate(adist % x(NP), adist % p(NP), adist % c(NP))
                adist % x(:) = get_real(NP)
                adist % p(:) = get_real(NP)
                adist % c(:) = get_real(NP)
              end select
            end do
            deallocate(LC)

          end select
        end do
      end associate
    end do

  end subroutine read_angular_dist

!===============================================================================
! READ_ENERGY_DIST parses the secondary energy distribution for each reaction
! with seconary neutrons (except elastic scattering)
!===============================================================================

  subroutine read_energy_dist(nuc)
    type(Nuclide), intent(inout) :: nuc

    integer :: i     ! loop index
    integer :: n
    integer :: IDAT  ! locator for distribution data
    integer :: LNW   ! location of next energy law
    integer :: LAW   ! Type of energy law

    ! Loop over all reactions
    do i = 1, NXS(5)
      ! Determine how many energy distributions are present for this reaction
      LNW = nint(XSS(JXS(10) + i - 1))
      n = 0
      do while (LNW > 0)
        n = n + 1
        LNW = nint(XSS(JXS(11) + LNW - 1))
      end do

      ! Allocate space for distributions and probability of validity
      associate (p => nuc % reactions(i + 1) % products(1))
        allocate(p % applicability(n))
        allocate(p % distribution(n))

        LNW = nint(XSS(JXS(10) + i - 1))
        n = 0
        do while (LNW > 0)
          n = n + 1

          ! Determine energy law and location of data
          LAW  = nint(XSS(JXS(11) + LNW))
          IDAT = nint(XSS(JXS(11) + LNW + 1))

          ! Read probability of law validity
          call p % applicability(n) % from_ace(XSS, JXS(11) + LNW + 2)

          ! Read energy law data
          call get_energy_dist(p % distribution(n) % obj, LAW, &
               JXS(11), IDAT, nuc % awr, nuc % reactions(i + 1) % Q_value)

          ! <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<
          ! Before the secondary distribution refactor, when the angle/energy
          ! distribution was uncorrelated, no angle was actually sampled. With
          ! the refactor, an angle is always sampled for an uncorrelated
          ! distribution even when no angle distribution exists in the ACE file
          ! (isotropic is assumed). To preserve the RNG stream, we explicitly
          ! mark fission reactions so that we avoid the angle sampling.
          if (any(nuc % reactions(i + 1) % MT == &
               [N_FISSION, N_F, N_NF, N_2NF, N_3NF])) then
            select type (aedist => p % distribution(n) % obj)
            type is (UncorrelatedAngleEnergy)
              aedist % fission = .true.
            end select
          end if
          ! <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<

          ! Get locator for next distribution
          LNW = nint(XSS(JXS(11) + LNW - 1))
        end do
      end associate
    end do

  end subroutine read_energy_dist

!===============================================================================
! GET_ENERGY_DIST reads in data for a single law for an energy distribution and
! calls itself recursively if there are multiple energy distributions for a
! single reaction
!===============================================================================

  recursive subroutine get_energy_dist(aedist, law, LDIS, IDAT, awr, Q_value)
    class(AngleEnergy), allocatable, intent(inout) :: aedist
    integer, intent(in) :: law
    integer, intent(in) :: LDIS
    integer, intent(in) :: IDAT
    real(8), intent(in) :: awr
    real(8), intent(in) :: Q_value

    integer :: i, j
    integer :: NR     ! number of interpolation regions
    integer :: NE     ! number of incoming energies
    integer :: NP     ! number of outgoing energies/angles
    integer :: interp
    integer, allocatable :: L(:)  ! locations of distributions for each Ein
    integer, allocatable :: LC(:)  ! locations of distributions for each Ein

    XSS_index = LDIS + IDAT - 1

    if (law == 44) then
      allocate(KalbachMann :: aedist)
    elseif (law == 61) then
      allocate(CorrelatedAngleEnergy :: aedist)
    elseif (law == 66) then
      allocate(NBodyPhaseSpace :: aedist)
    else
      allocate(UncorrelatedAngleEnergy :: aedist)
    end if

    select type (aedist)
    type is (UncorrelatedAngleEnergy)
      ! ========================================================================
      ! UNCORRELATED ENERGY DISTRIBUTIONS

      select case (law)
      case (1)
        allocate(TabularEquiprobable :: aedist % energy)
        select type (edist => aedist % energy)
        type is (TabularEquiprobable)
          NR = nint(XSS(XSS_index))
          NE = nint(XSS(XSS_index + 1 + 2*NR))
          if (NR > 0) then
            call fatal_error("Multiple interpolation regions not yet supported &
                 &for tabular equiprobable energy distributions.")
          end if
          edist % n_region = NR

          ! Read incoming energies for which outgoing energies are tabulated
          allocate(edist % energy_in(NE))
          XSS_index = XSS_index + 2 + 2*NR
          edist % energy_in(:) = get_real(NE)

          ! Read outgoing energy tables
          NP = nint(XSS(XSS_index))
          allocate(edist % energy_out(NP, NE))
          XSS_index = XSS_index + 1
          do i = 1, NE
            edist % energy_out(:, i) = get_real(NP)
          end do
        end select

      case (3)
        allocate(LevelInelastic :: aedist % energy)
        select type (edist => aedist % energy)
        type is (LevelInelastic)
          edist % threshold = XSS(XSS_index)
          edist % mass_ratio = XSS(XSS_index + 1)
        end select

      case (4)
        allocate(ContinuousTabular :: aedist % energy)
        select type (edist => aedist % energy)
        type is (ContinuousTabular)
          NR = nint(XSS(XSS_index))
          XSS_index = XSS_index + 1
          if (NR > 1) then
            call fatal_error("Multiple interpolation regions not yet supported &
                 &for continuous tabular energy distributions.")
          end if
          edist % n_region = NR

          ! Read breakpoints and interpolation parameters
          if (NR > 0) then
            allocate(edist % breakpoints(NR))
            allocate(edist % interpolation(NR))
            edist % breakpoints(:) = get_int(NR)
            edist % interpolation(:) = get_int(NR)
          end if

          ! Read incoming energies for which outgoing energies are tabulated and
          ! locators
          NE = nint(XSS(XSS_index))
          XSS_index = XSS_index + 1
          allocate(edist % energy(NE))
          allocate(L(NE))
          edist % energy(:) = get_real(NE)
          L(:) = get_int(NE)

          ! Read outgoing energy tables
          allocate(edist % distribution(NE))
          do i = 1, NE
            ! Determine interpolation and number of discrete points
            XSS_index = LDIS + L(i) - 1
            interp = nint(XSS(XSS_index))
            edist % distribution(i) % interpolation = mod(interp, 10)
            edist % distribution(i) % n_discrete = (interp - &
                 edist % distribution(i) % interpolation)/10

            ! check for discrete lines present
            if (edist % distribution(i) % n_discrete > 0) then
              call fatal_error("Discrete lines in continuous tabular &
                   &distribution not yet supported")
            end if

            ! Determine number of points and allocate space
            NP = nint(XSS(XSS_index + 1))
            allocate(edist % distribution(i) % e_out(NP))
            allocate(edist % distribution(i) % p(NP))
            allocate(edist % distribution(i) % c(NP))

            ! Read tabular PDF for outgoing energy
            XSS_index = XSS_index + 2
            edist % distribution(i) % e_out(:) = get_real(NP)
            edist % distribution(i) % p(:) = get_real(NP)
            edist % distribution(i) % c(:) = get_real(NP)
          end do

          deallocate(L)
        end select

      case (7)
        allocate(MaxwellEnergy :: aedist % energy)
        select type (edist => aedist % energy)
        type is (MaxwellEnergy)
          call edist % theta % from_ace(XSS, XSS_index)
          edist % u = XSS(XSS_index + 2 + 2*edist % theta % n_regions + &
               2*edist % theta % n_pairs)
        end select

      case (9)
        allocate(Evaporation :: aedist % energy)
        select type(edist => aedist % energy)
        type is (Evaporation)
          call edist % theta % from_ace(XSS, XSS_index)
          edist % u = XSS(XSS_index + 2 + 2*edist % theta % n_regions + &
               2*edist % theta % n_pairs)
        end select

      case (11)
        allocate(WattEnergy :: aedist % energy)
        select type(edist => aedist % energy)
        type is (WattEnergy)
          call edist % a % from_ace(XSS, XSS_index)
          XSS_index = XSS_index + 2 + 2*edist % a % n_regions + 2*edist % a % n_pairs
          call edist % b % from_ace(XSS, XSS_index)
          XSS_index = XSS_index + 2 + 2*edist % b % n_regions + 2*edist % b % n_pairs
          edist % u = XSS(XSS_index)
        end select

      end select

    type is (KalbachMann)
      ! ========================================================================
      ! CORRELATED KALBACH-MANN DISTRIBUTION

      NR = int(XSS(XSS_index))
      NE = int(XSS(XSS_index + 1 + 2*NR))
      if (NR > 0) then
        call fatal_error("Multiple interpolation regions not yet supported &
             &for Kalbach-Mann energy distributions.")
      end if
      aedist % n_region = NR

      ! Read incoming energies for which outgoing energies are tabulated and locators
      allocate(aedist % energy(NE))
      allocate(L(NE))
      XSS_index = XSS_index + 2 + 2*NR
      aedist % energy(:) = get_real(NE)
      L(:) = get_int(NE)

      ! Read outgoing energy tables
      allocate(aedist % distribution(NE))
      do i = 1, NE
        ! Determine interpolation and number of discrete points
        XSS_index = LDIS + L(i) - 1
        interp = nint(XSS(XSS_index))
        aedist % distribution(i) % interpolation = mod(interp, 10)
        aedist % distribution(i) % n_discrete = (interp - aedist % distribution(i) % interpolation)/10

        ! check for discrete lines present
        if (aedist % distribution(i) % n_discrete > 0) then
          call fatal_error("Discrete lines in Kalbach-Mann distribution not &
               &yet supported")
        end if

        ! Determine number of points and allocate space
        NP = nint(XSS(XSS_index + 1))
        allocate(aedist % distribution(i) % e_out(NP))
        allocate(aedist % distribution(i) % p(NP))
        allocate(aedist % distribution(i) % c(NP))
        allocate(aedist % distribution(i) % r(NP))
        allocate(aedist % distribution(i) % a(NP))

        ! Read tabular PDF for outgoing energy
        XSS_index = XSS_index + 2
        aedist % distribution(i) % e_out(:) = get_real(NP)
        aedist % distribution(i) % p(:) = get_real(NP)
        aedist % distribution(i) % c(:) = get_real(NP)
        aedist % distribution(i) % r(:) = get_real(NP)
        aedist % distribution(i) % a(:) = get_real(NP)
      end do

      deallocate(L)

    type is (CorrelatedAngleEnergy)
      ! ========================================================================
      ! CORRELATED ANGLE-ENERGY DISTRIBUTION

      NR = int(XSS(XSS_index))
      NE = int(XSS(XSS_index + 1 + 2*NR))
      if (NR > 0) then
        call fatal_error("Multiple interpolation regions not yet supported &
             &for correlated angle-energy distributions.")
      end if
      aedist % n_region = NR

      ! Read incoming energies for which outgoing energies are tabulated and
      ! locators
      allocate(aedist % energy(NE))
      allocate(L(NE))
      XSS_index = XSS_index + 2 + 2*NR
      aedist % energy(:) = get_real(NE)
      L(:) = get_int(NE)

      ! Read outgoing energy tables
      allocate(aedist % distribution(NE))
      do i = 1, NE
        ! Determine interpolation and number of discrete points
        XSS_index = LDIS + L(i) - 1
        interp = nint(XSS(XSS_index))
        aedist % distribution(i) % interpolation = mod(interp, 10)
        aedist % distribution(i) % n_discrete = (interp - aedist % distribution(i) % interpolation)/10

        ! check for discrete lines present
        if (aedist % distribution(i) % n_discrete > 0) then
          call fatal_error("Discrete lines in correlated angle-energy &
               &distribution not yet supported")
        end if

        ! Determine number of points and allocate space
        NP = nint(XSS(XSS_index + 1))
        allocate(aedist % distribution(i) % e_out(NP))
        allocate(aedist % distribution(i) % p(NP))
        allocate(aedist % distribution(i) % c(NP))
        allocate(LC(NP))

        ! Read tabular PDF for outgoing energy
        XSS_index = XSS_index + 2
        aedist % distribution(i) % e_out(:) = get_real(NP)
        aedist % distribution(i) % p(:) = get_real(NP)
        aedist % distribution(i) % c(:) = get_real(NP)
        LC(:) = get_int(NP)

        ! allocate angular distributions for each incoming/outgoing energy
        allocate(aedist % distribution(i) % angle(NP))
        do j = 1, NP
          if (LC(j) == 0) then
            ! isotropic
            allocate(Uniform :: aedist % distribution(i) % angle(j) % obj)
            select type (adist => aedist % distribution(i) % angle(j) % obj)
            type is (Uniform)
              adist % a = -ONE
              adist % b = ONE
            end select

          elseif (LC(j) > 0) then
            ! tabular distribution
            allocate(Tabular :: aedist % distribution(i) % angle(j) % obj)
          end if
        end do

        ! read angular distributions
        do j = 1, NP
          XSS_index = LDIS + abs(LC(j)) - 1
          select type(adist => aedist % distribution(i) % angle(j) % obj)
          type is (Tabular)
            ! determine interpolation and number of points
            interp = nint(XSS(XSS_index))
            NP = nint(XSS(XSS_index + 1))

            ! Get probability density data
            XSS_index = XSS_index + 2
            allocate(adist % x(NP), adist % p(NP), adist % c(NP))
            adist % x(:) = get_real(NP)
            adist % p(:) = get_real(NP)
            adist % c(:) = get_real(NP)
          end select
        end do
        deallocate(LC)

      end do

      deallocate(L)

    type is (NBodyPhaseSpace)
      ! ========================================================================
      ! N-BODY PHASE SPACE DISTRIBUTION

      aedist % n_bodies = int(XSS(XSS_index))
      aedist % mass_ratio = XSS(XSS_index + 1)
      aedist % A = awr
      aedist % Q = Q_value
    end select

  end subroutine get_energy_dist

!===============================================================================
! READ_UNR_RES reads in unresolved resonance probability tables if present.
!===============================================================================

  subroutine read_unr_res(nuc)
    type(Nuclide), intent(inout) :: nuc

    integer :: JXS23 ! location of URR data
    integer :: lc    ! locator
    integer :: N     ! # of incident energies
    integer :: M     ! # of probabilities
    integer :: i     ! index over incoming energies
    integer :: j     ! index over values
    integer :: k     ! index over probabilities

    ! determine locator for URR data
    JXS23 = JXS(23)

    ! check if URR data is present
    if (JXS23 /= 0) then
      nuc % urr_present = .true.
      allocate(nuc % urr_data)
      lc = JXS23
    else
      nuc % urr_present = .false.
      return
    end if

    ! read parameters
    nuc % urr_data % n_energy = int(XSS(lc))
    nuc % urr_data % n_prob = int(XSS(lc + 1))
    nuc % urr_data % interp = int(XSS(lc + 2))
    nuc % urr_data % inelastic_flag = int(XSS(lc + 3))
    nuc % urr_data % absorption_flag = int(XSS(lc + 4))
    if (int(XSS(lc + 5)) == 0) then
      nuc % urr_data % multiply_smooth = .false.
    else
      nuc % urr_data % multiply_smooth = .true.
    end if

    ! if the inelastic competition flag indicates that the inelastic cross
    ! section should be determined from a normal reaction cross section, we need
    ! to set up a pointer to that reaction
    nuc % urr_inelastic = NONE
    if (nuc % urr_data % inelastic_flag > 0) then
      do i = 1, nuc % n_reaction
        if (nuc % reactions(i) % MT == nuc % urr_data % inelastic_flag) then
          nuc % urr_inelastic = i
        end if
      end do

      ! Abort if no corresponding inelastic reaction was found
      if (nuc % urr_inelastic == NONE) then
        call fatal_error("Could not find inelastic reaction specified on &
             &unresolved resonance probability table.")
      end if
    end if

    ! allocate incident energies and probability tables
    N = nuc % urr_data % n_energy
    M = nuc % urr_data % n_prob
    allocate(nuc % urr_data % energy(N))
    allocate(nuc % urr_data % prob(N,6,M))

    ! read incident energies
    XSS_index = lc + 6
    nuc % urr_data % energy = get_real(N)

    ! read probability tables
    do i = 1, N
      do j = 1, 6
        do k = 1, M
          nuc % urr_data % prob(i,j,k) = XSS(XSS_index)
          XSS_index = XSS_index + 1
        end do
      end do
    end do

    ! Check for negative values
    if (any(nuc % urr_data % prob < ZERO)) then
      if (master) call warning("Negative value(s) found on probability table &
           &for nuclide " // nuc % name)
    end if

  end subroutine read_unr_res
!===============================================================================
! GENERATE_NU_FISSION precalculates the microscopic nu-fission cross section for
! a given nuclide. This is done so that the nu_total function does not need to
! be called during cross section lookups.
!===============================================================================

  subroutine generate_nu_fission(nuc)
    type(Nuclide), intent(inout) :: nuc

    integer :: i  ! index on nuclide energy grid

    do i = 1, size(nuc % energy)
      nuc % nu_fission(i) = nuc % nu(nuc % energy(i), EMISSION_TOTAL) * &
           nuc % fission(i)
    end do
  end subroutine generate_nu_fission

!===============================================================================
! READ_THERMAL_DATA reads elastic and inelastic cross sections and corresponding
! secondary energy/angle distributions derived from experimental S(a,b)
! data. Namely, this routine reads the ITIE, ITCE, ITXE, and ITCA blocks.
!===============================================================================

  subroutine read_thermal_data(table)
    type(SAlphaBeta), intent(inout) :: table

    integer :: i      ! index for incoming energies
    integer :: j      ! index for outgoing energies
    integer :: k      ! index for outoging angles
    integer :: lc     ! location in XSS array
    integer :: NE_in  ! number of incoming energies
    integer :: NE_out ! number of outgoing energies
    integer :: NMU    ! number of outgoing angles
    integer :: JXS4   ! location of elastic energy table
    integer(8), allocatable :: LOCC(:) ! Location of inelastic data

    ! read secondary energy mode for inelastic scattering
    table % secondary_mode = NXS(7)

    ! read number of inelastic energies and allocate arrays
    NE_in = int(XSS(JXS(1)))
    table % n_inelastic_e_in = NE_in
    allocate(table % inelastic_e_in(NE_in))
    allocate(table % inelastic_sigma(NE_in))

    ! read inelastic energies and cross-sections
    XSS_index = JXS(1) + 1
    table % inelastic_e_in = get_real(NE_in)
    table % inelastic_sigma = get_real(NE_in)

    ! set threshold value
    table % threshold_inelastic = table % inelastic_e_in(NE_in)

    ! allocate space for outgoing energy/angle for inelastic
    ! scattering
    if (table % secondary_mode == SAB_SECONDARY_EQUAL .or. &
         table % secondary_mode == SAB_SECONDARY_SKEWED) then
      NMU = NXS(3) + 1
      table % n_inelastic_mu = NMU
      NE_out = NXS(4)
      table % n_inelastic_e_out = NE_out
      allocate(table % inelastic_e_out(NE_out, NE_in))
      allocate(table % inelastic_mu(NMU, NE_out, NE_in))
    else if (table % secondary_mode == SAB_SECONDARY_CONT) then
      NMU = NXS(3) - 1
      table % n_inelastic_mu = NMU
      allocate(table % inelastic_data(NE_in))
      allocate(LOCC(NE_in))
      ! NE_out will be determined later
    end if

    ! read outgoing energy/angle distribution for inelastic scattering
    if (table % secondary_mode == SAB_SECONDARY_EQUAL .or. &
         table % secondary_mode == SAB_SECONDARY_SKEWED) then
      lc = JXS(3) - 1
      do i = 1, NE_in
        do j = 1, NE_out
          ! read outgoing energy
          table % inelastic_e_out(j,i) = XSS(lc + 1)

          ! read outgoing angles for this outgoing energy
          do k = 1, NMU
            table % inelastic_mu(k,j,i) = XSS(lc + 1 + k)
          end do

          ! advance pointer
          lc = lc + 1 + NMU
        end do
      end do
    else if (table % secondary_mode == SAB_SECONDARY_CONT) then
      ! Get the location pointers to each Ein's DistEnergySAB data
      LOCC = get_int(NE_in)
      ! Get the number of outgoing energies and allocate space accordingly
      do i = 1, NE_in
        NE_out = int(XSS(XSS_index + i - 1))
        table % inelastic_data(i) % n_e_out = NE_out
        allocate(table % inelastic_data(i) % e_out (NE_out))
        allocate(table % inelastic_data(i) % e_out_pdf (NE_out))
        allocate(table % inelastic_data(i) % e_out_cdf (NE_out))
        allocate(table % inelastic_data(i) % mu (NMU, NE_out))
      end do

      ! Now we can fill the inelastic_data(i) attributes
      do i = 1, NE_in
        XSS_index = int(LOCC(i))
        NE_out = table % inelastic_data(i) % n_e_out
        do j = 1, NE_out
          table % inelastic_data(i) % e_out(j) = XSS(XSS_index + 1)
          table % inelastic_data(i) % e_out_pdf(j) = XSS(XSS_index + 2)
          table % inelastic_data(i) % e_out_cdf(j) = XSS(XSS_index + 3)
          table % inelastic_data(i) % mu(:, j) = &
               XSS(XSS_index + 4: XSS_index + 4 + NMU - 1)
          XSS_index = XSS_index + 4 + NMU - 1
        end do
      end do
    end if

    ! read number of elastic energies and allocate arrays
    JXS4 = JXS(4)
    if (JXS4 /= 0) then
      NE_in = int(XSS(JXS4))
      table % n_elastic_e_in = NE_in
      allocate(table % elastic_e_in(NE_in))
      allocate(table % elastic_P(NE_in))

      ! read elastic energies and P
      XSS_index = JXS4 + 1
      table % elastic_e_in = get_real(NE_in)
      table % elastic_P    = get_real(NE_in)

      ! set threshold
      table % threshold_elastic = table % elastic_e_in(NE_in)

      ! determine whether sigma=P or sigma = P/E
      table % elastic_mode = NXS(5)
    else
      table % threshold_elastic = ZERO
      table % n_elastic_e_in = 0
    end if

    ! allocate space for outgoing energy/angle for elastic scattering
    NMU = NXS(6) + 1
    table % n_elastic_mu = NMU
    if (NMU > 0) then
      allocate(table % elastic_mu(NMU, NE_in))
    end if

    ! read equiprobable outgoing cosines for elastic scattering each
    ! incoming energy
    if (JXS4 /= 0 .and. NMU /= 0) then
      lc = JXS(6) - 1
      do i = 1, NE_in
        do j = 1, NMU
          table % elastic_mu(j,i) = XSS(lc + j)
        end do
        lc = lc + NMU
      end do
    end if

  end subroutine read_thermal_data

!===============================================================================
! GET_INT returns an array of integers read from the current position in the XSS
! array
!===============================================================================

  function get_int(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    integer             :: array(n_values) ! array of values

    array = int(XSS(XSS_index:XSS_index + n_values - 1))
    XSS_index = XSS_index + n_values

  end function get_int

!===============================================================================
! GET_REAL returns an array of real(8)s read from the current position in the
! XSS array
!===============================================================================

  function get_real(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    real(8)             :: array(n_values) ! array of values

    array = XSS(XSS_index:XSS_index + n_values - 1)
    XSS_index = XSS_index + n_values

  end function get_real

end module ace
