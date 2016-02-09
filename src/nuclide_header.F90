module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV

  use ace_header
  use constants
  use endf,        only: reaction_name
  use list_header, only: ListInt
  use math,        only: evaluate_legendre
  use string

  implicit none

!===============================================================================
! Nuclide contains the base nuclidic data for a nuclide, which does not depend
! upon how the nuclear data is represented (i.e., CE, or any variant of MG).
! The extended types, NuclideCE and NuclideMG deal with the rest
!===============================================================================

  type, abstract :: Nuclide
    character(12) :: name    ! name of nuclide, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    real(8)       :: awr     ! Atomic Weight Ratio
    integer       :: listing ! index in xs_listings
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Linked list of indices in nuclides array of instances of this same nuclide
    type(VectorInt) :: nuc_list

    ! Fission information
    logical :: fissionable         ! nuclide is fissionable?

  contains
    procedure(print_nuclide_),     deferred :: print ! Writes nuclide info
  end type Nuclide

  abstract interface

    subroutine print_nuclide_(this, unit)
      import Nuclide
      class(Nuclide),intent(in)     :: this
      integer, optional, intent(in) :: unit
    end subroutine print_nuclide_

  end interface

  type, extends(Nuclide) :: NuclideCE
    ! Energy grid information
    integer :: n_grid                     ! # of nuclide grid points
    integer, allocatable :: grid_index(:) ! log grid mapping indices
    real(8), allocatable :: energy(:)     ! energy values corresponding to xs

    ! Microscopic cross sections
    real(8), allocatable :: total(:)      ! total cross section
    real(8), allocatable :: elastic(:)    ! elastic scattering
    real(8), allocatable :: fission(:)    ! fission
    real(8), allocatable :: nu_fission(:) ! neutron production
    real(8), allocatable :: absorption(:) ! absorption (MT > 100)
    real(8), allocatable :: heating(:)    ! heating

    ! Resonance scattering info
    logical              :: resonant = .false. ! resonant scatterer?
    character(10)        :: name_0K = '' ! name of 0K nuclide, e.g. 92235.00c
    character(16)        :: scheme ! target velocity sampling scheme
    integer              :: n_grid_0K ! number of 0K energy grid points
    real(8), allocatable :: energy_0K(:)  ! energy grid for 0K xs
    real(8), allocatable :: elastic_0K(:) ! Microscopic elastic cross section
    real(8), allocatable :: xs_cdf(:) ! CDF of v_rel times cross section
    real(8)              :: E_min ! lower cutoff energy for res scattering
    real(8)              :: E_max ! upper cutoff energy for res scattering

    ! Fission information
    logical :: has_partial_fission ! nuclide has partial fission reactions?
    integer :: n_fission           ! # of fission reactions
    integer, allocatable :: index_fission(:) ! indices in reactions

    ! Total fission neutron emission
    integer :: nu_t_type
    real(8), allocatable :: nu_t_data(:)

    ! Prompt fission neutron emission
    integer :: nu_p_type
    real(8), allocatable :: nu_p_data(:)

    ! Delayed fission neutron emission
    integer :: nu_d_type
    integer :: n_precursor ! # of delayed neutron precursors
    real(8), allocatable :: nu_d_data(:)
    real(8), allocatable :: nu_d_precursor_data(:)
    type(AngleEnergyContainer), allocatable :: nu_d_edist(:)

    ! Unresolved resonance data
    logical                :: urr_present
    integer                :: urr_inelastic
    type(UrrData), pointer :: urr_data => null()

    ! Reactions
    integer :: n_reaction ! # of reactions
    type(Reaction), allocatable :: reactions(:)
    type(DictIntInt) :: reaction_index ! map MT values to index in reactions
                                       ! array; used at tally-time

  contains
    procedure :: clear => nuclidece_clear
    procedure :: print => nuclidece_print
  end type NuclideCE

  type, abstract, extends(Nuclide) :: NuclideMG
    ! Scattering Order Information
    integer :: order      ! Order of data (Scattering for NuclideIso,
                          ! Number of angles for all in NuclideAngle)
    integer :: scatt_type ! either legendre, histogram, or tabular.
    integer :: legendre_mu_points ! Number of tabular points to use to represent
                                  ! Legendre distribs, -1 if sample with the
                                  ! Legendres themselves
  contains
    procedure(nuclidemg_get_xs), deferred :: get_xs ! Get the xs
    procedure(nuclide_calc_f_), deferred  :: calc_f ! Calculates f, given mu
  end type NuclideMG

  abstract interface
    function nuclidemg_get_xs(this, g, xstype, gout, uvw, mu, i_azi, i_pol) &
         result(xs)
      import NuclideMG
      class(NuclideMG), intent(in) :: this
      integer, intent(in)           :: g      ! Incoming Energy group
      character(*), intent(in)      :: xstype ! Cross Section Type
      integer, optional, intent(in) :: gout   ! Outgoing Group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      integer, optional, intent(in) :: i_azi  ! Azimuthal Index
      integer, optional, intent(in) :: i_pol  ! Polar Index
      real(8)                       :: xs     ! Resultant xs
    end function nuclidemg_get_xs

    pure function nuclide_calc_f_(this, gin, gout, mu, uvw, i_azi, i_pol) result(f)
      import NuclideMG
      class(NuclideMG), intent(in) :: this
      integer, intent(in)           :: gin   ! Incoming Energy Group
      integer, intent(in)           :: gout  ! Outgoing Energy Group
      real(8), intent(in)           :: mu    ! Angle of interest
      real(8), intent(in), optional :: uvw(3) ! Direction vector
      integer, intent(in), optional :: i_azi ! Incoming Energy Group
      integer, intent(in), optional :: i_pol ! Outgoing Energy Group
      real(8)                       :: f     ! Return value of f(mu)

    end function nuclide_calc_f_
  end interface

!===============================================================================
! NuclideIso contains the base MGXS data for a nuclide specifically for
! isotropically weighted MGXS
!===============================================================================

  type, extends(NuclideMG) :: NuclideIso

    ! Microscopic cross sections
    real(8), allocatable :: total(:)        ! total cross section
    real(8), allocatable :: absorption(:)   ! absorption cross section
    real(8), allocatable :: scatter(:,:,:)  ! scattering information
    real(8), allocatable :: nu_fission(:,:) ! fission matrix (Gout x Gin)
    real(8), allocatable :: k_fission(:)    ! kappa-fission
    real(8), allocatable :: fission(:)      ! neutron production
    real(8), allocatable :: chi(:)          ! Fission Spectra
    real(8), allocatable :: mult(:,:)       ! Scatter multiplicity (Gout x Gin)

  contains
    procedure :: print  => nuclideiso_print  ! Writes nuclide info
    procedure :: get_xs => nuclideiso_get_xs ! Gets Size of Data w/in Object
    procedure :: calc_f => nuclideiso_calc_f ! Calcs f given mu
  end type NuclideIso

!===============================================================================
! NuclideAngle contains the base MGXS data for a nuclide specifically for
! explicit angle-dependent weighted MGXS
!===============================================================================

  type, extends(NuclideMG) :: NuclideAngle

    ! Microscopic cross sections. Dimensions are: (Npol, Nazi, Nl, Ng, Ng)
    real(8), allocatable :: total(:,:,:)        ! total cross section
    real(8), allocatable :: absorption(:,:,:)   ! absorption cross section
    real(8), allocatable :: scatter(:,:,:,:,:)  ! scattering information
    real(8), allocatable :: nu_fission(:,:,:,:) ! fission matrix (Gout x Gin)
    real(8), allocatable :: k_fission(:,:,:)    ! kappa-fission
    real(8), allocatable :: fission(:,:,:)      ! neutron production
    real(8), allocatable :: chi(:,:,:)          ! Fission Spectra
    real(8), allocatable :: mult(:,:,:,:)       ! Scatter multiplicity (Gout x Gin)

    ! In all cases, right-most indices are theta, phi
    integer              :: Npol         ! Number of polar angles
    integer              :: Nazi         ! Number of azimuthal angles
    real(8), allocatable :: polar(:)     ! polar angles
    real(8), allocatable :: azimuthal(:) ! azimuthal angles

  contains
    procedure :: print  => nuclideangle_print  ! Gets Size of Data w/in Object
    procedure :: get_xs => nuclideangle_get_xs ! Gets Size of Data w/in Object
    procedure :: calc_f => nuclideangle_calc_f ! Calcs f given mu
  end type NuclideAngle

!===============================================================================
! NUCLIDEMGCONTAINER pointer array for storing Nuclides
!===============================================================================

  type NuclideMGContainer
    class(NuclideMG), pointer :: obj
  end type NuclideMGContainer

!===============================================================================
! NUCLIDE0K temporarily contains all 0K cross section data and other parameters
! needed to treat resonance scattering before transferring them to NuclideCE
!===============================================================================

  type Nuclide0K
    character(10) :: nuclide             ! name of nuclide, e.g. U-238
    character(16) :: scheme = 'ares'     ! target velocity sampling scheme
    character(10) :: name                ! name of nuclide, e.g. 92235.03c
    character(10) :: name_0K             ! name of 0K nuclide, e.g. 92235.00c
    real(8)       :: E_min = 0.01e-6_8   ! lower cutoff energy for res scattering
    real(8)       :: E_max = 1000.0e-6_8 ! upper cutoff energy for res scattering
  end type Nuclide0K

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a
! particular nuclide at the current energy
!===============================================================================

  type NuclideMicroXS
    integer :: index_grid      ! index on nuclide energy grid
    integer :: index_temp      ! temperature index for nuclide
    real(8) :: last_E = ZERO   ! last evaluated energy
    real(8) :: interp_factor   ! interpolation factor on nuc. energy grid
    real(8) :: total           ! microscropic total xs
    real(8) :: elastic         ! microscopic elastic scattering xs
    real(8) :: absorption      ! microscopic absorption xs
    real(8) :: fission         ! microscopic fission xs
    real(8) :: nu_fission      ! microscopic production xs

    ! Information for S(a,b) use
    integer :: index_sab          ! index in sab_tables (zero means no table)
    integer :: last_index_sab = 0 ! index in sab_tables last used by this nuclide
    real(8) :: elastic_sab        ! microscopic elastic scattering on S(a,b) table

    ! Information for URR probability table use
    logical :: use_ptable  ! in URR range with probability tables?
    real(8) :: last_prn
  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type MaterialMacroXS
    real(8) :: total         ! macroscopic total xs
    real(8) :: elastic       ! macroscopic elastic scattering xs
    real(8) :: absorption    ! macroscopic absorption xs
    real(8) :: fission       ! macroscopic fission xs
    real(8) :: nu_fission    ! macroscopic production xs
  end type MaterialMacroXS

!===============================================================================
! XSLISTING contains data read from a CE or MG cross_sections.xml file
! (or equivalent)
!===============================================================================

  type XsListing
    character(12) :: name       ! table name, e.g. 92235.70c
    character(12) :: alias      ! table alias, e.g. U-235.70c
    integer       :: type       ! type of table (cont-E neutron, S(A,b), etc)
    integer       :: zaid       ! ZAID identifier = 1000*Z + A
    integer       :: filetype   ! ASCII or BINARY
    integer       :: location   ! location of table within library
    integer       :: recl       ! record length for library
    integer       :: entries    ! number of entries per record
    real(8)       :: awr        ! atomic weight ratio (# of neutron masses)
    real(8)       :: kT         ! Boltzmann constant * temperature (MeV)
    logical       :: metastable ! is this nuclide metastable?
    character(MAX_FILE_LEN) :: path ! path to library containing table
  end type XsListing

  contains

!===============================================================================
! NUCLIDECE_CLEAR resets and deallocates data in Nuclide, NuclideIso
! or NuclideAngle
!===============================================================================

    subroutine nuclidece_clear(this)

      class(NuclideCE), intent(inout) :: this ! The Nuclide object to clear

      integer :: i ! Loop counter

      if (associated(this % urr_data)) deallocate(this % urr_data)

      if (allocated(this % reactions)) then
        do i = 1, size(this % reactions)
          call this % reactions(i) % clear()
        end do
      end if

      call this % reaction_index % clear()

    end subroutine nuclidece_clear

!===============================================================================
! NUCLIDE*_PRINT displays information about a continuous-energy neutron
! cross_section table and its reactions and secondary angle/energy distributions
!===============================================================================

    subroutine nuclidece_print(this, unit)
      class(NuclideCE), intent(in) :: this
      integer, intent(in), optional :: unit

      integer :: i                 ! loop index over nuclides
      integer :: unit_             ! unit to write to
      integer :: size_xs           ! memory used for cross-sections (bytes)
      integer :: size_urr          ! memory used for probability tables (bytes)
      type(UrrData),  pointer :: urr

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Initialize totals
      size_urr = 0
      size_xs = 0

      ! Basic nuclide information
      write(unit_,*) 'Nuclide ' // trim(this % name)
      write(unit_,*) '  zaid = ' // trim(to_str(this % zaid))
      write(unit_,*) '  awr = ' // trim(to_str(this % awr))
      write(unit_,*) '  kT = ' // trim(to_str(this % kT))
      write(unit_,*) '  # of grid points = ' // trim(to_str(this % n_grid))
      write(unit_,*) '  Fissionable = ', this % fissionable
      write(unit_,*) '  # of fission reactions = ' // trim(to_str(this % n_fission))
      write(unit_,*) '  # of reactions = ' // trim(to_str(this % n_reaction))

      ! Information on each reaction
      write(unit_,*) '  Reaction     Q-value  COM    IE'
      do i = 1, this % n_reaction
        associate (rxn => this % reactions(i))
          write(unit_,'(3X,A11,1X,F8.3,3X,L1,3X,I6)') &
               reaction_name(rxn % MT), rxn % Q_value, rxn % scatter_in_cm, &
               rxn % threshold

          ! Accumulate data size
          size_xs = size_xs + (this % n_grid - rxn%threshold + 1) * 8
        end associate
      end do

      ! Add memory required for summary reactions (total, absorption, fission,
      ! nu-fission)
      size_xs = 8 * this % n_grid * 4

      ! Write information about URR probability tables
      size_urr = 0
      if (this % urr_present) then
        urr => this % urr_data
        write(unit_,*) '  Unresolved resonance probability table:'
        write(unit_,*) '    # of energies = ' // trim(to_str(urr % n_energy))
        write(unit_,*) '    # of probabilities = ' // trim(to_str(urr % n_prob))
        write(unit_,*) '    Interpolation =  ' // trim(to_str(urr % interp))
        write(unit_,*) '    Inelastic flag = ' // trim(to_str(urr % inelastic_flag))
        write(unit_,*) '    Absorption flag = ' // trim(to_str(urr % absorption_flag))
        write(unit_,*) '    Multiply by smooth? ', urr % multiply_smooth
        write(unit_,*) '    Min energy = ', trim(to_str(urr % energy(1)))
        write(unit_,*) '    Max energy = ', trim(to_str(urr % energy(urr % n_energy)))

        ! Calculate memory used by probability tables and add to total
        size_urr = urr % n_energy * (urr % n_prob * 6 + 1) * 8
      end if

      ! Write memory used
      write(unit_,*) '  Memory Requirements'
      write(unit_,*) '    Cross sections = ' // trim(to_str(size_xs)) // ' bytes'
      write(unit_,*) '    Probability Tables = ' // &
           trim(to_str(size_urr)) // ' bytes'

      ! Blank line at end of nuclide
      write(unit_,*)
    end subroutine nuclidece_print

    subroutine nuclidemg_print(this, unit_)
      class(NuclideMG), intent(in) :: this
      integer, intent(in)           :: unit_

      character(MAX_LINE_LEN) :: temp_str

      ! Basic nuclide information
      write(unit_,*) 'Nuclide ' // trim(this % name)
      if (this % zaid > 0) then
        ! Dont print if data was macroscopic and thus zaid & AWR would be nonsense
        write(unit_,*) '  zaid = ' // trim(to_str(this % zaid))
        write(unit_,*) '  awr = ' // trim(to_str(this % awr))
      end if
      write(unit_,*) '  kT = ' // trim(to_str(this % kT))
      if (this % scatt_type == ANGLE_LEGENDRE) then
        temp_str = "Legendre"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        write(unit_,*) '  # of Scatter Moments = ' // &
             trim(to_str(this % order))
      else if (this % scatt_type == ANGLE_HISTOGRAM) then
        temp_str = "Histogram"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        write(unit_,*) '  # of Scatter Bins = ' // &
             trim(to_str(this % order))
      else if (this % scatt_type == ANGLE_TABULAR) then
        temp_str = "Tabular"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        write(unit_,*) '  # of Scatter Points = ' // trim(to_str(this % order))
      end if
      write(unit_,*) '  Fissionable = ', this % fissionable

    end subroutine nuclidemg_print

    subroutine nuclideiso_print(this, unit)

      class(NuclideIso), intent(in) :: this
      integer, optional, intent(in)  :: unit

      integer :: unit_             ! unit to write to
      integer :: size_total, size_scattmat, size_mgxs

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Write Basic Nuclide Information
      call nuclidemg_print(this, unit_)

      ! Determine size of mgxs and scattering matrices
      size_scattmat = (size(this % scatter) + size(this % mult)) * 8
      size_mgxs = size(this % total) + size(this % absorption) + &
           size(this % nu_fission) + size(this % k_fission) + &
           size(this % fission) + size(this % chi)
      size_mgxs = size_mgxs * 8

      ! Calculate total memory
      size_total = size_scattmat + size_mgxs

      ! Write memory used
      write(unit_,*) '  Memory Requirements'
      write(unit_,*) '    Cross sections = ' // trim(to_str(size_mgxs)) // ' bytes'
      write(unit_,*) '    Scattering Matrices = ' // &
           trim(to_str(size_scattmat)) // ' bytes'
      write(unit_,*) '    Total = ' // trim(to_str(size_total)) // ' bytes'

      ! Blank line at end of nuclide
      write(unit_,*)

    end subroutine nuclideiso_print

    subroutine nuclideangle_print(this, unit)

      class(NuclideAngle), intent(in) :: this
      integer, optional, intent(in)    :: unit

      integer :: unit_             ! unit to write to
      integer :: size_total, size_scattmat, size_mgxs

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Write Basic Nuclide Information
      call nuclidemg_print(this, unit_)
      write(unit_,*) '  # of Polar Angles = ' // trim(to_str(this % Npol))
      write(unit_,*) '  # of Azimuthal Angles = ' // trim(to_str(this % Nazi))

      ! Determine size of mgxs and scattering matrices
      size_scattmat = (size(this % scatter) + size(this % mult)) * 8
      size_mgxs = size(this % total) + size(this % absorption) + &
           size(this % nu_fission) + size(this % k_fission) + &
           size(this % fission) + size(this % chi)
      size_mgxs = size_mgxs * 8

      ! Calculate total memory
      size_total = size_scattmat + size_mgxs

      ! Write memory used
      write(unit_,*) '  Memory Requirements'
      write(unit_,*) '    Cross sections = ' // trim(to_str(size_mgxs)) // ' bytes'
      write(unit_,*) '    Scattering Matrices = ' // &
           trim(to_str(size_scattmat)) // ' bytes'
      write(unit_,*) '    Total = ' // trim(to_str(size_total)) // ' bytes'

      ! Blank line at end of nuclide
      write(unit_,*)


    end subroutine nuclideangle_print

!===============================================================================
! NUCLIDE*_GET_XS Returns the requested data type
!===============================================================================

    function nuclideiso_get_xs(this, g, xstype, gout, uvw, mu, i_azi, i_pol) &
         result(xs)
      class(NuclideIso), intent(in) :: this
      integer, intent(in)            :: g      ! Incoming Energy group
      character(*), intent(in)       :: xstype ! Cross Section Type
      integer, optional, intent(in)  :: gout   ! Outgoing Group
      real(8), optional, intent(in)  :: uvw(3) ! Requested Angle
      real(8), optional, intent(in)  :: mu     ! Change in angle
      integer, optional, intent(in)  :: i_azi  ! Azimuthal Index
      integer, optional, intent(in)  :: i_pol  ! Polar Index
      real(8)                        :: xs     ! Resultant xs

      xs = ZERO

      if ((xstype == 'nu_fission' .or. xstype == 'fission' .or. xstype =='chi' &
           .or. xstype =='k_fission') .and. (.not. this % fissionable)) then
        return
      end if

      if (present(gout)) then
        select case(xstype)
        case('mult')
          xs = this % mult(gout,g)
        case('nu_fission')
          xs = this % nu_fission(gout,g)
        case('f_mu', 'f_mu/mult')
          xs = this % calc_f(g, gout, mu)
          if (xstype == 'f_mu/mult') then
            xs = xs / this % mult(gout,g)
          end if
        end select
      else
        select case(xstype)
        case('total')
          xs = this % total(g)
        case('absorption')
          xs = this % absorption(g)
        case('fission')
          xs = this % fission(g)
        case('k_fission')
          if (allocated(this % k_fission)) then
            xs = this % k_fission(g)
          end if
        case('chi')
          xs = this % chi(g)
        case('scatter')
          xs = this % total(g) - this % absorption(g)
        end select
      end if
    end function nuclideiso_get_xs

    function nuclideangle_get_xs(this, g, xstype, gout, uvw, mu, i_azi, i_pol) &
         result(xs)
      class(NuclideAngle), intent(in) :: this
      integer, intent(in)              :: g      ! Incoming Energy group
      character(*), intent(in)         :: xstype ! Cross Section Type
      integer, optional, intent(in)    :: gout   ! Outgoing Group
      real(8), optional, intent(in)    :: mu     ! Change in angle
      real(8), optional, intent(in)    :: uvw(3) ! Requested Angle
      integer, optional, intent(in)    :: i_azi  ! Azimuthal Index
      integer, optional, intent(in)    :: i_pol  ! Polar Index
      real(8)                          :: xs     ! Resultant xs

      integer :: i_azi_, i_pol_

      xs = ZERO

      if ((xstype == 'nu_fission' .or. xstype == 'fission' .or. xstype =='chi' &
           .or. xstype =='k_fission') .and. (.not. this % fissionable)) then
        return
      end if

      if (present(i_azi) .and. present(i_pol)) then
        i_azi_ = i_azi
        i_pol_ = i_pol
      else
        call find_angle(this % polar, this % azimuthal, uvw, i_azi_, i_pol_)
      end if

      if (present(gout)) then
        select case(xstype)
        case('mult')
          xs = this % mult(gout,g,i_azi_,i_pol_)
        case('nu_fission')
          xs = this % nu_fission(gout,g,i_azi_,i_pol_)
        case('chi')
          xs = this % chi(gout,i_azi_,i_pol_)
        case('f_mu', 'f_mu/mult')
          xs = this % calc_f(g, gout, mu, I_AZI=i_azi_, I_POL=i_pol_)
          if (xstype == 'f_mu/mult') then
            xs = xs / this % mult(gout,g,i_azi_,i_pol_)
          end if
        end select
      else
        select case(xstype)
        case('total')
          xs = this % total(g,i_azi_,i_pol_)
        case('absorption')
          xs = this % absorption(g,i_azi_,i_pol_)
        case('fission')
          xs = this % fission(g,i_azi_,i_pol_)
        case('k_fission')
          if (allocated(this % k_fission)) then
            xs = this % k_fission(g,i_azi_,i_pol_)
          end if
        case('chi')
          xs = this % chi(g,i_azi_,i_pol_)
        case('scatter')
          xs = this % total(g,i_azi_,i_pol_) - this % absorption(g,i_azi_,i_pol_)
        end select
      end if

    end function nuclideangle_get_xs

!===============================================================================
! NUCLIDE*_CALC_F Finds the value of f(mu), the scattering angle probability,
! given mu
!===============================================================================

    pure function nuclideiso_calc_f(this, gin, gout, mu, uvw, i_azi, i_pol) &
         result(f)
      class(NuclideIso), intent(in) :: this
      integer, intent(in)            :: gin  ! Incoming Energy Group
      integer, intent(in)            :: gout ! Outgoing Energy Group
      real(8), intent(in)            :: mu   ! Angle of interest
      real(8), intent(in), optional  :: uvw(3) ! Direction vector
      integer, intent(in), optional  :: i_azi ! Incoming Energy Group
      integer, intent(in), optional  :: i_pol ! Outgoing Energy Group
      real(8)                        :: f    ! Return value of f(mu)

      real(8) :: dmu, r
      integer :: imu

      if (this % scatt_type == ANGLE_LEGENDRE) then
        f = evaluate_legendre(this % scatter(gout,gin,:), mu)
      else if (this % scatt_type == ANGLE_TABULAR) then
        dmu = TWO / real(this % order - 1)
        ! Find mu bin algebraically, knowing that the spacing is equal
        f   = (mu + ONE) / dmu + ONE
        imu = floor(f)
        ! But save the amount that mu is past the previous index
        ! so we can use interpolation later.
        f = f - real(imu)
        ! Adjust so interpolation works on the last bin if necessary
        if (imu == size(this % scatter, dim=3)) then
          imu = imu - 1
        end if

        ! Now intepolate to find f(mu)
        r  = f / dmu
        f = (ONE - r) * this % scatter(gout,gin,imu) + &
             r * this % scatter(gout,gin,imu+1)
      else ! (ANGLE_HISTOGRAM)
        dmu = TWO / real(this % order)
        ! Find mu bin algebraically, knowing that the spacing is equal
        imu   = floor((mu + ONE) / dmu + ONE)
        ! Adjust so interpolation works on the last bin if necessary
        if (imu == size(this % scatter, dim=3)) then
          imu = imu - 1
        end if
        f = this % scatter(gout, gin, imu)

      end if

    end function nuclideiso_calc_f

    pure function nuclideangle_calc_f(this, gin, gout, mu, uvw, i_azi, &
                                          i_pol) result(f)
      class(NuclideAngle), intent(in) :: this
      integer, intent(in)              :: gin  ! Incoming Energy Group
      integer, intent(in)              :: gout ! Outgoing Energy Group
      real(8), intent(in)              :: mu   ! Angle of interest
      real(8), intent(in), optional    :: uvw(3) ! Direction vector
      integer, intent(in), optional    :: i_azi ! Incoming Energy Group
      integer, intent(in), optional    :: i_pol ! Outgoing Energy Group
      real(8)                          :: f    ! Return value of f(mu)

      real(8) :: dmu, r
      integer :: imu
      integer :: i_azi_, i_pol_
      if (present(i_azi) .and. present(i_pol)) then
        i_azi_ = i_azi
        i_pol_ = i_pol
      else if (present(uvw)) then
        call find_angle(this % polar, this % azimuthal, uvw, i_azi_, i_pol_)
      end if

      if (this % scatt_type == ANGLE_LEGENDRE) then
        f = evaluate_legendre(this % scatter(gout,gin,:,i_azi_,i_pol_), mu)
      else if (this % scatt_type == ANGLE_TABULAR) then
        dmu = TWO / real(this % order - 1)
        ! Find mu bin algebraically, knowing that the spacing is equal
        f   = (mu + ONE) / dmu + ONE
        imu = floor(f)
        ! But save the amount that mu is past the previous index
        ! so we can use interpolation later.
        f = f - real(imu)
        ! Adjust so interpolation works on the last bin if necessary
        if (imu == size(this % scatter, dim=3)) then
          imu = imu - 1
        end if

        ! Now intepolate to find f(mu)
        r  = f / dmu
        f = (ONE - r) * this % scatter(gout,gin,imu,i_azi_,i_pol_) + &
             r * this % scatter(gout,gin,imu+1,i_azi_,i_pol_)
      else ! (ANGLE_HISTOGRAM)
        dmu = TWO / real(this % order)
        ! Find mu bin algebraically, knowing that the spacing is equal
        imu   = floor((mu + ONE) / dmu + ONE)
        ! Adjust so interpolation works on the last bin if necessary
        if (imu == size(this % scatter, dim=3)) then
          imu = imu - 1
        end if
        f = this % scatter(gout, gin, imu,i_azi_,i_pol_)

      end if

    end function nuclideangle_calc_f

!===============================================================================
! find_angle finds the closest angle on the data grid and returns that index
!===============================================================================

    pure subroutine find_angle(polar, azimuthal, uvw, i_azi, i_pol)
      real(8), intent(in) :: polar(:)     ! Polar angles [0,pi]
      real(8), intent(in) :: azimuthal(:) ! Azi. angles [-pi,pi]
      real(8), intent(in) :: uvw(3)       ! Direction of motion
      integer, intent(inout) :: i_pol     ! Closest polar bin
      integer, intent(inout) :: i_azi     ! Closest azi bin

      real(8) my_pol, my_azi, dangle

      ! Convert uvw to polar and azi

      my_pol = acos(uvw(3))
      my_azi = atan2(uvw(2), uvw(1))

      ! Search for equi-binned angles
      dangle = PI / real(size(polar),8)
      i_pol  = floor(my_pol / dangle + ONE)
      dangle = TWO * PI / real(size(azimuthal),8)
      i_azi  = floor((my_azi + PI) / dangle + ONE)

    end subroutine find_angle

end module nuclide_header
