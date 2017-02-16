module URR_isotope

  use URR_constants!TODO
  use URR_cross_sections,    only: CrossSections,&
                                   xs_samples_tmp
  use URR_error,             only: EXIT_SUCCESS,&
                                   EXIT_FAILURE,&
                                   ERROR,&
                                   exit_status,&
                                   log_message
  use URR_faddeeva,          only: faddeeva_w,&
                                   quickw
  use URR_interpolate,       only: interp_factor,&
                                   interpolate
  use URR_openmc_wrapper,    only: ListInt, ListReal, master, prn, binary_search
  use URR_probability_table, only: ProbabilityTable
  use URR_resonance,         only: BreitWignerResonanceListVector1D,&
                                   BreitWignerResonanceVector1D,&
                                   BreitWignerResonanceVector2D,&
                                   ReichMooreResonanceVector1D,&
                                   Resonance,&
                                   wigner_level_spacing
  use URR_settings!TODO
  use URR_vector,            only: VectorReal1D,&
                                   VectorReal2D,&
                                   VectorInt1D

  implicit none
  private
  public :: Isotope,&
            isotopes

!> Type containing data and procedures for processing the URR of an isotope
  type Isotope

    real(8) :: AWR ! weight of nucleus in neutron masses
    real(8) :: T   ! current isotope temperature [K]
    type(ListReal) :: ace_T_list ! list of temperatures isotope has ACE data at
    type(ListInt) :: ace_index_list ! list of indices for different temperatures
    logical :: metastable = .false. ! is isotope metastable?

    ! Current quantum mechanical variables
    real(8) :: k_n       ! wavenumber for neutron energy
    real(8) :: k_lam     ! wavenumber for resonance energy
    real(8) :: P_l_n     ! penetration for neutron energy
    real(8) :: P_l_lam   ! penetration for resonance energy
    real(8) :: S_l_n     ! shift for neutron energy
    real(8) :: phi_l_n   ! phase shift for neutron energy

    ! ENDF-6 nuclear data
    logical :: been_read = .false. ! has ENDF-6 data already been read in?
    integer, allocatable :: NLS(:)  ! number of l values in each energy region
    integer, allocatable :: NJS(:)  ! number of J values for each URR l
    integer, allocatable :: LRU(:)  ! resolved (1) or unresolved (2) parameters
    integer, allocatable :: LRF(:)  ! ENDF-6 resonance formalism number
    integer, allocatable :: NRO(:)  ! AP energy-dependence flag
    integer, allocatable :: NAPS(:) ! channel radius handling flag
    integer :: MAT     ! isotope MAT number
    integer :: ZAI     ! isotope ZAID number
    integer :: LRP     ! resonance parameter flag
    integer :: NER     ! number of resonance energy regions
    integer :: LSSF    ! self-shielding factor flag
    integer :: INT     ! ENDF-6 interpolation law for parameters and xs
    integer :: MF3_INT ! ENDF-6 interpolation law for File 3 xs
! TODO: check for interpolation law consistency
! TODO: use LRX if ZERO's aren't given for competitive width in ENDF when LRX = 0
    integer :: LFW     ! energy-dependent fission width flag
    integer :: LRX     ! competitive width flag
    integer :: NE      ! number of URR tabulated data energies
    real(8) :: E_ex1 = INF ! first level inelastic scattering excitation energy
    real(8) :: E_ex2 = INF ! second level inelastic scattering excitation energy
    real(8), allocatable :: EL(:)  ! lower energy bound of energy region
    real(8), allocatable :: EH(:)  ! upper energy bound of energy region
    real(8), allocatable :: AP(:)  ! scattering radius
    real(8), allocatable :: ac(:)  ! channel radius
    real(8), allocatable :: SPI(:) ! total spin
    real(8), allocatable :: ES(:)  ! URR tabulated data energies

    ! mean values
    type(VectorReal2D), allocatable :: D_mean(:)   ! level spacing
    type(VectorReal2D), allocatable :: GN0_mean(:) ! reduced neutron width
    type(VectorReal2D), allocatable :: GG_mean(:)  ! radiative capture width
    type(VectorReal2D), allocatable :: GF_mean(:)  ! fission width
    type(VectorReal2D), allocatable :: GX_mean(:)  ! competitive width
    type(VectorReal1D), allocatable :: AJ(:)   ! total angular momentum
    type(VectorReal1D), allocatable :: DOFN(:) ! # neutron channels
    type(VectorReal1D), allocatable :: DOFG(:) ! # capture channels
    type(VectorReal1D), allocatable :: DOFF(:) ! # fission channels
    type(VectorReal1D), allocatable :: DOFX(:) ! # competitive channels

    ! current values
    real(8) :: E_last = ZERO ! last neutron energy
    real(8) :: E   ! neutron energy
    real(8) :: J   ! total angular momentum
    real(8) :: g_J ! statistical spin factor
    real(8) :: D   ! level spacing
    real(8) :: GN0 ! reduced neutron width
    real(8) :: GG  ! radiative capture width
    real(8) :: GF  ! fission width
    real(8) :: GX  ! competitive width
    integer :: L    ! orbital quantum number
    integer :: AMUN ! number of neutron channels (degrees of freedom)
    integer :: AMUG ! number of capture channels (degrees of freedom)
    integer :: AMUF ! number of fission channels (degrees of freedom)
    integer :: AMUX ! number of competitive channels (degrees of freedom)

    ! infinite-dilute cross sections computed by averaging
    ! over resonance ladders generated from MF2 and the energy grid
    ! on which they are computed
    integer :: num_avg_xs_grid
    real(8), allocatable :: E_avg_xs(:)
    type(CrossSections), allocatable :: avg_xs(:)

    ! MF3 cross sections and energies
    real(8), allocatable :: MF3_n_e(:)
    real(8), allocatable :: MF3_f_e(:)
    real(8), allocatable :: MF3_g_e(:)
    real(8), allocatable :: MF3_x_e(:)
    real(8), allocatable :: MF3_n(:)
    real(8), allocatable :: MF3_f(:)
    real(8), allocatable :: MF3_g(:)
    real(8), allocatable :: MF3_x(:)

    ! URR treatment parameters and indices
    logical :: otf_urr_xs   = .false. ! calculate URR xs on-the-fly?
    logical :: prob_bands   = .false. ! calculate probability tables?
    logical :: point_urr_xs = .false. ! calculate pointwise URR cross sections?
    integer :: i_urr ! index of URR energy range
    real(8) :: max_E_urr ! max energy for URR treatment

    ! probability tables for given (energy, temperature) pairs
    integer :: nE_tabs ! number of probability table energies
    integer :: nT_tabs ! number of probability table temperatures
    integer :: n_bands ! number of probability table bands
    real(8), allocatable :: E_tabs(:) ! probability table energies
    real(8), allocatable :: T_tabs(:) ! probability table temperatures
    type(ProbabilityTable), allocatable :: prob_tables(:,:)
    type(CrossSections), allocatable :: xs_samples(:,:)

    ! for each (l, realization), a vector of lists of (l, J) resonances
    type(BreitWignerResonanceListVector1D), allocatable :: urr_resonances_tmp(:,:)

    ! for each (l, realization), a vector of vectors of (l, J) resonances
    type(BreitWignerResonanceVector2D), allocatable :: urr_resonances(:,:)

    ! max resonances for a given (l, realization)
    type(VectorInt1D), allocatable :: n_lam(:,:)

    ! vector of Reich-Moore resonances for each l
    type(ReichMooreResonanceVector1D), allocatable :: rm_resonances(:)

    ! vector of BW resonances for each l
    type(BreitWignerResonanceVector1D), allocatable :: bw_resonances(:)

    ! for each (l), a vector of vectors of local (l, J) resonances
    type(BreitWignerResonanceVector2D), allocatable :: local_realization(:)

    ! pointwise URR cross section data
    type(ListReal) :: E_tmp ! scratch energy grid values
    type(ListReal) :: n_tmp ! scratch elastic xs
    type(ListReal) :: g_tmp ! scratch capture xs
    type(ListReal) :: f_tmp ! scratch fission xs
    type(ListReal) :: x_tmp ! scratch competitive xs
    type(ListReal) :: t_tmp ! scratch total xs
    real(8), allocatable :: urr_E(:) ! energy grid values
    type(CrossSections), allocatable :: point_xs(:) ! pointwise, CE cross sections !TODO

  ! Type-Bound procedures
  contains

    ! allocate resonance energy range variables
    procedure :: alloc_energy_range => alloc_energy_range

    ! deallocate resonance energy range variables
    procedure :: dealloc_energy_ranges => dealloc_energy_ranges

    ! deallocate File 3 average (infinite-dilute) cross sections
    procedure :: dealloc_MF3 => dealloc_MF3

    ! allocate URR resonance ensemble realization
    procedure :: alloc_ensemble => alloc_ensemble

    ! deallocate URR resonance ensemble realization
    procedure :: dealloc_ensemble => dealloc_ensemble

    ! allocate local URR resonance realization
    procedure :: alloc_local_realization => alloc_local_realization

    ! deallocate local URR resonance realization
    procedure :: dealloc_local_realization => dealloc_local_realization

    ! allocate probability tables
    procedure :: alloc_prob_tables => alloc_prob_tables

    ! deallocate probability tables
    procedure :: dealloc_prob_tables => dealloc_prob_tables

    ! zeros out statistics accumulators
    procedure :: flush_prob_table_stats => flush_prob_table_stats

    ! zeros out batch accumulators
    procedure :: flush_batches => flush_batches

    ! zeros out xs value accumulator
    procedure :: flush_histories => flush_histories

    ! calculate xs means and standard errors
    procedure :: statistics => statistics

    ! add the contribution of an additional resonance
    procedure :: resonance_contribution => resonance_contribution

    ! update batch accumulators with history result
    procedure :: accum_history => accum_history

    ! accumulate values for a single batch
    procedure :: accum_batch => accum_batch

    ! generate all resonance parameters for URR realizations
    procedure :: resonance_ladder_realization => resonance_ladder_realization

    ! generate pointwise, continuous-energy cross sections
    procedure :: generate_pointwise_xs => generate_pointwise_xs

    ! set channel radius
    procedure :: channel_radius => channel_radius

    ! deallocate isotope
    procedure :: dealloc => dealloc_isotope

    ! add the resonance parameters for a single resonance to the realization
    procedure :: set_parameters => set_parameters

    ! get the URR resonance parameters for a single resonance
    procedure :: get_parameters => get_parameters

    ! get the mean resonance parameters at an energy
    procedure :: get_mean_parameters => get_mean_parameters

    ! sample energy spacing between adjacent resonances
    procedure :: level_spacing => level_spacing

    ! sample channel partial widths
    procedure :: channel_widths => channel_widths

    ! sample resonance parameters for a resonance
    procedure :: resonance_parameters => resonance_parameters

    ! wrapper for calculation of partial cross sections at E_n
    procedure :: xs => calc_xs

    ! calculate SLBW resonance partial cross sections
    procedure :: slbw_xs => slbw_xs

    ! calculate MLBW resonance partial cross sections
    procedure :: mlbw_xs => mlbw_xs

    ! compute G function for MLBW xs
    procedure :: G_func => G_func

    ! compute H function for MLBW xs
    procedure :: H_func => H_func

    ! add potential scattering component to elastic and total xs
    procedure :: potential_xs => potential_xs

    ! add resonance xs component to the evaluator-supplied background xs in MF3
    procedure :: add_mf3 => add_mf3

    ! add resonance xs component to MF3 or multiply MF3 by self-shielding factor
    procedure :: self_shielding => self_shielding ! make sure this isn't being done both when generating prob tables and then again, inline in prob_band_xs

    ! interpolate computed average cross sections
    procedure :: interpolate_avg_xs => interpolate_avg_xs

    ! find energy of last RRR resonance in a given (l,J) spin sequence
    procedure :: last_resolved_resonance_energy => last_resolved_resonance_energy

    ! find index of RRR resonance that contributes to URR xs
    procedure :: resolved_resonance_index => resolved_resonance_index

    ! calculate xs from pre-computed probability tables
    procedure :: prob_band_xs => prob_band_xs

    ! calculate xs from a new realization of resonances on-the-fly
    procedure :: new_realization_otf_xs => new_realization_otf_xs

    ! calculate xs from a fixed realization of resonances on-the-fly
    procedure :: fixed_realization_otf_xs => fixed_realization_otf_xs

    ! generate probability tables for the URR
    procedure :: generate_prob_tables => generate_prob_tables

  end type Isotope

  type(Isotope), allocatable, target :: isotopes(:)

contains


!> Generate resonance parameters for independent realizations of a complete set
!! (ladder) of URR resonances
  subroutine resonance_ladder_realization(this)

    class(Isotope), intent(inout), target :: this ! isotope object

    type(Resonance) :: res ! resonance object
    integer :: iso         ! isotope index
    integer :: i_l         ! orbital quantum number index
    integer :: i_J         ! total angular momentum quantum number
    integer :: i_E         ! tabulated URR parameters energy index
    integer :: i_ens       ! ensemble index
    integer :: i_res       ! l-wave resonance counter
    integer :: n_res       ! number of l-wave resonances to include
    integer :: n_above_urr ! number of resonances abover upper URR energy
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    real(8) :: m     ! energy interpolation factor

    ! allocate vector of linked lists of (l,J) resonances for (l, realization)
    allocate(this % urr_resonances_tmp(this % NLS(this % i_urr), num_urr_realizations))
    allocate(this % n_lam(this % NLS(this % i_urr), num_urr_realizations))

    res % i_res = 0

    ! loop over independent realizations
    ENSEMBLE_LOOP: do i_ens = 1, num_urr_realizations

      ! loop over orbital angular momenta
      ORBITAL_ANG_MOM_LOOP: do i_l = 1, this % NLS(this % i_urr)

        ! alloc vector of NJS(l) linked lists of resonances for (l, realization)
        allocate(this % urr_resonances_tmp(i_l, i_ens) % J(this % NJS(i_l)))
        allocate(this % n_lam(i_l, i_ens) % dim1(this % NJS(i_l)))

        ! set current orbital angular momentum quantum number
        this % L = i_l - 1

        ! get the number of contributing l-wave resonances for this l
        n_res = num_contributing_resonances(this % L)

        ! loop over total angular momenta
        TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)

          ! set current total angular momentum quantum number
          this % J = this % AJ(i_l) % dim1(i_J)

          ! set current partial width degrees of freedom, mean level spacing
          this % AMUX = int(this % DOFX(i_l) % dim1(i_J))
          this % AMUN = int(this % DOFN(i_l) % dim1(i_J))
          this % AMUG = int(this % DOFG(i_l) % dim1(i_J))
          this % AMUF = int(this % DOFF(i_l) % dim1(i_J))
          this % D = this % D_mean(i_l) % dim2(i_J) % dim1(1)

          ! set energy of the lowest-lying contributing URR resonance
          if (this % i_urr == 1) then

            ! the URR is the first resonance energy region so place
            ! resonance energy randomly about lower URR energy bound
            E_res = this % EL(this % i_urr)&
                 + (ONE - TWO * prn()) * wigner_level_spacing(this % D, prn())

          else if (i_l > this % NLS(this % i_urr - 1)) then

            ! the URR has more l-states than the RRR so place resonance energy
            ! randomly about lower URR energy bound
            E_res = this % EL(this % i_urr)&
                 + (ONE - TWO * prn()) * wigner_level_spacing(this % D, prn())

          else

            ! offset first URR resonance energy from the highest-energy RRR
            ! resonance with the same (l,J) spin sequence 
            E_res = this % last_resolved_resonance_energy(this % L, this % J)&
                 + wigner_level_spacing(this % D, prn())

          end if

          ! point to first resonance for this l-wave, realization, J
          this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
               => this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % first

          ! zero resonance counters
          i_res = 0
          n_above_urr = 0
          RESONANCE_LOOP: do while(n_above_urr < n_res + 1)

            i_res = i_res + 1

            ! compute interpolation factor for mean parameters and
            ! don't extrapolate paramters outside of URR because
            ! non-physical quantities can result
            if (E_res < this % ES(1)) then
              i_E = 1
              m = 0
            else if (E_res > this % ES(this % NE)) then
              i_E = this % NE - 1
              m = 1
            else
              i_E = binary_search(this % ES, this % NE, E_res)
              m = interp_factor(&
                   E_res, this % ES(i_E), this % ES(i_E + 1), this % INT)
            end if

            ! set current mean unresolved resonance parameters
            this % D = interpolate(m, &
                 this % D_mean(i_l) % dim2(i_J) % dim1(i_E), &
                 this % D_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
            this % GN0 = interpolate(m, &
                 this % GN0_mean(i_l) % dim2(i_J) % dim1(i_E), &
                 this % GN0_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
            this % GG = interpolate(m, &
                 this % GG_mean(i_l) % dim2(i_J) % dim1(i_E), &
                 this % GG_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
            if (this % INT == LINEAR_LINEAR &
                 .or. this % GF_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
              this % GF = interpolate(m, &
                   this % GF_mean(i_l) % dim2(i_J) % dim1(i_E), &
                   this % GF_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
            else
              this % GF = ZERO
            end if
            if (this % INT == LINEAR_LINEAR &
              .or. this % GX_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
              this % GX = interpolate(m, &
                   this % GX_mean(i_l) % dim2(i_J) % dim1(i_E), &
                   this % GX_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
            else
              this % GX = ZERO
            end if

            ! sample unresolved resonance parameters for this spin
            ! sequence, at this energy
            res % E_lam = E_res
            this % E = E_res
            this % k_n   = wavenumber(this % AWR, abs(this % E))
            this % k_lam = wavenumber(this % AWR, abs(this % E))
            this % P_l_n =&
                 penetration(this % L, this % k_n * this % ac(this % i_urr))
            this % P_l_lam =&
                 penetration(this % L, this % k_lam * this % ac(this % i_urr))
            this % S_l_n =&
                 energy_shift(this % L, this % k_n * this % ac(this % i_urr))
            this % phi_l_n =&
                 phase_shift(this % L, this % k_n * this % AP(this % i_urr))
            call this % channel_widths(res, i_l, i_J)
            call this % set_parameters(res, i_ens, i_l, i_J)

            ! add an additional resonance
            E_res = E_res + wigner_level_spacing(this % D, prn())

            if (E_res > this % EH(this % i_urr)) n_above_urr = n_above_urr + 1

          end do RESONANCE_LOOP
  
          this % n_lam(i_l, i_ens) % dim1(i_J) = i_res

        end do TOTAL_ANG_MOM_LOOP
      end do ORBITAL_ANG_MOM_LOOP
    end do ENSEMBLE_LOOP

    ! allocate URR resonance ensemble realizations
    call this % alloc_ensemble()

    ! transfer linked list URR ensembles to permanent array
    do i_ens = 1, num_urr_realizations
      do i_l = 1, this % NLS(this % i_urr)
        do i_J = 1, this % NJS(i_l)

          ! point to first resonance for this l-wave, realization, J
          this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
               => this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % first
        
          do i_res = 1, this % n_lam(i_l, i_ens) % dim1(i_J)

            ! resonance energy
            this % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % E_lam&
                 = this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % E_lam
          
            ! total angular momentum, J
            this % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % AJ&
                 = this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % AJ
          
            ! neutron width
            this % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GN&
                 = this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GN
          
            ! gamma width
            this % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GG&
                 = this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GG
          
            ! fission width
            this % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GF&
                 = this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GF
          
            ! competitive width
            this % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GX&
                 = this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GX
          
            ! total width
            this % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GT&
                 = this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GT
       
            ! point to next resonance
            this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
                 => this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % next

          end do

          call this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % dealloc()

        end do

        deallocate(this % urr_resonances_tmp(i_l, i_ens) % J)

      end do
    end do

    ! deallocate temporary linked list URR resonance ensemble realizations
    deallocate(this % urr_resonances_tmp)

  end subroutine resonance_ladder_realization


!> Generate pointwise, continuous-energy cross section data in the URR
  subroutine generate_pointwise_xs(this, T)

    class(Isotope), intent(inout) :: this ! isotope object
    real(8), intent(in) :: T ! isotope temperature [K]

    type(Resonance) :: res ! resonance object
    type(CrossSections) :: xs ! object containing all partial xs values
    integer :: i_l       ! orbital quantum number
    integer :: i_J       ! total angular momentum quantum number
    integer :: n_res     ! number of contributing l-state resonances
    integer :: n_rrr_res ! number of RRR resonances we need to grab
    integer :: i_low     ! index of lowest-lying resonance
    integer :: i_res     ! resonance counter
    integer :: i_rrr_res ! RRR resonance index
    integer :: i_ES      ! index of current URR tabulated energy
    integer :: iE        ! pointwise energy grid index
    integer :: i_list       ! index in list of resonances for all (l,J)
    integer :: i_list_start ! index where to start searching for E_lam
    real(8) :: E_lo         ! lower energy grid point for interpolation
    real(8) :: E_hi         ! upper energy grid point for interpolation
    real(8) :: E_mid        ! trial interpolated energy grid point
    real(8) :: xs_lo        ! total xs at E_lo
    real(8) :: xs_hi        ! total xs at E_hi
    real(8) :: xs_mid       ! total xs value at E_mid
    real(8) :: xs_mid_trial ! interpolated xs at E_mid
    real(8) :: rel_err      ! relative error between interpolated and exact xs
    real(8) :: E_0          ! energy at current list index
    real(8) :: E_1          ! energy at next list index
    real(8) :: E_last       ! energy at last grid index

    this % T = T

    ! only one realization allowed when using a pointwise representation
    i_realization = 1

    ! initialize linked list of energies to lower URR bound
    call this % E_tmp % append(this % EL(this % i_urr))

    ! loop over resonances for all spin sequences and add energies
    do i_l = 1, this % NLS(this % i_urr)
      do i_J = 1, this % NJS(i_l)
        if (master)&
             write(*,'(I7,A48,I1,A1,F4.1)') this % ZAI,&
             ': Generating resonances for (l,J) spin sequence ',&
             i_l - 1, ',', this % AJ(i_l) % dim1(i_J)
        i_list_start = 1
        do i_res = 1, this % n_lam(i_l, i_realization) % dim1(i_J)
          do i_list = i_list_start, this % E_tmp % size()
            if (this%urr_resonances(i_l,i_realization)%J(i_J)%res(i_res)%E_lam <=&
                 this % E_tmp % get_item(i_list)) exit
          end do
          call this % E_tmp % insert(i_list,&
               this % urr_resonances(i_l,i_realization) % J(i_J) % res(i_res) % E_lam)
          i_list_start = i_list
        end do
      end do
    end do

    ! clean energy grid of duplicates, values outside URR
    i_list = 1
    E_last = min(this % EH(this % i_urr), this % max_E_urr)
    do
      if (i_list == this % E_tmp % size()) exit
      E_0 = this % E_tmp % get_item(i_list)
      E_1 = this % E_tmp % get_item(i_list+1)
      if (E_0 > E_1) then
        call exit_status(EXIT_FAILURE, 'Pointwise URR energy grid not monotonically increasing')
        return
      end if
      if ((E_0 == E_1)&
           .or. (E_0 > E_last)&
           .or. (E_0 < this % EL(this % i_urr))) then
        call this % E_tmp % remove(E_0)
        cycle
      end if
      i_list = i_list + 1
    end do
    E_0 = this % E_tmp % get_item(i_list)
    if ((E_0 >= E_last) .or. (E_0 < this % EL(this % i_urr)))&
         call this % E_tmp % remove(E_0)
    call this % E_tmp % append(E_last)

    ! set first two energy points
    iE = 1
    E_lo = this % E_tmp % get_item(iE)
    E_hi = this % E_tmp % get_item(iE+1)

    ! calculate xs vals at the first energy
    this % E = E_lo
    this % k_n = wavenumber(this % AWR, abs(this % E))

    ! reset xs accumulators
    call xs % flush()

    ! Get resonance parameters for a local realization about E_n
    ! ----------------------------------------------------------
    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP_1: do i_l = 1, this % NLS(this % i_urr)

      ! set current orbital angular momentum quantum number
      this % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = num_contributing_resonances(this % L)

      ! loop over total quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP_1: do i_J = 1, this % NJS(i_l)

        if (this % E&
             < this % urr_resonances(i_l, i_realization) % J(i_J) % res(1)%E_lam) then
          i_low = 1
        else
          i_low = binary_search(&
               this % urr_resonances(i_l,i_realization) % J(i_J) % res(:) % E_lam,&
               this % n_lam(i_l, i_realization) % dim1(i_J), this % E)
        end if

        ! set current total angular momentum quantum number
        this % J = this % AJ(i_l) % dim1(i_J)

        ! add resonances to this ladder
        res % i_res = 0

        ! if we're near the lower end of the URR, need to incorporate
        ! resolved resonance region resonances in order to fix-up
        ! (i.e. smooth out) cross sections at the RRR-URR crossover
        ! energy
        if (i_low - n_res/2 + 1 < 1) then

          ! if there is an RRR to get resonances from
          if (this % i_urr > 1) then

            ! if the RRR has resonances with this l-state
            if (i_l <= this % NLS(this % i_urr - 1)) then

              ! how many RRR resonances are contributing
              n_rrr_res = abs(i_low - n_res/2)

              ! loop over contributing resolved resonance region resonances
              LOC_RRR_RESONANCES_LOOP_1: do i_res = n_rrr_res, 1, -1

                i_rrr_res = this % resolved_resonance_index(this % L, this % J, i_res)

                ! fewer RRR resonances w/ this J value then needed;
                ! just grab the URR 'edge' resonances that were generated in the RRR
!TODO: take however many RRR resonances there actually are, even if too few
                if (i_rrr_res == 0) exit

                ! add this resolved resonance
                res % i_res = res % i_res + 1
                call this % get_parameters(&
                     res, i_rrr_res, i_l, i_J, this % i_urr - 1)

              end do LOC_RRR_RESONANCES_LOOP_1
            end if
          end if

          ! loop over contributing unresolved resonance region resonances
          LOC_URR_RESONANCES_LOOP_1: do i_res = 1, n_res - res % i_res
            res % i_res = res % i_res + 1
            call this % get_parameters(res, i_res, i_l, i_J, this % i_urr)
          end do LOC_URR_RESONANCES_LOOP_1

        else
          ! we're firmly in the URR and can ignore anything going on in
          ! the upper resolved resonance region energies
          LOC_URR_LOOP_1: do i_res = i_low - n_res/2 + 1, i_low + n_res/2
            res % i_res = res % i_res + 1
            call this % get_parameters(res, i_res, i_l, i_J, this % i_urr)
          end do LOC_URR_LOOP_1
        end if
      end do LOC_TOTAL_ANG_MOM_LOOP_1
    end do LOC_ORBITAL_ANG_MOM_LOOP_1

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP_1: do i_l = 1, this % NLS(this % i_urr)

      ! set current orbital angular momentum quantum number
      this % L = i_l - 1

      ! penetration
      this % P_l_n = penetration(this % L, this % k_n * this % ac(this % i_urr))

      ! resonance energy shift factor
      this % S_l_n = energy_shift(this % L, this % k_n * this % ac(this % i_urr))

      ! hard-sphere phase shift
      this % phi_l_n = phase_shift(this % L, this % k_n * this % AP(this%i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = num_contributing_resonances(this % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP_1: do i_J = 1, this % NJS(i_l)

        ! set current total angular momentum quantum number
        this % J = this % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        this % g_J = (TWO*this % J + ONE) / (FOUR*this % SPI(this%i_urr) + TWO)

        ! loop over resonances localized about E_n
        RESONANCES_LOOP_1: do i_res = 1, n_res

          res % i_res = i_res
          res % E_lam = this%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
          res % Gam_n = this%local_realization(i_l) % J(i_J) % res(i_res)%GN
          res % Gam_g = this%local_realization(i_l) % J(i_J) % res(i_res)%GG
          res % Gam_f = this%local_realization(i_l) % J(i_J) % res(i_res)%GF
          res % Gam_x = this%local_realization(i_l) % J(i_J) % res(i_res)%GX
          res % Gam_t = this%local_realization(i_l) % J(i_J) % res(i_res)%GT

          call this % xs(res)
          call xs % accum_resonance(res % xs_contribution)

        end do RESONANCES_LOOP_1
      end do TOTAL_ANG_MOM_LOOP_1
    end do ORBITAL_ANG_MOM_LOOP_1

    ! add potential scattering contribution
    call this % potential_xs(xs)

    ! combine resonance and File 3 components
    call this % self_shielding(xs)

    ! initialize xs linked lists
    call this % n_tmp % append(xs % n)
    call this % g_tmp % append(xs % g)
    call this % f_tmp % append(xs % f)
    call this % x_tmp % append(xs % x)
    call this % t_tmp % append(xs % t)

    xs_lo = this % t_tmp % get_item(iE)

    ! calculate xs vals at the second energy
    this % E = E_hi
    this % k_n = wavenumber(this % AWR, abs(this % E))

    ! reset xs accumulators
    call xs % flush()

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP_2: do i_l = 1, this % NLS(this % i_urr)

      ! set current orbital angular momentum quantum number
      this % L = i_l - 1

      ! penetration
      this % P_l_n = penetration(this % L, this % k_n * this % ac(this % i_urr))

      ! resonance energy shift factor
      this % S_l_n = energy_shift(this % L, this % k_n * this % ac(this % i_urr))

      ! hard-sphere phase shift
      this % phi_l_n = phase_shift(this % L, this % k_n * this % AP(this%i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = num_contributing_resonances(this % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP_2: do i_J = 1, this % NJS(i_l)

        ! set current total angular momentum quantum number
        this % J = this % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        this % g_J = (TWO*this % J + ONE) / (FOUR*this % SPI(this%i_urr) + TWO)

        ! loop over resonances localized about E_n
        RESONANCES_LOOP_2: do i_res = 1, n_res

          res % i_res = i_res
          res % E_lam = this%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
          res % Gam_n = this%local_realization(i_l) % J(i_J) % res(i_res)%GN
          res % Gam_g = this%local_realization(i_l) % J(i_J) % res(i_res)%GG
          res % Gam_f = this%local_realization(i_l) % J(i_J) % res(i_res)%GF
          res % Gam_x = this%local_realization(i_l) % J(i_J) % res(i_res)%GX
          res % Gam_t = this%local_realization(i_l) % J(i_J) % res(i_res)%GT

          call this % xs(res)
          call xs % accum_resonance(res % xs_contribution)

        end do RESONANCES_LOOP_2
      end do TOTAL_ANG_MOM_LOOP_2
    end do ORBITAL_ANG_MOM_LOOP_2

    ! add potential scattering contribution
    call this % potential_xs(xs)

    ! combine resonance and File 3 components
    call this % self_shielding(xs)

    xs_hi = xs % t

    ! refine energy grid until total cross section convergence
    i_ES = 2
    do
      if (i_ES <= this % NE) then
        if (this % E > this % ES(i_ES - 1)) then
          if (master) write(*,'(I7,A29,ES12.5,A3)') this % ZAI,&
               ': Reconstructing URR xs below', this % ES(i_ES), ' eV'
          i_ES = i_ES + 1
        end if
      end if

      ! compute interpolated energy and xs value
      E_mid = (E_lo + E_hi) / TWO
      xs_mid_trial = (xs_lo + xs_hi) / TWO

      ! calculate xs vals at the interpolated energy
      this % E = E_mid
      this % k_n = wavenumber(this % AWR, abs(this % E))

      ! reset xs accumulators
      call xs % flush()

      ! loop over orbital quantum numbers
      ORBITAL_ANG_MOM_LOOP_3: do i_l = 1, this % NLS(this % i_urr)

        ! set current orbital angular momentum quantum number
        this % L = i_l - 1

        ! penetration
        this % P_l_n = penetration(this % L, this % k_n * this % ac(this % i_urr))

        ! resonance energy shift factor
        this % S_l_n = energy_shift(this % L, this % k_n * this % ac(this % i_urr))

        ! hard-sphere phase shift
        this % phi_l_n = phase_shift(this % L, this % k_n * this % AP(this%i_urr))

        ! get the number of contributing l-wave resonances for this l
        n_res = num_contributing_resonances(this % L)

        ! loop over total angular momentum quantum numbers
        TOTAL_ANG_MOM_LOOP_3: do i_J = 1, this % NJS(i_l)

          ! set current total angular momentum quantum number
          this % J = this % AJ(i_l) % dim1(i_J)

          ! compute statistical spin factor
          this % g_J = (TWO*this % J + ONE) / (FOUR*this % SPI(this%i_urr) + TWO)

          ! loop over resonances localized about E_n
          RESONANCES_LOOP_3: do i_res = 1, n_res

            res % i_res = i_res
            res % E_lam = this%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
            res % Gam_n = this%local_realization(i_l) % J(i_J) % res(i_res)%GN
            res % Gam_g = this%local_realization(i_l) % J(i_J) % res(i_res)%GG
            res % Gam_f = this%local_realization(i_l) % J(i_J) % res(i_res)%GF
            res % Gam_x = this%local_realization(i_l) % J(i_J) % res(i_res)%GX
            res % Gam_t = this%local_realization(i_l) % J(i_J) % res(i_res)%GT

            call this % xs(res)
            call xs % accum_resonance(res % xs_contribution)

          end do RESONANCES_LOOP_3
        end do TOTAL_ANG_MOM_LOOP_3
      end do ORBITAL_ANG_MOM_LOOP_3

      ! add potential scattering contribution
      call this % potential_xs(xs)

      ! combine resonance and File 3 components
      call this % self_shielding(xs)

      xs_mid = xs % t

      ! compute relative error
      if (xs_mid <= ZERO) then
        rel_err = INF
      else
        rel_err = abs((xs_mid_trial - xs_mid) / xs_mid)
      end if

      ! refine energy mesh or accept the point
      if (rel_err <= rel_err_tolerance_pointwise .or. (E_mid - E_lo <= min_delta_E_pointwise)) then
        iE = iE + 1
        if (.not. this % E_tmp % contains(E_mid))&
             call this % E_tmp % insert(iE, E_mid)
        call this % n_tmp % insert(iE, xs % n)
        call this % g_tmp % insert(iE, xs % g)
        call this % f_tmp % insert(iE, xs % f)
        call this % x_tmp % insert(iE, xs % x)
        call this % t_tmp % insert(iE, xs % t)

        ! proceed from the midpoint energy just added
        E_lo = E_mid
        E_hi = this % E_tmp % get_item(iE+1)
        xs_lo = xs_mid

        ! calculate xs vals at the new upper energy
        this % E = E_hi
        this % k_n = wavenumber(this % AWR, abs(this % E))

        ! reset xs accumulators
        call xs % flush()

        ! Get resonance parameters for a local realization about E_n
        ! ----------------------------------------------------------
        ! loop over orbital quantum numbers
        LOC_ORBITAL_ANG_MOM_LOOP_4: do i_l = 1, this % NLS(this % i_urr)

          ! set current orbital angular momentum quantum number
          this % L = i_l - 1

          ! get the number of contributing l-wave resonances for this l
          n_res = num_contributing_resonances(this % L)

          ! loop over total quantum numbers
          LOC_TOTAL_ANG_MOM_LOOP_4: do i_J = 1, this % NJS(i_l)

            if (this % E&
                 < this % urr_resonances(i_l, i_realization) % J(i_J) % res(1)%E_lam) then
              i_low = 1
            else
              i_low = binary_search(&
                   this % urr_resonances(i_l,i_realization) % J(i_J) % res(:) % E_lam,&
                   this % n_lam(i_l, i_realization) % dim1(i_J), this % E)
            end if

            ! set current total angular momentum quantum number
            this % J = this % AJ(i_l) % dim1(i_J)



            ! add resonances to this ladder
            res % i_res = 0

            ! if we're near the lower end of the URR, need to incorporate
            ! resolved resonance region resonances in order to fix-up
            ! (i.e. smooth out) cross sections at the RRR-URR crossover
            ! energy
            if (i_low - n_res/2 + 1 < 1) then

              ! if there is an RRR to get resonances from
              if (this % i_urr > 1) then

                ! if the RRR has resonances with this l-state
                if (i_l <= this % NLS(this % i_urr - 1)) then

                  ! how many RRR resonances are contributing
                  n_rrr_res = abs(i_low - n_res/2)

                  ! loop over contributing resolved resonance region resonances
                  LOC_RRR_RESONANCES_LOOP_4: do i_res = n_rrr_res, 1, -1

                    i_rrr_res = this % resolved_resonance_index(&
                         this % L, this % J, i_res)

                    ! fewer RRR resonances w/ this J value then needed;
                    ! just grab the URR 'edge' resonances that were generated
                    ! in the RRR
!TODO: take however many RRR resonances there actually are, even if too few
                    if (i_rrr_res == 0) exit

                    ! add this resolved resonance
                    res % i_res = res % i_res + 1
                    call this % get_parameters(&
                         res, i_rrr_res, i_l, i_J, this % i_urr - 1)

                  end do LOC_RRR_RESONANCES_LOOP_4
                end if
              end if

              ! loop over contributing unresolved resonance region resonances
              LOC_URR_RESONANCES_LOOP_4: do i_res = 1, n_res - res % i_res
                res % i_res = res % i_res + 1
                call this % get_parameters(res, i_res, i_l, i_J, this % i_urr)
              end do LOC_URR_RESONANCES_LOOP_4

            else
              ! we're firmly in the URR and can ignore anything going on in
              ! the upper resolved resonance region energies
              LOC_URR_LOOP_4: do i_res = i_low - n_res/2 + 1, i_low + n_res/2
                res % i_res = res % i_res + 1
                call this % get_parameters(res, i_res, i_l, i_J, this % i_urr)
              end do LOC_URR_LOOP_4
            end if
          end do LOC_TOTAL_ANG_MOM_LOOP_4
        end do LOC_ORBITAL_ANG_MOM_LOOP_4

        ! loop over orbital quantum numbers
        ORBITAL_ANG_MOM_LOOP_4: do i_l = 1, this % NLS(this % i_urr)

          ! set current orbital angular momentum quantum number
          this % L = i_l - 1

          ! penetration
          this % P_l_n = penetration(this % L, this % k_n * this % ac(this % i_urr))

          ! resonance energy shift factor
          this % S_l_n = energy_shift(this % L, this % k_n * this % ac(this % i_urr))

          ! hard-sphere phase shift
          this % phi_l_n = phase_shift(this % L, this % k_n * this % AP(this%i_urr))

          ! get the number of contributing l-wave resonances for this l
          n_res = num_contributing_resonances(this % L)

          ! loop over total angular momentum quantum numbers
          TOTAL_ANG_MOM_LOOP_4: do i_J = 1, this % NJS(i_l)

            ! set current total angular momentum quantum number
            this % J = this % AJ(i_l) % dim1(i_J)

            ! compute statistical spin factor
            this % g_J = (TWO*this % J + ONE) / (FOUR*this % SPI(this%i_urr) + TWO)

            ! loop over resonances localized about E_n
            RESONANCES_LOOP_4: do i_res = 1, n_res

              res % i_res = i_res
              res % E_lam = this%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
              res % Gam_n = this%local_realization(i_l) % J(i_J) % res(i_res)%GN
              res % Gam_g = this%local_realization(i_l) % J(i_J) % res(i_res)%GG
              res % Gam_f = this%local_realization(i_l) % J(i_J) % res(i_res)%GF
              res % Gam_x = this%local_realization(i_l) % J(i_J) % res(i_res)%GX
              res % Gam_t = this%local_realization(i_l) % J(i_J) % res(i_res)%GT

              call this % xs(res)
              call xs % accum_resonance(res % xs_contribution)

            end do RESONANCES_LOOP_4
          end do TOTAL_ANG_MOM_LOOP_4
        end do ORBITAL_ANG_MOM_LOOP_4

        ! add potential scattering contribution
        call this % potential_xs(xs)

        ! combine resonance and File 3 components
        call this % self_shielding(xs)

        ! add the point if it's the last, otherwise cycle the loop
        if (E_hi >= E_last) then
          if (.not. this % E_tmp % contains(E_hi))&
               call this % E_tmp % insert(iE+1, E_hi)
          call this % n_tmp % insert(iE+1, xs % n)
          call this % g_tmp % insert(iE+1, xs % g)
          call this % f_tmp % insert(iE+1, xs % f)
          call this % x_tmp % insert(iE+1, xs % x)
          call this % t_tmp % insert(iE+1, xs % t)
          exit

        else
          xs_hi = xs % t

        end if
 
      else
        E_hi = E_mid
        xs_hi = xs_mid

      end if
    end do

    ! pass temporary, energy-xs linked lists to dynamic vectors
    allocate(this % urr_E(this % E_tmp % size()))
    allocate(this % point_xs(this % E_tmp % size()))
    do iE = 1, this % E_tmp % size()
      this % urr_E(iE) = this % E_tmp % get_item(iE)
      this % point_xs(iE) % n = this % n_tmp % get_item(iE)
      this % point_xs(iE) % g = this % g_tmp % get_item(iE)
      this % point_xs(iE) % f = this % f_tmp % get_item(iE)
      this % point_xs(iE) % x = this % x_tmp % get_item(iE)
      this % point_xs(iE) % t = this % t_tmp % get_item(iE)
    end do
    call this % E_tmp % clear()
    call this % n_tmp % clear()
    call this % g_tmp % clear()
    call this % f_tmp % clear()
    call this % x_tmp % clear()
    call this % t_tmp % clear()

  end subroutine generate_pointwise_xs


!> Calculate unresolved resonance region cross sections, at a single energy,
!! on-the-fly from resonance parameters for a single realization
  subroutine fixed_realization_otf_xs(this, E, T, xs_in, xs_out)

    class(Isotope), intent(inout) :: this ! isotope object
    real(8), intent(in) :: E ! neutron energy
    real(8), intent(in) :: T ! isotope temperature [K]
    type(CrossSections), intent(in)  :: xs_in  ! application code cross sections
    type(CrossSections), intent(out) :: xs_out ! library output cross sections

    type(Resonance) :: res ! resonance object
    type(CrossSections) :: avg_xs ! computed average cross sections object
    integer :: iavg      ! average cross section index
    integer :: i_l       ! orbital quantum number
    integer :: i_J       ! total angular momentum quantum number
    integer :: n_res     ! number of contributing l-state resonances
    integer :: n_rrr_res ! number of RRR resonances we need to grab
    integer :: i_low     ! index of lowest-lying resonance
    integer :: i_res     ! resonance counter
    integer :: i_rrr_res ! RRR resonance index
    real(8) :: favg ! average cross section interpolation factor

    this % T = T
    this % E = E
    this % k_n = wavenumber(this % AWR, abs(this % E))

    ! reset xs accumulators
    call xs_out % flush()

    ! select which resonance structure realization to use
    if (i_realization_user == 0) then
      ! random realization
      i_realization = 1 + floor(prn() * num_urr_realizations)
      if (i_realization < 1) then
        call exit_status(EXIT_FAILURE, 'i_realization is sampled to be < 1')
        return
      end if
      if (i_realization > num_urr_realizations) then
        call exit_status(EXIT_FAILURE, 'i_realization is sampled to be > num_urr_realizations')
        return
      end if

    else
      ! user-specified realization
      i_realization = i_realization_user

    end if
    
    ! Get resonance parameters for a local realization about E_n
    ! ----------------------------------------------------------
    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, this % NLS(this % i_urr)

      ! set current orbital angular momentum quantum number
      this % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = num_contributing_resonances(this % L)

      ! loop over total angular momentum quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)

        ! set current total angular momentum quantum number
        this % J = this % AJ(i_l) % dim1(i_J)

        ! find the nearest lower resonance
        if (this % E&
             < this % urr_resonances(i_l,i_realization) % J(i_J) % res(1) % E_lam) then
          i_low = 1
        else
          i_low = binary_search(&
               this % urr_resonances(i_l, i_realization) % J(i_J) % res(:) % E_lam,&
               this % n_lam(i_l, i_realization) % dim1(i_J), this % E)
        end if

        ! add resonances to this ladder
        res % i_res = 0

        ! if we're near the lower end of the URR, need to incorporate
        ! resolved resonance region resonances in order to fix-up
        ! (i.e. smooth out) cross sections at the RRR-URR crossover
        ! energy
        if (i_low - n_res/2 + 1 < 1) then

          ! if there is an RRR to get resonances from
          if (this % i_urr > 1) then

            ! if the RRR has resonances with this l-state
            if (i_l <= this % NLS(this % i_urr - 1)) then

              ! how many RRR resonances are contributing
              n_rrr_res = abs(i_low - n_res/2)

              ! loop over contributing resolved resonance region resonances
              LOC_RRR_RESONANCES_LOOP: do i_res = n_rrr_res, 1, -1

                i_rrr_res = this % resolved_resonance_index(this % L, this % J, i_res)

                ! fewer RRR resonances w/ this J value then needed;
                ! just grab the URR 'edge' resonances that were generated in the RRR
!TODO: take however many RRR resonances there actually are, even if too few
                if (i_rrr_res == 0) exit

                ! add this resolved resonance
                res % i_res = res % i_res + 1
                call this % get_parameters(&
                     res, i_rrr_res, i_l, i_J, this % i_urr - 1)

              end do LOC_RRR_RESONANCES_LOOP
            end if
          end if

          ! loop over contributing unresolved resonance region resonances
          LOC_URR_RESONANCES_LOOP: do i_res = 1, n_res - res % i_res
            res % i_res = res % i_res + 1
            call this % get_parameters(res, i_res, i_l, i_J, this % i_urr)
          end do LOC_URR_RESONANCES_LOOP

        else

          ! we're firmly in the URR and can ignore anything going on in
          ! the upper resolved resonance region energies
          LOC_URR_LOOP: do i_res = i_low - n_res/2 + 1, i_low + n_res/2
            res % i_res = res % i_res + 1
            call this % get_parameters(res, i_res, i_l, i_J, this % i_urr)
          end do LOC_URR_LOOP

        end if
      end do LOC_TOTAL_ANG_MOM_LOOP
    end do LOC_ORBITAL_ANG_MOM_LOOP

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP: do i_l = 1, this % NLS(this % i_urr)

      ! set current orbital angular momentum quantum number
      this % L = i_l - 1

      ! penetration
      this % P_l_n = penetration(this % L,&
           this % k_n * this % ac(this % i_urr))
      
      ! resonance energy shift factor
      this % S_l_n = energy_shift(this % L,&
           this % k_n * this % ac(this % i_urr))
      
      ! hard-sphere phase shift
      this % phi_l_n = phase_shift(this % L,&
           this % k_n * this % AP(this % i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = num_contributing_resonances(this % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)

        ! set current total angular momentum quantum number
        this % J = this % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        this % g_J = (TWO * this % J + ONE) &
             / (FOUR * this % SPI(this % i_urr) + TWO)

        ! loop over resonances localized about e_n
        RESONANCES_LOOP: do i_res = 1, n_res

          res % i_res = i_res
          res % E_lam&
               = this % local_realization(i_l) % J(i_J) % res(i_res) % E_lam
          res % Gam_n&
               = this % local_realization(i_l) % J(i_J) % res(i_res) % GN
          res % Gam_g&
               = this % local_realization(i_l) % J(i_J) % res(i_res) % GG
          res % Gam_f&
               = this % local_realization(i_l) % J(i_J) % res(i_res) % GF
          res % Gam_x&
               = this % local_realization(i_l) % J(i_J) % res(i_res) % GX
          res % Gam_t&
               = this % local_realization(i_l) % J(i_J) % res(i_res) % GT

          call this % xs(res)
          call xs_out % accum_resonance(res % xs_contribution)

        end do RESONANCES_LOOP
      end do TOTAL_ANG_MOM_LOOP
    end do ORBITAL_ANG_MOM_LOOP

    ! add potential scattering contribution
    call this % potential_xs(xs_out)

    ! interpret MF3 data according to ENDF-6 LSSF flag
    if (this % LSSF == 0) then
      ! add resonance xs component to MF3 background xs

      call this % add_mf3(xs_out)

    else if (this % LSSF == 1) then
      ! multipy the self-shielding factors by the MF3 background xs

      ! determine index in average xs energy grid
      if (this % E < this % E_avg_xs(1)) then
        iavg = 1

      else if (this % E > this % E_avg_xs(this % num_avg_xs_grid - 1)) then
        iavg = this % num_avg_xs_grid - 1

      else
        iavg = binary_search(this % E_avg_xs, this % num_avg_xs_grid, this % E)

      end if

      favg = interp_factor(&
           E, this % E_avg_xs(iavg), this % E_avg_xs(iavg + 1), this % INT)

      ! interpolate computed averaged URR cross sections
      call this % interpolate_avg_xs(favg, iavg, avg_xs)

      ! competitive xs
      if (avg_xs % x > ZERO&
           .and. this % E <= (ONE + ENDF_PRECISION) * this % E_ex2&
           .and. competitive_structure) then
        ! resonance structure via self-shielding factor if requested
        xs_out % x = xs_out % x / avg_xs % x * xs_in % x
      else
        xs_out % x = xs_in % x
      end if

      ! elastic xs
      xs_out % n = xs_out % n / avg_xs % n * xs_in % n

      ! set negative elastic xs and/or competitive xs to zero
      if (xs_out % n < ZERO) xs_out % n = ZERO
      if (xs_out % x < ZERO) xs_out % x = ZERO

      ! capture xs
      if (avg_xs % g > ZERO) xs_out % g = xs_out % g / avg_xs % g * xs_in % g

      ! fission xs
      if (avg_xs % f > ZERO) xs_out % f = xs_out % f / avg_xs % f * xs_in % f

      ! total xs
      xs_out % t = xs_out % n + xs_out % g + xs_out % f + xs_out % x

    else
      call exit_status(EXIT_FAILURE, 'ENDF-6 LSSF not allowed - must be 0 or 1')
      return

    end if

  end subroutine fixed_realization_otf_xs


!> Generate URR probability tables
  subroutine generate_prob_tables(this, i_isotope)

    class(Isotope), intent(inout), target :: this ! isotope object
    integer, intent(in) :: i_isotope ! URR isotope index

    type(ProbabilityTable), pointer :: ptable ! prob. table pointer
    type(Resonance) :: res ! resonance object
    type(CrossSections) :: xs_potential ! potential cross sections object
    character(6) :: zaid_str ! ZAID number as a string
    integer :: i_b    ! batch index
    integer :: i_band ! probability band index
    integer :: i_E    ! energy grid index
    integer :: i_grid ! File 3 energy grid index
    integer :: i_h    ! history index
    integer :: i_J    ! total angular momentum quantum number index
    integer :: i_l    ! orbital quantum number index
    integer :: i_r    ! resonance index
    integer :: i_T    ! temperature index
    integer :: n_res  ! number of resonances to include for a given l-wave
    integer :: i_mag  ! index in array of xs samples based on total xs magnitude
    integer :: hits_per_band ! number of hits in each equiprobable band
    integer :: avg_unit = 98 ! avg xs output file unit
    integer :: tab_unit = 99 ! tables output file unit
    real(8) :: E        ! neutron lab energy
    real(8) :: fmf3     ! File 3 energy grid interpolation factor
    real(8) :: xs_t_min ! min realized total xs
    real(8) :: xs_t_max ! max realized total xs

    xs_t_min = 1.0e6_8
    xs_t_max = XS_CUTOFF

    write(zaid_str, '(I6)') this % ZAI
    if (master)&
      write(*,*) 'Generating probability tables for ZA '//&
      trim(adjustl(zaid_str))

    if (write_prob_tables) then
      open(unit = tab_unit, file = trim(adjustl(zaid_str))//'-prob-tables.dat')
      write(tab_unit, '("ENDF-6 File:")', advance='no')
      write(tab_unit, *) trim(adjustl(path_endf_files))//trim(adjustl(endf_filenames(i_isotope)))
      write(tab_unit, '("Resonance Formalism:")', advance='no')
      write(tab_unit, *) get_formalism_name(formalism)
      write(tab_unit, '("Contributing s-wave Resonances:")', advance='no')
      write(tab_unit, *) num_l_waves(1)
      write(tab_unit, '("Contributing p-wave Resonances:")', advance='no')
      write(tab_unit, *) num_l_waves(2)
      write(tab_unit, '("Contributing d-wave Resonances:")', advance='no')
      write(tab_unit, *) num_l_waves(3)
      write(tab_unit, '("Contributing f-wave Resonances:")', advance='no')
      write(tab_unit, *) num_l_waves(4)
      write(tab_unit, '("Model Competitive Reaction Resonance Structure:",L2)') competitive_structure
      write(tab_unit, '("Parameter Energy Dependence:")', advance='no')
      write(tab_unit, *) get_energy_dependence(parameter_energy_dependence)
      write(tab_unit, '("Faddeeva Evaluation:")', advance='no')
      write(tab_unit, *) get_faddeeva_method(faddeeva_method)
      write(tab_unit, '("Target Relative Tolerance on Average Partial Cross Sections:")', advance='no')
      write(tab_unit, '(ES24.16)') rel_err_tolerance_avg_xs
      write(tab_unit, '("Energies:")', advance='no')
      write(tab_unit, *) this % nE_tabs
      write(tab_unit, '("Temperatures:")', advance='no')
      write(tab_unit, *) this % nT_tabs
      write(tab_unit, '("Bands:")', advance='no')
      write(tab_unit, *) this % n_bands
    end if

    if (write_avg_xs) then
      open(unit = avg_unit, file = trim(adjustl(zaid_str))//'-avg-urr-xs.dat')
      write(avg_unit, '("ENDF-6 File:")', advance='no')
      write(avg_unit, *) trim(adjustl(path_endf_files))//trim(adjustl(endf_filenames(i_isotope)))
      write(avg_unit, '("Resonance Formalism:")', advance='no')
      write(avg_unit, *) get_formalism_name(formalism)
      write(avg_unit, '("Contributing s-wave Resonances:")', advance='no')
      write(avg_unit, *) num_l_waves(1)
      write(avg_unit, '("Contributing p-wave Resonances:")', advance='no')
      write(avg_unit, *) num_l_waves(2)
      write(avg_unit, '("Contributing d-wave Resonances:")', advance='no')
      write(avg_unit, *) num_l_waves(3)
      write(avg_unit, '("Contributing f-wave Resonances:")', advance='no')
      write(avg_unit, *) num_l_waves(4)
      write(avg_unit, '("Model Competitive Reaction Resonance Structure:",L2)') competitive_structure
      write(avg_unit, '("Parameter Energy Dependence:")', advance='no')
      write(avg_unit, *) get_energy_dependence(parameter_energy_dependence)
      write(avg_unit, '("Faddeeva Evaluation:")', advance='no')
      write(avg_unit, *) get_faddeeva_method(faddeeva_method)
      write(avg_unit, '("Target Relative Tolerance on Average Partial Cross Sections:")', advance='no')
      write(avg_unit, '(ES24.16)') rel_err_tolerance_avg_xs
      write(avg_unit, '("Energies:")', advance='no')
      write(avg_unit, *) this % nE_tabs
      write(avg_unit, '(6A24)')&
           'Energy [eV]', 'Total [b]', 'Elastic [b]', 'Capture [b]', 'Fission [b]','Competitive [b]'
    end if

    ! loop over energy mesh
    ENERGY_LOOP: do i_E = 1, this % nE_tabs

      this % E = this % E_tabs(i_E)
      E = this % E
      this % k_n = wavenumber(this % AWR, abs(this % E))
      this % k_lam = this % k_n 
 
      ! reset accumulator of statistics
      call this % flush_prob_table_stats(i_E)

      i_b = 0

      allocate(this % xs_samples(num_histories_prob_tables, this % nT_tabs))
      this % xs_samples(:,:) % t = ZERO
      this % xs_samples(:,:) % n = ZERO
      this % xs_samples(:,:) % g = ZERO
      this % xs_samples(:,:) % f = ZERO
      this % xs_samples(:,:) % x = ZERO

      ! loop over batches until convergence
      BATCH_LOOP: do

        i_b = i_b + 1

        ! reset batch accumulators
        call this % flush_batches()

        ! loop over realizations
        HISTORY_LOOP: do i_h = 1, num_histories_prob_tables

          ! reset accumulator of histories
          call this % flush_histories()

          ! Get resonance parameters for a local realization about E_n
          ! ----------------------------------------------------------
          ! loop over orbital quantum numbers
          LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, this % NLS(this % i_urr)

            ! set current orbital angular momentum quantum number
            this % L = i_l - 1

            ! penetration
            this % P_l_n = penetration(this % L,&
                 this % k_n * this % ac(this % i_urr))

            ! resonance energy shift factor
            this % S_l_n = energy_shift(this % L,&
                 this % k_n * this % ac(this % i_urr))
            
            ! hard-sphere phase shift
            this % phi_l_n = phase_shift(this % L,&
                 this % k_n * this % AP(this % i_urr))

           ! get the number of contributing l-wave resonances for this l
            n_res = num_contributing_resonances(this % L)

            ! loop over total angular momentum quantum numbers
            LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)

              ! set current total angular momentum quantum number
              this % J = this % AJ(i_l) % dim1(i_J)

              ! set current partial width degrees of freedom
              this % AMUX = int(this % DOFX(i_l) % dim1(i_J))
              this % AMUN = int(this % DOFN(i_l) % dim1(i_J))
              this % AMUG = int(this % DOFG(i_l) % dim1(i_J))
              this % AMUF = int(this % DOFF(i_l) % dim1(i_J))

              ! zero the resonance counter
              res % i_res = 0

              ! set mean URR parameters to neutron energy
              call this % get_mean_parameters(E, i_l, i_J)

              ! sample unresolved resonance parameters for this spin
              ! sequence, at this energy
              this % k_lam = this % k_n
              this % P_l_lam = penetration(this % L,&
                   this % k_lam * this % ac(this % i_urr)) 
              call this % resonance_parameters(res, i_l, i_J)

              ! loop over the addition of resonances to this ladder
              LOC_RESONANCES_LOOP: do i_r = 1, n_res

                if (parameter_energy_dependence == E_RESONANCE) then
                  ! interpolate mean URR parameters to current resonance energy
                  call this % get_mean_parameters(res % E_lam, i_l, i_J)
                  this % k_lam = wavenumber(this % AWR, abs(res % E_lam))
                  this % P_l_lam = penetration(this % L,&
                       this % k_lam * this % ac(this % i_urr))
                end if

                res % i_res = i_r

                ! sample unresolved resonance parameters for this spin
                ! sequence, at this energy
                call this % resonance_parameters(res, i_l, i_J)

              end do LOC_RESONANCES_LOOP
            end do LOC_TOTAL_ANG_MOM_LOOP
          end do LOC_ORBITAL_ANG_MOM_LOOP

          ! loop over orbital quantum numbers
          ORBITAL_ANG_MOM_LOOP: do i_l = 1, this % NLS(this % i_urr)

            ! set current orbital angular momentum quantum number
            this % L = i_l - 1

            ! penetration
            this % P_l_n = penetration(this % L,&
                 this % k_n * this % ac(this % i_urr))
      
            ! resonance energy shift factor
            this % S_l_n = energy_shift(this % L,&
                 this % k_n * this % ac(this % i_urr))
      
            ! hard-sphere phase shift
            this % phi_l_n = phase_shift(this % L,&
                 this % k_n * this % AP(this % i_urr))

            ! get the number of contributing l-wave resonances for this l
            n_res = num_contributing_resonances(this % L)

            ! loop over total angular momentum quantum numbers
            TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)

              ! set current total angular momentum quantum number
              this % J = this % AJ(i_l) % dim1(i_J)

              ! compute statistical spin factor
              this % g_J = (TWO * this % J + ONE) &
                   / (FOUR * this % SPI(this % i_urr) + TWO)

              ! loop over resonances localized about e_n
              RESONANCES_LOOP: do i_r = 1, n_res

                res % i_res = i_r
                res % E_lam&
                     = this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
                res % Gam_n&
                     = this % local_realization(i_l) % J(i_J) % res(i_r) % GN
                res % Gam_g&
                     = this % local_realization(i_l) % J(i_J) % res(i_r) % GG
                res % Gam_f&
                     = this % local_realization(i_l) % J(i_J) % res(i_r) % GF
                res % Gam_x&
                     = this % local_realization(i_l) % J(i_J) % res(i_r) % GX
                res % Gam_t&
                     = this % local_realization(i_l) % J(i_J) % res(i_r) % GT

                TEMPERATURES_LOOP: do i_T = 1, this % nT_tabs

                  ! set current temperature
                  this % T = this % T_tabs(i_T)

                  ! calculate the contribution to the partial cross sections,
                  ! at this energy, from an additional resonance
                  call this % xs(res)

                  ! add this contribution to the accumulated partial cross
                  ! section values built up from all resonances
! TODO: consider moving t outside of loop
                  call this % resonance_contribution(res, i_E, i_T)

                end do TEMPERATURES_LOOP
              end do RESONANCES_LOOP
            end do TOTAL_ANG_MOM_LOOP
          end do ORBITAL_ANG_MOM_LOOP

          TEMPERATURES_LOOPb: do i_T = 1, this % nT_tabs

            ! add potential scattering contribution
            call xs_potential % flush()
            call this % potential_xs(xs_potential)
            this % prob_tables(i_E, i_T) % avg_n % xs&
                 = this % prob_tables(i_E, i_T) % avg_n % xs + xs_potential % n
            this % prob_tables(i_E, i_T) % avg_t % xs&
                 = this % prob_tables(i_E, i_T) % avg_t % xs + xs_potential % t

            ! set negative elastic xs to zero
            if (this % prob_tables(i_E, i_T) % avg_n % xs < ZERO) then
              this % prob_tables(i_E, i_T) % avg_t % xs &
                   = this % prob_tables(i_E, i_T) % avg_t % xs &
                   + abs(this % prob_tables(i_E, i_T) % avg_n % xs)
              this % prob_tables(i_E, i_T) % avg_n % xs = ZERO
            end if

            ! use MF3 fission cross section if requested
            if (background_xs_treatment == FALSE) then
              continue
            else
              this % prob_tables(i_E, i_T) % avg_t % xs &
                   = this % prob_tables(i_E, i_T) % avg_t % xs &
                   - this % prob_tables(i_E, i_T) % avg_f % xs
              this % prob_tables(i_E, i_T) % avg_f % xs = ZERO
              if (background_xs_treatment == ENDFFILE .and. allocated(this % MF3_f_e)) then
                if (E >= this % MF3_f_e(1)) then
                  i_grid = binary_search(this % MF3_f_e,size(this % MF3_f_e),E)
                  if (this % INT == LINEAR_LINEAR &
                       .or. (this % MF3_f(i_grid) > XS_CUTOFF &
                       .and. this % MF3_f(i_grid + 1) > XS_CUTOFF)) then
                    fmf3 = interp_factor(E, this % MF3_f_e(i_grid), &
                         this % MF3_f_e(i_grid + 1), this % INT)
                    this % prob_tables(i_E, i_T) % avg_f % xs &
                         = interpolate(fmf3, this % MF3_f(i_grid), &
                         this % MF3_f(i_grid + 1), this % INT)
                    this % prob_tables(i_E, i_T) % avg_t % xs &
                         = this % prob_tables(i_E, i_T) % avg_t % xs &
                         + this % prob_tables(i_E, i_T) % avg_f % xs
                  end if
                end if
              end if
            end if
            if (this % prob_tables(i_E, i_T) % avg_f % xs < ZERO) then
              call exit_status(EXIT_FAILURE, 'Negative fission xs encountered')
              return
            end if
            ! use MF3 competitive cross section if requested
            if (background_xs_treatment == FALSE) then
              continue
            else
              this % prob_tables(i_E, i_T) % avg_t % xs &
                   = this % prob_tables(i_E, i_T) % avg_t % xs &
                   - this % prob_tables(i_E, i_T) % avg_x % xs
              this % prob_tables(i_E, i_T) % avg_x % xs = ZERO
              if (background_xs_treatment == ENDFFILE .and. allocated(this % MF3_x_e)) then
                if (E >= this % MF3_x_e(1)) then
                  i_grid = binary_search(this % MF3_x_e,size(this % MF3_x_e),E)
                  if (this % INT == LINEAR_LINEAR &
                       .or. (this % MF3_x(i_grid) > XS_CUTOFF &
                       .and. this % MF3_x(i_grid + 1) > XS_CUTOFF)) then
                    fmf3 = interp_factor(E, this % MF3_x_e(i_grid), &
                         this % MF3_x_e(i_grid + 1), this % INT)
                    this % prob_tables(i_E, i_T) % avg_x % xs &
                         = interpolate(fmf3, this % MF3_x(i_grid), &
                         this % MF3_x(i_grid + 1), this % INT)
                    this % prob_tables(i_E, i_T) % avg_t % xs &
                         = this % prob_tables(i_E, i_T) % avg_t % xs &
                         + this % prob_tables(i_E, i_T) % avg_x % xs
                  end if
                end if
              end if
            end if

            ! set negative competitive xs to zero
            if (this % prob_tables(i_E, i_T) % avg_x % xs < ZERO) then
              this % prob_tables(i_E, i_T) % avg_t % xs &
                   = this % prob_tables(i_E, i_T) % avg_t % xs &
                   + abs(this % prob_tables(i_E, i_T) % avg_x % xs)
              this % prob_tables(i_E, i_T) % avg_x % xs = ZERO
            end if

            ! Set min and max xs values encountered
            if (this % prob_tables(i_E, i_T) % avg_t % xs < xs_t_min)&
                 xs_t_min = this % prob_tables(i_E, i_T) % avg_t % xs
            if (this % prob_tables(i_E, i_T) % avg_t % xs > xs_t_max)&
                 xs_t_max = this % prob_tables(i_E, i_T) % avg_t % xs

            ! accumulate the result of this history
            call this % accum_history(i_E, i_T)

            ! find index where this sample belongs based on total cross section
            if (i_b == 1 .and. i_h == 1) then
              i_mag = num_histories_prob_tables
            else
              if (this % prob_tables(i_E, i_T) % avg_t % xs&
                   > this % xs_samples(i_b * num_histories_prob_tables, i_T) % t) then
                i_mag = i_b * num_histories_prob_tables
              else
               i_mag = binary_search(this % xs_samples(:, i_T) % t,&
                     i_b * num_histories_prob_tables,&
                     this % prob_tables(i_E, i_T) % avg_t % xs)
              end if
            end if

            ! insert total and conditional partial cross sections
            this%xs_samples(1 : i_mag - 1, i_T) = this % xs_samples(2:i_mag,i_T)
            this%xs_samples(i_mag,i_T) % t = this%prob_tables(i_E,i_T)% avg_t%xs
            this%xs_samples(i_mag,i_T) % n = this%prob_tables(i_E,i_T)% avg_n%xs
            this%xs_samples(i_mag,i_T) % g = this%prob_tables(i_E,i_T)% avg_g%xs
            this%xs_samples(i_mag,i_T) % f = this%prob_tables(i_E,i_T)% avg_f%xs
            this%xs_samples(i_mag,i_T) % x = this%prob_tables(i_E,i_T)% avg_x%xs
          
          end do TEMPERATURES_LOOPb
        end do HISTORY_LOOP

        call move_alloc(this % xs_samples, xs_samples_tmp)
        allocate(this % xs_samples((i_b + 1) * num_histories_prob_tables, this % nT_tabs))
        this % xs_samples(:,:) % t = ZERO
        this % xs_samples(:,:) % n = ZERO
        this % xs_samples(:,:) % g = ZERO
        this % xs_samples(:,:) % f = ZERO
        this % xs_samples(:,:) % x = ZERO
        this % xs_samples(num_histories_prob_tables+1 : (i_b+1)*num_histories_prob_tables, :)&
             = xs_samples_tmp
        deallocate(xs_samples_tmp)
        
        ! accumulate the result of this batch
        call this % accum_batch(i_E)

        ! calculate statistics for this batch
        call this % statistics(i_b)

        if ((i_b > min_num_batches_prob_tables &
             .and. max(maxval(this % prob_tables(i_E, :) % avg_t % rel_unc), &
                       maxval(this % prob_tables(i_E, :) % avg_n % rel_unc), &
                       maxval(this % prob_tables(i_E, :) % avg_g % rel_unc), &
                       maxval(this % prob_tables(i_E, :) % avg_f % rel_unc), &
                       maxval(this % prob_tables(i_E, :) % avg_x % rel_unc)) &
             < rel_err_tolerance_avg_xs) .or. i_b == max_num_batches_prob_tables) exit

      end do BATCH_LOOP

      call move_alloc(this % xs_samples, xs_samples_tmp)
      allocate(this % xs_samples((i_b) * num_histories_prob_tables, this % nT_tabs))
      this % xs_samples(:,:)&
           = xs_samples_tmp(num_histories_prob_tables+1:(i_b+1)*num_histories_prob_tables,:)
      deallocate(xs_samples_tmp)
      do i_T = 1, this % nT_tabs
        hits_per_band = nint(i_b * num_histories_prob_tables / dble(this % n_bands))
        do i_band = 1, this % n_bands - 1
          this % prob_tables(i_E, i_T) % t(i_band) % cnt_mean&
               = dble(hits_per_band) / dble(i_b * num_histories_prob_tables)
          this % prob_tables(i_E, i_T) % t(i_band) % xs_mean&
               = sum(this % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % t) / dble(hits_per_band)
          this % prob_tables(i_E, i_T) % n(i_band) % xs_mean&
               = sum(this % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % n) / dble(hits_per_band)
          this % prob_tables(i_E, i_T) % g(i_band) % xs_mean&
               = sum(this % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % g) / dble(hits_per_band)
          this % prob_tables(i_E, i_T) % f(i_band) % xs_mean&
               = sum(this % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % f) / dble(hits_per_band)
          this % prob_tables(i_E, i_T) % x(i_band) % xs_mean&
               = sum(this % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % x) / dble(hits_per_band)
        end do
        this % prob_tables(i_E, i_T) % t(this % n_bands) % cnt_mean&
             = dble(i_b * num_histories_prob_tables - (this % n_bands - 1) * hits_per_band)&
             / dble(i_b * num_histories_prob_tables)
        this % prob_tables(i_E, i_T) % t(this % n_bands) % xs_mean&
             = sum(this % xs_samples((this % n_bands - 1) * hits_per_band + 1&
             : i_b * num_histories_prob_tables, i_T) % t)&
             /dble(i_b * num_histories_prob_tables - (this % n_bands - 1) * hits_per_band) 
        this % prob_tables(i_E, i_T) % n(this % n_bands) % xs_mean&
             = sum(this % xs_samples((this % n_bands - 1) * hits_per_band + 1&
             : i_b * num_histories_prob_tables, i_T) % n)&
             /dble(i_b * num_histories_prob_tables - (this % n_bands - 1) * hits_per_band) 
        this % prob_tables(i_E, i_T) % g(this % n_bands) % xs_mean&
             = sum(this % xs_samples((this % n_bands - 1) * hits_per_band + 1&
             : i_b * num_histories_prob_tables, i_T) % g)&
             /dble(i_b * num_histories_prob_tables - (this % n_bands - 1) * hits_per_band) 
        this % prob_tables(i_E, i_T) % f(this % n_bands) % xs_mean&
             = sum(this % xs_samples((this % n_bands - 1) * hits_per_band + 1&
             : i_b * num_histories_prob_tables, i_T) % f)&
             /dble(i_b * num_histories_prob_tables - (this % n_bands - 1) * hits_per_band) 
        this % prob_tables(i_E, i_T) % x(this % n_bands) % xs_mean&
             = sum(this % xs_samples((this % n_bands - 1) * hits_per_band + 1&
             : i_b * num_histories_prob_tables, i_T) % x)&
             /dble(i_b * num_histories_prob_tables - (this % n_bands - 1) * hits_per_band) 
      end do

      deallocate(this % xs_samples)

      ! write probability tables out to a file
      if (write_prob_tables) then
        if (master)&
          write(*,'(A32,ES10.3,A3)') 'Generated probability tables at', &
               this % E_tabs(i_E), ' eV'
        do i_T = 1, this % nT_tabs
          this % T = this % T_tabs(i_T)
          write(tab_unit, '(A24,ES24.16)') 'E [eV]', this % E_tabs(i_E)
          write(tab_unit, '(A24,ES24.16)') 'T [K]', this % T
          write(tab_unit, '(8A24)') 'Min Total [b]', 'Max Total [b]',&
               'Probability', 'Total [b]', 'Elastic [b]', 'Capture [b]', 'Fission [b]', 'Competitive [b]'
          ptable => this % prob_tables(i_E, i_T)
          do i_band = 1, this % n_bands
            if (i_band == 1) then
              write(tab_unit, '(ES24.16,A24,6ES24.16)')&
                   xs_t_min, '',&
                   ptable % t(i_band) % cnt_mean,&
                   ptable % t(i_band) % xs_mean,&
                   ptable % n(i_band) % xs_mean,&
                   ptable % g(i_band) % xs_mean,&
                   ptable % f(i_band) % xs_mean,&
                   ptable % x(i_band) % xs_mean
            else if (i_band == this % n_bands) then
              write(tab_unit, '(A24,7ES24.16)')&
                   '', xs_t_max,&
                   ptable % t(i_band) % cnt_mean,&
                   ptable % t(i_band) % xs_mean,&
                   ptable % n(i_band) % xs_mean,&
                   ptable % g(i_band) % xs_mean,&
                   ptable % f(i_band) % xs_mean,&
                   ptable % x(i_band) % xs_mean
            else
              write(tab_unit, '(2A24,6ES24.16)')&
                   '', '',&
                   ptable % t(i_band) % cnt_mean,&
                   ptable % t(i_band) % xs_mean,&
                   ptable % n(i_band) % xs_mean,&
                   ptable % g(i_band) % xs_mean,&
                   ptable % f(i_band) % xs_mean,&
                   ptable % x(i_band) % xs_mean
            end if
          end do
          write(tab_unit,'(A24,I24)') 'Batches', i_b
          write(tab_unit,'(24X,A24,24X,5ES24.16)')&
               'Mean [b]',&
               ptable % avg_t % xs_mean,&
               ptable % avg_n % xs_mean,&
               ptable % avg_g % xs_mean,&
               ptable % avg_f % xs_mean,&
               ptable % avg_x % xs_mean
          write(tab_unit,'(24X,A24,24X,5ES24.16)')&
               '1sigma [b]',&
               ptable % avg_t % xs_sem,&
               ptable % avg_n % xs_sem,&
               ptable % avg_g % xs_sem,&
               ptable % avg_f % xs_sem,&
               ptable % avg_x % xs_sem
          write(tab_unit,'(24X,A24,24X,5ES24.16)')&
               '1sigma/Mean',&
               ptable % avg_t % rel_unc,&
               ptable % avg_n % rel_unc,&
               ptable % avg_g % rel_unc,&
               ptable % avg_f % rel_unc,&
               ptable % avg_x % rel_unc
        end do
      end if

      ! write averaged URR cross sections out to a file
      if (write_avg_xs) then
        ptable => this % prob_tables(i_E, 1)
        write(avg_unit, '(6ES24.16)')&
             this % E_tabs(i_E),&
             ptable % avg_t % xs_mean,&
             ptable % avg_n % xs_mean,&
             ptable % avg_g % xs_mean,&
             ptable % avg_f % xs_mean,&
             ptable % avg_x % xs_mean
      end if

    end do ENERGY_LOOP

    ! write max 1sigma/mean values
    if (write_avg_xs) then
      write(avg_unit, '(A24,5ES24.16)')&
           'Max 1sigma/Mean',&
           maxval(this % prob_tables(:, 1) % avg_t % rel_unc),&
           maxval(this % prob_tables(:, 1) % avg_n % rel_unc),&
           maxval(this % prob_tables(:, 1) % avg_g % rel_unc),&
           maxval(this % prob_tables(:, 1) % avg_f % rel_unc),&
           maxval(this % prob_tables(:, 1) % avg_x % rel_unc)
    end if

    if (write_prob_tables) close(tab_unit)
    if (write_avg_xs) close(avg_unit)

    nullify(ptable)

  end subroutine generate_prob_tables


!> Allocate variables for the resonance energy ranges
  subroutine alloc_energy_range(this)

    class(Isotope), intent(inout) :: this ! isotope object

    allocate(this % NLS(this % NER))
    allocate(this % LRU(this % NER))
    allocate(this % LRF(this % NER))
    allocate(this % NRO(this % NER))
    allocate(this % NAPS(this % NER))
    allocate(this % EL(this % NER))
    allocate(this % EH(this % NER))
    allocate(this % AP(this % NER))
    allocate(this % ac(this % NER))
    allocate(this % SPI(this % NER))

  end subroutine alloc_energy_range


!> Deallocate variables for the resonance energy ranges
  subroutine dealloc_energy_ranges(this)

    class(Isotope), intent(inout) :: this ! isotope object

    deallocate(this % NLS)
    deallocate(this % LRU)
    deallocate(this % LRF)
    deallocate(this % NRO)
    deallocate(this % NAPS)
    deallocate(this % EL)
    deallocate(this % EH)
    deallocate(this % AP)
    deallocate(this % ac)
    deallocate(this % SPI)

  end subroutine dealloc_energy_ranges


!> Allocate URR resonance ensemble realizations
  subroutine alloc_ensemble(this)

    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_ens ! realization/ensemble index
    integer :: i_l   ! orbital angular momentum index
    integer :: i_J   ! total angular momentum index

    ! allocate URR resonance realizations
    allocate(this % urr_resonances(this % NLS(this % i_urr), num_urr_realizations))

    ! loop over realizations
    do i_ens = 1, num_urr_realizations

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        allocate(this % urr_resonances(i_l, i_ens) % J(this % NJS(i_l)))

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)

          ! allocate URR resonances
          allocate(this % urr_resonances(i_l, i_ens) % J(i_J) %&
               res(this % n_lam(i_l, i_ens) % dim1(i_J)))

        end do
      end do
    end do

  end subroutine alloc_ensemble


!> Deallocate a URR resonance ensemble realization
  subroutine dealloc_ensemble(this)

    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_ens ! realization/ensemble index
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum index

    ! loop over realizations
    do i_ens = 1, num_urr_realizations

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)
! TODO: deallocate resonances of the proper formalism once MLBW, RM allowed
          ! deallocate resonance parameters
          deallocate(this % urr_resonances(i_l, i_ens) % J(i_J) % res)
      
        end do

        deallocate(this % urr_resonances(i_l, i_ens) % J)

      end do
    end do

    ! deallocate URR resonance realizations
    deallocate(this % urr_resonances)

  end subroutine dealloc_ensemble


!> Allocate probability tables for this isotope
  subroutine alloc_prob_tables(this)

    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    if (E_grid_scheme_prob_tables == ENDF6) then
      this % nE_tabs = this % NE
      allocate(this % E_tabs(this % nE_tabs))
      this % E_tabs(:) = this % ES
    else if (E_grid_scheme_prob_tables == LOGARITHMIC) then
      this % nE_tabs = num_energies_prob_tables
      allocate(this % E_tabs(this % nE_tabs))
      this % E_tabs(:) = (/(this % EL(this % i_urr) &
           * exp(i_E * log(this % EH(this % i_urr) / this % EL(this % i_urr)) &
           / (this % nE_tabs - 1)), i_E = 0, this % nE_tabs - 1)/)
    else if (E_grid_scheme_prob_tables == USER) then
      this % nE_tabs = num_energies_prob_tables
      allocate(this % E_tabs(this % nE_tabs))
      this % E_tabs(:) = E_grid_prob_tables
    end if

    this % nT_tabs = num_temperatures_prob_tables
    allocate(this % T_tabs(this % nT_tabs))
    this % T_tabs(:) = T_grid_prob_tables
    this % n_bands = num_bands_prob_tables

    allocate(this % prob_tables(this % nE_tabs, this % nT_tabs))

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
        allocate(this % prob_tables(i_E, i_T) % t(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % n(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % g(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % f(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % x(this % n_bands))
      end do
    end do

  end subroutine alloc_prob_tables


!> Deallocate probability tables for this isotope
  subroutine dealloc_prob_tables(this)

    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
        deallocate(this % prob_tables(i_E, i_T) % t)
        deallocate(this % prob_tables(i_E, i_T) % n)
        deallocate(this % prob_tables(i_E, i_T) % g)
        deallocate(this % prob_tables(i_E, i_T) % f)
        deallocate(this % prob_tables(i_E, i_T) % x)
      end do
    end do

    deallocate(this % E_tabs)
    deallocate(this % prob_tables)

  end subroutine dealloc_prob_tables


!> Zero out statistics accumulators for an isotope's probability tables
  subroutine flush_prob_table_stats(this, i_E)

    class(Isotope), intent(inout) :: this ! isotope object
    integer, intent(in) :: i_E ! tabulated URR energy index

    integer :: i_T ! temperature index
    integer :: i_b ! probability band index

    do i_T = 1, this % nT_tabs
      do i_b = 1, this % n_bands
        call this % prob_tables(i_E, i_T) % t(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % n(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % g(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % f(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % x(i_b) % flush_xs_stats()
      end do
      call this % prob_tables(i_E, i_T) % avg_t % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_n % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_g % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_f % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_x % flush_xs_stats()
    end do

  end subroutine flush_prob_table_stats


!> Compute or set the channel radius depending on ENDF flags
  subroutine channel_radius(this, i_ER)

    class(Isotope), intent(inout) :: this ! isotope object
    integer, intent(in) :: i_ER ! resonance energy range index

    select case (this % NRO(i_ER))

    ! scattering radius is independent of energy
    case (0)

      select case (this % NAPS(i_ER))

      ! use channel radius for penetrabilities and shift factors but
      ! scattering radius for phase shifts
      case (0)
        this % ac(i_ER) = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

      ! use scattering radius for penetrabilities, shift factors and phase
      ! shifts
      case (1)
        this % ac(i_ER) = this % AP(i_ER)

      ! invalid scattering radius treatment flag
      case default
        call exit_status(EXIT_FAILURE, 'ENDF-6 NAPS flag must be 0 or 1 when NRO is 0')
        return

      end select

    ! scattering radius is energy dependent
    case (1)

      select case (this % NAPS(i_ER))

      ! use channel radius for penetrabilities and shift factors but
      ! scattering radius for phase shifts
      case (0)
        this % ac(i_ER) = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

      ! use scattering radius for penetrabilities, shift factors and phase
      ! shifts
      case (1)
        this % ac(i_ER) = this % AP(i_ER)

! TODO: understand this and implement it correctly
      ! use energy dependent scattering radius in phase shifts but the energy
      ! independent AP scattering radius value for penetrabilities and shift
      ! factors
      case (2)
        this % ac(i_ER) = this % AP(i_ER)

      ! invalid scattering radius treatment flag
      case default
        call exit_status(EXIT_FAILURE, 'ENDF-6 NAPS flag must be 0, 1, or 2 when NRO is 1')
        return
      
      end select

    ! invalid energy dependence of scattering radius flag
    case default
      call exit_status(EXIT_FAILURE, 'ENDF-6 NRO flag must be 0 or 1')
      return
    end select

  end subroutine channel_radius


!> Allocate an NLS-length vector of NJS(l)-length vectors of NRS(l,J)-length
!! vectors of spin group resonances
  subroutine alloc_local_realization(this)
  
    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_l ! orbital angular momentum index 
    integer :: i_J ! total angular momentum index

    ! allocate NLS-length vec of NJS(l)-length vec of NRS(l,J)-length vec
    allocate(this % local_realization(this % NLS(this % i_urr)))

    ! loop over orbital quantum numbers
    do i_l = 1, this % NLS(this % i_urr)

      ! allocate each vector of (l,J) vectors
      allocate(this % local_realization(i_l) % J(this % NJS(i_l)))

      ! loop of total quantum numbers
      do i_J = 1, this % NJS(i_l)

        ! allocate each vector of (l,J) resonances
        allocate(this%local_realization(i_l)%J(i_J)%res(num_contributing_resonances(i_l-1)))

      end do
    end do

  end subroutine alloc_local_realization


!> Deallocate an NLS-length vector of NJS(l)-length vectors of NRS(l,J)-length
!! vectors of spin group resonances
  subroutine dealloc_local_realization(this)

    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_l ! orbital angular momentum index 
    integer :: i_J ! total angular momentum index

    ! loop over orbital quantum numbers
    do i_l = 1, this % NLS(this % i_urr)

      ! loop of total quantum numbers
      do i_J = 1, this % NJS(i_l)

        ! deallocate each vector of (l,J) resonances
        call this % local_realization(i_l) % J(i_J) % dealloc()

      end do
      
      ! deallocate each vector of (l,J) vectors
      deallocate(this % local_realization(i_l) % J)

    end do

    ! deallocate NLS-length vec of NJS(l)-length vec of NRS(l,J)-length vec
    deallocate(this % local_realization)

  end subroutine dealloc_local_realization


!> Deallocate isotope object
  subroutine dealloc_isotope(this)

    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_l   ! orbital angular momentum index
    integer :: i_J   ! total angular momentum index
    integer :: i_ens ! realization index

    ! deallocate mean parameters
    do i_l = 1, this % NLS(this % i_urr)
      do i_J = 1, this % NJS(i_l)
        deallocate(this % D_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GN0_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GG_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GF_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GX_mean(i_l) % dim2(i_J) % dim1)
      end do
      if (allocated(this % n_lam)) then
        do i_ens = 1, num_urr_realizations
          deallocate(this % n_lam(i_l, i_ens) % dim1)
        end do
      end if
      deallocate(this % D_mean(i_l) % dim2)
      deallocate(this % GN0_mean(i_l) % dim2)
      deallocate(this % GG_mean(i_l) % dim2)
      deallocate(this % GF_mean(i_l) % dim2)
      deallocate(this % GX_mean(i_l) % dim2)
      deallocate(this % AJ(i_l) % dim1)
      deallocate(this % DOFX(i_l) % dim1)
      deallocate(this % DOFN(i_l) % dim1)
      deallocate(this % DOFG(i_l) % dim1)
      deallocate(this % DOFF(i_l) % dim1)
    end do
    deallocate(this % D_mean)
    deallocate(this % GN0_mean)
    deallocate(this % GG_mean)
    deallocate(this % GF_mean)
    deallocate(this % GX_mean)
    deallocate(this % AJ)
    deallocate(this % DOFX)
    deallocate(this % DOFN)
    deallocate(this % DOFG)
    deallocate(this % DOFF)
    deallocate(this % ES)

    ! deallocate probability tables
    if (allocated(this % prob_tables)) call this % dealloc_prob_tables()

    ! deallocate URR resonances
    if (allocated(this % n_lam)) deallocate(this % n_lam)
    if (allocated(this % urr_resonances)) call this % dealloc_ensemble()

    ! deallocate RRR resonances
    if (allocated(this % bw_resonances)) then
      do i_l = 1, this % NLS(this % i_urr - 1)
        call this % bw_resonances(i_l) % dealloc()
      end do
      deallocate(this % bw_resonances)
    end if
    if (allocated(this % rm_resonances)) then
      do i_l = 1, this % NLS(this % i_urr - 1)
        call this % rm_resonances(i_l) % dealloc()
      end do
      deallocate(this % rm_resonances)
    end if

    ! deallocate pointwise data
    if (allocated(this % urr_E)) deallocate(this % urr_E)
    if (allocated(this % point_xs)) deallocate(this % point_xs)
       
    ! deallocate averaged, infinite-dilute URR cross sections
    if (allocated(this % E_avg_xs)) deallocate(this % E_avg_xs)
    if (allocated(this % avg_xs)) deallocate(this % avg_xs)

    ! deallocate ENDF-6 File 3 cross sections
    call this % dealloc_MF3()

    ! deallocate energy range variables and total angular momenta counts
    call this % dealloc_energy_ranges()
    deallocate(this % NJS)

  end subroutine dealloc_isotope
  

!> Calculate URR cross section values on-the-fly, generating a new realization
!! about each new E_n OR from pre-computed pointwise values reconstructed at
!! simulation initialization
!! @tTODO: fix-up the RRR-URR energy crossover as in xs_otf; probably need
!! to utilize mlbw_resonances, slbw_resonances, in place of rm_resonances, where
!! appropriate
  subroutine new_realization_otf_xs(this, E, T, xs_in, xs_out)

    class(Isotope), intent(inout) :: this ! isotope object
    real(8) :: E ! neutron energy
    real(8) :: T ! isotope temperature [K]
    type(CrossSections), intent(in)  :: xs_in  ! application code cross sections
    type(CrossSections), intent(out) :: xs_out ! library output cross sections

    type(Resonance) :: res ! resonance object
    type(CrossSections) :: avg_xs ! computed average cross sections object
    integer :: i_E      ! first URR parameters energy mesh index
    integer :: iavg     ! average cross section index
    integer :: i_l      ! orbital quantum number index
    integer :: i_J      ! total angular momentum quantum number index
    integer :: i_r      ! resonance index
    integer :: n_res    ! number of resonances to include for a given l-wave
    real(8) :: m            ! pointwise xs energy interpolation factor
    real(8) :: favg         ! average cross section interpolation factor

    if (this % point_urr_xs) then
      i_E = binary_search(this % urr_E, size(this % urr_E), E)
      m = interp_factor(E, this % urr_E(i_E), this % urr_E(i_E + 1), &
           LINEAR_LINEAR)
      xs_out % n = interpolate(&
           m, this % point_xs(i_E) % n, this % point_xs(i_E + 1) % n, LINEAR_LINEAR)
      xs_out % f = interpolate(&
           m, this % point_xs(i_E) % f, this % point_xs(i_E + 1) % f, LINEAR_LINEAR)
      xs_out % g = interpolate(&
           m, this % point_xs(i_E) % g, this % point_xs(i_E + 1) % g, LINEAR_LINEAR)
      if (this % point_xs(i_E) % x > ZERO) then
        xs_out % x = interpolate(&
             m, this % point_xs(i_E) % x, this % point_xs(i_E + 1) % x, LINEAR_LINEAR)
      else
        xs_out % x = ZERO
      end if
      xs_out % t = xs_out % n + xs_out % g + xs_out % f + xs_out % x
      return
    end if

    this % E = E
    this % T = T
    this % k_n = wavenumber(this % AWR, abs(this % E))
    this % k_lam = this % k_n

    ! reset xs accumulators
    call xs_out % flush()

    ! Get resonance parameters for a local realization about E_n
    ! ----------------------------------------------------------
    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, this % NLS(this % i_urr)

      ! set current orbital angular momentum quantum number
      this % L = i_l - 1

      ! penetration
      this % P_l_n = penetration(this % L,&
           this % k_n * this % ac(this % i_urr))
      
      ! resonance energy shift factor
      this % S_l_n = energy_shift(this % L,&
           this % k_n * this % ac(this % i_urr))
          
      ! hard-sphere phase shift
      this % phi_l_n = phase_shift(this % L,&
           this % k_n * this % AP(this % i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = num_contributing_resonances(this % L)

      ! loop over total angular momentum quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)

        ! set current total angular momentum quantum number
        this % J = this % AJ(i_l) % dim1(i_J)

        ! set current partial width degrees of freedom
        this % AMUX = int(this % DOFX(i_l) % dim1(i_J))
        this % AMUN = int(this % DOFN(i_l) % dim1(i_J))
        this % AMUG = int(this % DOFG(i_l) % dim1(i_J))
        this % AMUF = int(this % DOFF(i_l) % dim1(i_J))

        ! zero the resonance counter
        res % i_res = 0

        ! set mean URR parameters to neutron energy
        call this % get_mean_parameters(E, i_l, i_J)

        ! sample resonance parameters for this spin sequence, at this energy
        this % k_lam = this % k_n
        this % P_l_lam = penetration(&
             this % L, this % k_lam * this % ac(this % i_urr)) 
        call this % resonance_parameters(res, i_l, i_J)

        ! loop over the addition of resonances to this ladder
        LOC_RESONANCES_LOOP: do i_r = 1, n_res

          if (parameter_energy_dependence == E_RESONANCE) then
            ! interpolate mean URR parameters to current resonance energy
            call this % get_mean_parameters(res % E_lam, i_l, i_J)
            this % k_lam = wavenumber(this % AWR, abs(res % E_lam))
            this % P_l_lam = penetration(&
                 this % L, this % k_lam * this % ac(this % i_urr))
          end if

          res % i_res = i_r

          ! sample unresolved resonance parameters for this spin
          ! sequence, at this energy
          call this % resonance_parameters(res, i_l, i_J)

        end do LOC_RESONANCES_LOOP
      end do LOC_TOTAL_ANG_MOM_LOOP
    end do LOC_ORBITAL_ANG_MOM_LOOP

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP: do i_l = 1, this % NLS(this % i_urr)

      ! set current orbital angular momentum quantum number
      this % L = i_l - 1

      ! penetration
      this % P_l_n = penetration(&
           this % L, this % k_n * this % ac(this % i_urr))
      
      ! resonance energy shift factor
      this % S_l_n = energy_shift(&
           this % L, this % k_n * this % ac(this % i_urr))
      
      ! hard-sphere phase shift
      this % phi_l_n = phase_shift(&
           this % L, this % k_n * this % AP(this % i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = num_contributing_resonances(this % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)

        ! set current total angular momentum quantum number
        this % J = this % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        this % g_J = (TWO * this % J + ONE) &
             / (FOUR * this % SPI(this % i_urr) + TWO)

        ! loop over resonances localized about E_n
        RESONANCES_LOOP: do i_r = 1, n_res
          res % i_res = i_r
          res % E_lam&
               = this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
          res % Gam_n&
               = this % local_realization(i_l) % J(i_J) % res(i_r) % GN
          res % Gam_g&
               = this % local_realization(i_l) % J(i_J) % res(i_r) % GG
          res % Gam_f&
               = this % local_realization(i_l) % J(i_J) % res(i_r) % GF
          res % Gam_x&
               = this % local_realization(i_l) % J(i_J) % res(i_r) % GX
          res % Gam_t&
               = this % local_realization(i_l) % J(i_J) % res(i_r) % GT

          ! calculate the contribution to the partial cross sections,
          ! at this energy, from an additional resonance
          call this % xs(res)

          ! add this contribution to the accumulated partial cross
          ! section values built up from all resonances
! TODO: move t outside of loop
          call xs_out % accum_resonance(res % xs_contribution)

        end do RESONANCES_LOOP
      end do TOTAL_ANG_MOM_LOOP
    end do ORBITAL_ANG_MOM_LOOP

    ! add potential scattering contribution
    call this % potential_xs(xs_out)

    ! interpret MF3 data according to ENDF-6 LSSF flag
    if (this % LSSF == 0) then
      ! add resonance xs component to MF3 background xs

      call this % add_mf3(xs_out)

    else if (this % LSSF == 1) then
      ! multipy the self-shielding factors by the MF3 background xs

      ! determine index in average xs energy grid
      if (this % E < this % E_avg_xs(1)) then
        iavg = 1

      else if (this % E > this % E_avg_xs(this % num_avg_xs_grid - 1)) then
        iavg = this % num_avg_xs_grid - 1

      else
        iavg = binary_search(this % E_avg_xs, this % num_avg_xs_grid, this % E)

      end if

      favg = interp_factor(&
           E, this % E_avg_xs(iavg), this % E_avg_xs(iavg + 1), this % INT)

      ! interpolate computed averaged URR cross sections
      call this % interpolate_avg_xs(favg, iavg, avg_xs)

      if (avg_xs % x > ZERO&
           .and. this % E <= (ONE + ENDF_PRECISION) * this % E_ex2&
           .and. competitive_structure) then
        ! resonance structure via self-shielding factor if requested
        xs_out % x = xs_out % x / avg_xs % x * xs_in % x
      else
        xs_out % x = xs_in % x
      end if

      ! elastic xs
      xs_out % n = xs_out % n / avg_xs % n * xs_in % n

      ! set negative elastic xs and/or competitive xs to zero
      if (xs_out % n < ZERO) xs_out % n = ZERO
      if (xs_out % x < ZERO) xs_out % x = ZERO

      ! capture xs
      if (avg_xs % g > ZERO) xs_out % g = xs_out % g / avg_xs % g * xs_in % g

      ! fission xs
      if (avg_xs % f > ZERO) xs_out % f = xs_out % f / avg_xs % f * xs_in % f

      ! total xs
      xs_out % t = xs_out % n + xs_out % g + xs_out % f + xs_out % x

    else
      call exit_status(EXIT_FAILURE, 'ENDF-6 LSSF not allowed - must be 0 or 1')
      return

    end if

    this % E_last = this % E

  end subroutine new_realization_otf_xs


!> Calculate a URR cross section from pre-computed probability tables
  subroutine prob_band_xs(this, E, T, xs_in, xs_out, r)

    class(Isotope), intent(inout) :: this ! isotope object
    real(8), intent(in) :: E ! energy
    real(8), intent(in) :: T ! temperature [K]
    type(CrossSections), intent(in)  :: xs_in  ! application code cross sections
    type(CrossSections), intent(out) :: xs_out ! library output cross sections
    real(8), intent(in) :: r ! pseudo-random number

    type(CrossSections) :: avg_xs ! computed average cross sections object
    integer :: i      ! loop index
    integer :: i_E    ! index for energy
    integer :: i_Tlow ! index for lower temperature bound
    integer :: i_Tup  ! index for upper temperature bound
    integer :: i_low  ! band index at lower bounding energy
    integer :: i_up   ! band index at upper bounding energy
    integer :: iavg   ! average cross section index
    real(8) :: fE     ! energy interpolation factor
    real(8) :: fT     ! temperature interpolation factor
    real(8) :: favg   ! average cross section interpolation factor
    real(8) :: xsTlow ! energy-interpolated xs at lower temperature
    real(8) :: xsTup  ! energy-interpolated xs at upper temperature

    call xs_out % flush()

    this % E = E
    this % T = T

    ! determine energy table
    i_E = 1
    do
      if (E < this % E_tabs(i_E + 1)) exit
      i_E = i_E + 1
    end do

    fE = interp_factor(E, this % E_tabs(i_E), this % E_tabs(i_E+1), this % INT)

    if (this % nT_tabs == 1) then
      i_Tlow = 1
      i_Tup  = 1
      fT = ZERO

    else
      if (T < this % T_tabs(1)) then
        call log_message(ERROR, 'Encountered temperature below probability tables.\&
             &  Extrapolating off temperature grid.')
        i_Tlow = 1

      else if (T > this % T_tabs(this % nT_tabs)) then
        call log_message(ERROR, 'Encountered temperature above probability tables.\&
             &  Extrapolating off temperature grid.')
        i_Tlow = this % nT_tabs - 1

      else
        i_Tlow = binary_search(this % T_tabs, this % nT_tabs, T)

      end if
      i_Tup  = i_Tlow + 1
      fT = interp_factor(T, this % T_tabs(i_Tlow), this % T_tabs(i_Tup), temperature_interp_scheme)

    end if

    ! equiprobable band widths in xs-space implies:
    i_low = 1 + floor(r * this % n_bands)

    ! equiprobable band widths across energy implies:
    i_up  = i_low

    ! elastic xs
    xsTlow = interpolate(fE, &
         this % prob_tables(i_E, i_Tlow) % n(i_low) % xs_mean, &
         this % prob_tables(i_E + 1, i_Tlow) % n(i_up) % xs_mean, this % INT)
    xsTup = interpolate(fE, &
         this % prob_tables(i_E, i_Tup) % n(i_low) % xs_mean, &
         this % prob_tables(i_E + 1, i_Tup) % n(i_up) % xs_mean, this % INT)
    xs_out % n = interpolate(fT, xsTlow, xsTup, temperature_interp_scheme)

    ! fission xs
    if ((this % INT == LINEAR_LINEAR .and. (temperature_interp_scheme == LINEAR_LINEAR .or. &
         temperature_interp_scheme == STATISTICAL .or. temperature_interp_scheme == LOW_NEIGHBOR)) .or. &
         (this % prob_tables(i_E, i_Tlow) % f(i_low) % xs_mean > ZERO .and. &
         this % prob_tables(i_E+1, i_Tlow) % f(i_up) % xs_mean > ZERO)) then
      xsTlow = interpolate(fE, &
           this % prob_tables(i_E, i_Tlow) % f(i_low) % xs_mean, &
           this % prob_tables(i_E + 1, i_Tlow) % f(i_up) % xs_mean, this % INT)
      xsTup = interpolate(fE, &
           this % prob_tables(i_E, i_Tup) % f(i_low) % xs_mean, &
           this % prob_tables(i_E + 1, i_Tup) % f(i_up) % xs_mean, this % INT)
      xs_out % f = interpolate(fT, xsTlow, xsTup, temperature_interp_scheme)

    else
      xs_out % f = ZERO

    end if

    ! capture xs
    if ((this % INT == LINEAR_LINEAR .and. (temperature_interp_scheme == LINEAR_LINEAR .or. &
         temperature_interp_scheme == STATISTICAL .or. temperature_interp_scheme == LOW_NEIGHBOR)) .or. &
         (this % prob_tables(i_E, i_Tlow) % g(i_low) % xs_mean > ZERO .and. &
         this % prob_tables(i_E+1, i_Tlow) % g(i_up) % xs_mean > ZERO)) then
      xsTlow = interpolate(fE, &
           this % prob_tables(i_E, i_Tlow) % g(i_low) % xs_mean, &
           this % prob_tables(i_E + 1, i_Tlow) % g(i_up) % xs_mean, this % INT)
      xsTup = interpolate(fE, &
           this % prob_tables(i_E, i_Tup) % g(i_low) % xs_mean, &
           this % prob_tables(i_E + 1, i_Tup) % g(i_up) % xs_mean, this % INT)
      xs_out % g = interpolate(fT, xsTlow, xsTup, temperature_interp_scheme)

    else
      xs_out % g = ZERO

    end if

    ! competitive xs
    if ((this % INT == LINEAR_LINEAR .and. (temperature_interp_scheme == LINEAR_LINEAR .or. &
         temperature_interp_scheme == STATISTICAL .or. temperature_interp_scheme == LOW_NEIGHBOR)) .or. &
         (this % prob_tables(i_E, i_Tlow) % x(i_low) % xs_mean > ZERO .and. &
         this % prob_tables(i_E+1, i_Tlow) % x(i_up) % xs_mean > ZERO)) then
      xsTlow = interpolate(fE, &
           this % prob_tables(i_E, i_Tlow) % x(i_low) % xs_mean, &
           this % prob_tables(i_E + 1, i_Tlow) % x(i_up) % xs_mean, this % INT)
      xsTup = interpolate(fE, &
           this % prob_tables(i_E, i_Tup) % x(i_low) % xs_mean, &
           this % prob_tables(i_E + 1, i_Tup) % x(i_up) % xs_mean, this % INT)
      xs_out % x = interpolate(fT, xsTlow, xsTup, temperature_interp_scheme)

    else
      xs_out % x = ZERO

    end if

    if (this % LSSF == 0) then
! TODO: consider moving addition of File 3 background to prob band calculation
      if (competitive_structure) then
        xs_out % x = xs_in % x + xs_out % x
      else
        xs_out % x = xs_in % x
      end if 
      xs_out % g = xs_in % g + xs_out % g
      xs_out % n = xs_in % n + xs_out % n
      xs_out % f = xs_in % f + xs_out % f

    else if (this % LSSF == 1) then
! TODO: consider moving multiplication by File 3 background to prob band calculation
      iavg = binary_search(this % E_avg_xs, this % num_avg_xs_grid, this % E)
      favg = interp_factor(&
           E, this % E_avg_xs(iavg), this % E_avg_xs(iavg + 1), this % INT)

      ! interpolate computed averaged URR cross sections
      call this % interpolate_avg_xs(favg, iavg, avg_xs)

      ! competitive xs
      if (avg_xs % x > ZERO&
           .and. this % E <= (ONE + ENDF_PRECISION) * this % E_ex2&
           .and. competitive_structure) then
        ! resonance structure via self-shielding factor if requested
        xs_out % x = xs_out % x / avg_xs % x * xs_in % x
      else
        xs_out % x = xs_in % x
      end if

      ! elastic xs
      xs_out % n = xs_out % n / avg_xs % n * xs_in % n

      ! set negative elastic xs and competitive xs to zero
      if (xs_out % n < ZERO) xs_out % n = ZERO
      if (xs_out % x < ZERO) xs_out % x = ZERO

      ! capture xs
      if (avg_xs % g > ZERO) xs_out % g = xs_out % g / avg_xs % g * xs_in % g

      ! fission xs
      if (avg_xs % f > ZERO) xs_out % f = xs_out % f / avg_xs % f * xs_in % f

      ! total xs
      xs_out % t = xs_out % n + xs_out % g + xs_out % f + xs_out % x

    end if

    this % E_last = E

  end subroutine prob_band_xs


!> Calculate a value for the G-function appearing in the NJOY-2012 form
!! of the MLBW resonance formalae
  function G_func(this, E_res, G_n, G_t, i_res) result(G_val)

    class(Isotope), intent(inout) :: this ! isotope pointer
    real(8), intent(in) :: E_res ! energy of resonance contributing to xs at E_n
    real(8), intent(in) :: G_n   ! neutron width of the resonance at E_res
    real(8), intent(in) :: G_t   ! total width of the resonance at E_res
    integer, intent(in) :: i_res ! index of current nearest resonance
    real(8) :: G_val ! G-function value

    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index
    integer :: i_r   ! spin sequence resonance index
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam   ! center-of-mass neutron wavenumber at E_lam
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: P_l_lam ! penetration at resonance energy
    real(8) :: S_l_lam ! resonance energy shift factor

    G_val = ZERO

    i_l = this % L + 1

    TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)
      RESONANCE_LOOP: do i_r = 1, num_contributing_resonances(this % L)

        if (i_r == i_res .and. this % J == this % AJ(i_l) % dim1(i_J)) cycle

        ! if URR parameters have resonance energy dependence
        if (parameter_energy_dependence == E_RESONANCE) then

          ! absolute value of energy in order to handle bound levels which have
          ! negative resonance energies
          k_lam = wavenumber(this % AWR,&
               abs(this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam))
          P_l_lam = penetration(this % L, k_lam * this % ac(this % i_urr))
          S_l_lam = energy_shift(this % L, k_lam * this % ac(this % i_urr))

          E_shift = this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
               + this % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * (S_l_lam - this % S_l_n) / (TWO * P_l_lam)

          Gam_n_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * this % P_l_n / P_l_lam

          if (this % E > (ONE + ENDF_PRECISION) * this % E_ex2) then

            ! two competitive reactions possible;
            ! can't calculate an energy-dependent width because it depends on
            ! the two (unprovided) reaction partial widths
            Gam_x_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GX

          else if (this % E >= this % E_ex1) then

            ! compute an energy-dependent width for the one competitive reaction
            k_n_x = wavenumber(this % AWR, this % E - this % E_ex1)
            k_lam_x = wavenumber(this % AWR,&
                 abs(this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
                 - this % E_ex1))
            Gam_x_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GX&
                 * penetration(this % L, k_n_x   * this % ac(this % i_urr))&
                 / penetration(this % L, k_lam_x * this % ac(this % i_urr))

          else

            Gam_x_n = ZERO

          end if

        else

          ! assume all URR parameters already have neutron energy dependence
          E_shift = this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
          Gam_n_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GN
          Gam_x_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GX

        end if

        Gam_t_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GT&
             - this % local_realization(i_l) % J(i_J) % res(i_r) % GN&
             - this % local_realization(i_l) % J(i_J) % res(i_r) % GX&
             + Gam_n_n&
             + Gam_x_n

        G_val = G_val + G_n * Gam_n_n * (G_t + Gam_t_n)&
             / ((E_res - E_shift) * (E_res - E_shift)&
             + (G_t + Gam_t_n) * (G_t + Gam_t_n) / FOUR)

      end do RESONANCE_LOOP
    end do TOTAL_ANG_MOM_LOOP

  end function G_func


!> Calculate a value for the H-function appearing in the NJOY-2012 form
!! of the MLBW resonance formalae
  function H_func(this, E_res, G_n, G_t, i_res) result(H_val)

    class(Isotope), intent(inout) :: this ! isotope pointer
    real(8), intent(in) :: E_res ! energy of resonance contributing to xs at E_n
    real(8), intent(in) :: G_n   ! neutron width of the resonance at E_res
    real(8), intent(in) :: G_t   ! total width of the resonance at E_res
    integer, intent(in) :: i_res ! index of current nearest resonance
    real(8) :: H_val ! H-function value

    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index
    integer :: i_r   ! spin sequence resonance index
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam   ! center-of-mass neutron wavenumber at E_lam
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: P_l_lam ! penetration at resonance energy
    real(8) :: S_l_lam ! resonance energy shift factor

    H_val = ZERO

    i_l = this % L + 1

    TOTAL_ANG_MOM_LOOP: do i_J = 1, this % NJS(i_l)
      RESONANCE_LOOP: do i_r = 1, num_contributing_resonances(this % L)

        if (i_r == i_res .and. this % J == this % AJ(i_l) % dim1(i_J)) cycle

        ! if URR parameters have resonance energy dependence
        if (parameter_energy_dependence == E_RESONANCE) then

          ! absolute value of energy in order to handle bound levels which have
          ! negative resonance energies
          k_lam = wavenumber(this % AWR,&
               abs(this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam))
          P_l_lam = penetration(this % L, k_lam * this % ac(this % i_urr))
          S_l_lam = energy_shift(this % L, k_lam * this % ac(this % i_urr))

          E_shift = this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
               + this % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * (S_l_lam - this % S_l_n) / (TWO * P_l_lam)

          Gam_n_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * this % P_l_n / P_l_lam

          if (this % E > (ONE + ENDF_PRECISION) * this % E_ex2) then

            ! two competitive reactions possible;
            ! can't calculate an energy-dependent width because it depends on
            ! the two (unprovided) reaction partial widths
            Gam_x_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GX

          else if (this % E >= this % E_ex1) then

            ! compute an energy-dependent width for the one competitive reaction
            k_n_x = wavenumber(this % AWR, this % E - this % E_ex1)
            k_lam_x = wavenumber(this % AWR,&
                 abs(this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
                 - this % E_ex1))
            Gam_x_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GX&
                 * penetration(this % L, k_n_x   * this % ac(this % i_urr))&
                 / penetration(this % L, k_lam_x * this % ac(this % i_urr))

          else

            Gam_x_n = ZERO

          end if

        else

          ! assume all URR parameters already have neutron energy dependence
          E_shift = this % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
          Gam_n_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GN
          Gam_x_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GX

        end if

        Gam_t_n = this % local_realization(i_l) % J(i_J) % res(i_r) % GT&
             - this % local_realization(i_l) % J(i_J) % res(i_r) % GN&
             - this % local_realization(i_l) % J(i_J) % res(i_r) % GX&
             + Gam_n_n&
             + Gam_x_n

        H_val = H_val + G_n * Gam_n_n * (E_res - E_shift)&
             / ((E_res - E_shift) * (E_res - E_shift)&
             + (G_t + Gam_t_n) * (G_t + Gam_t_n) / FOUR)

      end do RESONANCE_LOOP
    end do TOTAL_ANG_MOM_LOOP

  end function H_func


!> Calculate hard sphere penetrability factors
  function penetration(L, rho) result(P)

    integer, intent(in) :: L ! current orbital quantum number
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: P    ! penetration factor

    ! pre-compute exponentiation
    rho2 = rho * rho

    ! calculate penetrability for the current orbital quantum number
    select case(L)

    case(0)
      P = rho

    case(1)
      P = rho * rho2 / (ONE + rho2)

    case(2)
      P = rho * rho2 * rho2 / (9.0_8 + rho2 * (THREE + rho2))

    case(3)
      P = rho * rho2 * rho2 * rho2 &
           / (225.0_8 + rho2 * (45.0_8 + rho2 * (6.0_8 + rho2)))

    case(4)
      P = rho * rho2 * rho2 * rho2 * rho2 &
           / (11025.0_8 + rho2 * (1575.0_8 + rho2 &
           * (135.0_8 + rho2 * (10.0_8 + rho2))))

    case default
      call exit_status(EXIT_FAILURE, 'Orbital quantum number not allowed.')
      return

    end select

  end function penetration


!> Calculate hard sphere phase shifts
  function phase_shift(L, rho) result(phi)

    integer :: L    ! current orbital quantum number
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: rho3 ! rho**3
    real(8) :: rho4 ! rho**4
    real(8) :: phi  ! hard sphere phase shift

    ! pre-compute exponentiations
    rho2 = rho * rho
    rho3 = rho * rho2
    rho4 = rho * rho3

    ! calculate phase shift for the current orbital quantum number
    select case(L)

    case(0)
      phi = rho

    case(1)
      phi = rho - atan(rho)

    case(2)
      phi = rho - atan(THREE * rho / (THREE - rho2))

    case(3)
      phi = rho - atan((15.0_8 * rho - rho3) / (15.0_8 - 6.0_8 * rho2))

    case(4)
      phi = rho - atan((105.0_8 * rho - 10.0_8 * rho3) &
           / (105.0_8 - 45.0_8 * rho2 + rho4))

    case default
      call exit_status(EXIT_FAILURE, 'Orbital quantum number not allowed.')
      return

    end select

  end function phase_shift


!> Calculate resonance energy shift factors
  function energy_shift(L, rho) result(S)

    integer :: L    ! current orbital quantum number
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: S    ! shift factor (for shifting the resonance energy)

    ! pre-compute exponentiation
    rho2 = rho * rho

    ! calculate shift factor for current orbital quantum number
    select case(L)

    case(0)
      S = ZERO

    case(1)
      S = -ONE / (ONE + rho2)

    case(2)
      S = -(18.0_8 + THREE * rho2) / (9.0_8 + rho2 * (THREE + rho2))

    case(3)
      S = -(675.0_8 + rho2 * (90.0_8 + 6.0_8 * rho2)) &
           / (225.0_8 + rho2 * (45.0_8 + rho2 * (6.0_8 + rho2)))

    case(4)
      S = -(44100.0_8 + rho2 * (4725.0_8 + rho2 * (270.0_8 + 10.0_8 * rho2)))&
           / (11025.0_8 + rho2 * (1575.0_8 + rho2 * (135.0_8 &
           + rho2 * (10.0_8 + rho2))))

    case default
      call exit_status(EXIT_FAILURE, 'Orbital quantum number not allowed')
      return

    end select

  end function energy_shift


!> Compute Doppler integral function, psi
  function psi(T, theta, x) result(psi_val)

    real(8)    :: T       ! temperature [K]
    real(8)    :: theta   ! psi argument
    real(8)    :: x       ! psi argument
    real(8)    :: psi_val ! calculated value of psi
    real(8)    :: relerr  ! relative error of the Faddeeva evaluation
    complex(8) :: w_val   ! complex return value of the Faddeeva evaluation

    if (T > ZERO) then
      select case (faddeeva_method)
      case (MIT_W)
        ! S.G. Johnson's Faddeeva evaluation
        relerr = 1.0e-6
        w_val = faddeeva_w(cmplx(theta * x * HALF, theta * HALF, 8), relerr)
        psi_val = SQRT_PI * HALF * theta&
             * real(real(w_val, 8), 8)

      case (QUICK_W)
        ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY)
        psi_val = SQRT_PI * HALF * theta&
             * real(real(quickw(cmplx(theta * x * HALF, theta * HALF, 8)),8),8)

      case default
        call exit_status(EXIT_FAILURE, 'Unrecognized W function evaluation method')
        return

      end select

    else
      psi_val = ONE / (ONE + x*x)

    end if

  end function psi


!> Compute Doppler integral function, chi
  function chi(T, theta, x) result(chi_val)

    real(8)    :: T       ! temperature [K]
    real(8)    :: theta   !
    real(8)    :: x       !
    real(8)    :: chi_val ! calculated value of chi
    complex(8) :: w_val   ! complex return value of the Faddeeva evaluation
    real(8)    :: relerr  ! relative error of the Faddeeva evaluation

    if (T > ZERO) then
      ! evaluate the W (Faddeeva) function
      select case (faddeeva_method)

      case (MIT_W)
        ! S.G. Johnson's Faddeeva evaluation
        relerr = 1.0e-6
        w_val = faddeeva_w(cmplx(theta * x * HALF, theta * HALF, 8), relerr)
        chi_val = SQRT_PI * HALF * theta&
             * real(aimag(w_val), 8)

      case (QUICK_W)
        ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY)
        chi_val = SQRT_PI * HALF * theta&
             * real(aimag(quickw(cmplx(theta * x * HALF, theta * HALF, 8))), 8)

      case default
        call exit_status(EXIT_FAILURE, 'Unrecognized W function evaluation method')
        return

      end select

    else
      chi_val = x / (ONE + x*x)

    end if

  end function chi


!> Compute wavenumber in the center-of-mass reference frame
  function wavenumber(A, E) result(k_val)

    real(8) :: A     ! atomic weight ratio
    real(8) :: E     ! evaluation energy
    real(8) :: k_val ! computed wavenumber

    ! compute center-of-mass neutron wavenumber evaluated at some energy
    k_val = C_1 * A / (A + ONE) * sqrt(E)

  end function wavenumber


!> Add the contribution of potential scattering to the elastic and total cross
!! sections
  subroutine potential_xs(this, xs)

    class(Isotope), intent(in) :: this ! isotope object
    type(CrossSections), intent(inout) :: xs ! partial cross sections object
    integer :: i_l ! orbital quantum number index
    real(8) :: sig_pot ! potential scattering cross section
    real(8) :: phi_l_n ! hard-sphere phase shift

    sig_pot = ZERO
    do i_l = 0, this % NLS(this % i_urr) - 1
      phi_l_n = phase_shift(i_l, this % k_n * this % AP(this % i_urr)) 
      sig_pot = sig_pot&
           + FOUR * PI / (this % k_n * this % k_n) * (TWO * dble(i_l) + ONE)&
           * sin(phi_l_n) * sin(phi_l_n)
    end do

    xs % n = xs % n + sig_pot
    xs % t = xs % t + sig_pot

  end subroutine potential_xs


!> Add the single-history xs realization to the single-batch accumulator
  subroutine accum_history(this, i_E, i_T)

    class(Isotope), intent(inout) :: this ! cross section object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    ! accumulate infinite-dilute xs values for this realization
    this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_t % xs
    this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_n % xs
    this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_g % xs
    this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_f % xs
    this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_x % xs

  end subroutine accum_history


!> Flush the single-history probability table cross section values
  subroutine flush_histories(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
        this % prob_tables(i_E, i_T) % avg_t % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_n % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_g % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_f % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_x % xs = ZERO
      end do
    end do

  end subroutine flush_histories


!> Add the single-batch results to the overall accumulators
  subroutine accum_batch(this, i_E)

    class(Isotope), intent(inout) :: this ! cross section object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    if (this % E /= this % E_tabs(i_E)) return
    do i_T = 1, this % nT_tabs
      ! accumulate single-batch infinite-dilute means
      this % prob_tables(i_E, i_T) % avg_t % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_t % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
           / dble(num_histories_prob_tables)
      this % prob_tables(i_E, i_T) % avg_n % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_n % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
           / dble(num_histories_prob_tables)
      this % prob_tables(i_E, i_T) % avg_g % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_g % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
           / dble(num_histories_prob_tables)
      this % prob_tables(i_E, i_T) % avg_f % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_f % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
           / dble(num_histories_prob_tables)
      this % prob_tables(i_E, i_T) % avg_x % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_x % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
           / dble(num_histories_prob_tables)

      ! accumulate squared single-batch infinite-dilute means
      this % prob_tables(i_E, i_T) % avg_t % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_t % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
           / dble(num_histories_prob_tables))&
           * (this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
           / dble(num_histories_prob_tables))
      this % prob_tables(i_E, i_T) % avg_n % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_n % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
           / dble(num_histories_prob_tables))&
           * (this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
           / dble(num_histories_prob_tables))
      this % prob_tables(i_E, i_T) % avg_g % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_g % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
           / dble(num_histories_prob_tables))&
           * (this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
           / dble(num_histories_prob_tables))
      this % prob_tables(i_E, i_T) % avg_f % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_f % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
           / dble(num_histories_prob_tables))&
           * (this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
           / dble(num_histories_prob_tables))
      this % prob_tables(i_E, i_T) % avg_x % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_x % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
           / dble(num_histories_prob_tables))&
           * (this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
           / dble(num_histories_prob_tables))

    end do

  end subroutine accum_batch


!> Compute batch-based means and standard errors of those means
  subroutine statistics(this, i_bat_int)

    class(Isotope), intent(inout), target :: this ! cross section object
    type(ProbabilityTable), pointer :: ptable ! prob table pointer
    integer :: i_E       ! tabulated URR energy index
    integer :: i_T       ! temperature index
    integer :: i_bat_int ! current batch index
    real(8) :: i_bat ! current batch index for calcs

    i_bat = dble(i_bat_int)

    do i_E = 1, this % nE_tabs
      if (this % E /= this % E_tabs(i_E)) cycle
      do i_T = 1, this % nT_tabs
        ptable => this % prob_tables(i_E, i_T)
        ptable % avg_t % xs_mean = ptable % avg_t % xs_sum / i_bat
        ptable % avg_n % xs_mean = ptable % avg_n % xs_sum / i_bat
        ptable % avg_g % xs_mean = ptable % avg_g % xs_sum / i_bat
        ptable % avg_f % xs_mean = ptable % avg_f % xs_sum / i_bat
        ptable % avg_x % xs_mean = ptable % avg_x % xs_sum / i_bat

        if (ptable % avg_t % xs_mean > ZERO .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_t % xs_sem&
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_t % xs_sum2 / i_bat&
               - (ptable % avg_t % xs_mean) * (ptable % avg_t % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_t % rel_unc &
               = ptable % avg_t % xs_sem / ptable % avg_t % xs_mean
        end if
        if (ptable % avg_n % xs_mean > ZERO .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_n % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_n % xs_sum2 / i_bat&
               - (ptable % avg_n % xs_mean) * (ptable % avg_n % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_n % rel_unc &
               = ptable % avg_n % xs_sem / ptable % avg_n % xs_mean
        end if
        if (ptable % avg_g % xs_mean > ZERO .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_g % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_g % xs_sum2 / i_bat&
               - (ptable % avg_g % xs_mean) * (ptable % avg_g % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_g % rel_unc &
               = ptable % avg_g % xs_sem / ptable % avg_g % xs_mean
        end if
        if (ptable % avg_f % xs_mean > ZERO .and. i_bat_int > 1&
             .and. (ptable % avg_f % xs_sum2 / i_bat &
             - ptable % avg_f % xs_mean * ptable % avg_f % xs_mean) > ZERO) then
          ! compute standard errors of mean xs values
          ptable % avg_f % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_f % xs_sum2 / i_bat&
               - (ptable % avg_f % xs_mean) * (ptable % avg_f % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_f % rel_unc &
               = ptable % avg_f % xs_sem / ptable % avg_f % xs_mean
        end if
        if (ptable % avg_x % xs_mean > ZERO &
             .and. i_bat_int > 1 .and. competitive_structure) then
          ! compute standard errors of mean xs values
          ptable % avg_x % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_x % xs_sum2 / i_bat&
               - (ptable % avg_x % xs_mean) * (ptable % avg_x % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_x % rel_unc &
               = ptable % avg_x % xs_sem / ptable % avg_x % xs_mean
        end if
      end do
    end do

    nullify(ptable)

  end subroutine statistics


!> Determine the energy of the highest-energy resolved resonance
!! region resonance for a given (l,J) spin sequence
  function last_resolved_resonance_energy(this, l_val, J_val) result(E_val)

    class(Isotope), intent(in) :: this ! isotope object
    integer, intent(in) :: l_val ! orbital quantum number
    real(8), intent(in) :: J_val ! total angular momentum quantum number
    real(8) :: E_val ! highest-energy RRR resonance energy for (l,J)

    integer :: iso   ! isotope index
    integer :: i_res ! RRR resonance index for a given l

    select case (this % LRF(this % i_urr - 1))
    case (SLBW)
      do i_res = size(this % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (this % bw_resonances(l_val + 1) % res(i_res) % AJ== J_val&
             .and. this % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < this % EL(this % i_urr)) then
          E_val = this % bw_resonances(l_val + 1) % res(i_res) % E_lam
          exit
        end if
      end do
    
    case (MLBW)
      do i_res = size(this % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (this % bw_resonances(l_val + 1) % res(i_res) % AJ == J_val&
             .and. this % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < this % EL(this % i_urr)) then
          E_val = this % bw_resonances(l_val + 1) % res(i_res) % E_lam
          exit
        end if
      end do
    
    case (REICH_MOORE)
      do i_res = size(this % rm_resonances(l_val + 1) % res(:)), 1, -1
        if (this % rm_resonances(l_val + 1) % res(i_res) % AJ == J_val&
             .and. this % rm_resonances(l_val + 1) % res(i_res) % E_lam&
             < this % EL(this % i_urr)) then
          E_val = this % rm_resonances(l_val + 1) % res(i_res) % E_lam
          exit
        end if
      end do
    
    case default
      call exit_status(EXIT_FAILURE, 'Unrecognized/unsupported RRR formalism')
      return

    end select

  end function last_resolved_resonance_energy


!> Find the index of the RRR resonance which we need to add the
!! contribution of to a URR xs
  function resolved_resonance_index(this, l_val, J_val, n_rrr_res) result(i_res)

    class(Isotope), intent(in) :: this ! isotope object
    integer, intent(in) :: l_val     ! orbital quantum number
    real(8), intent(in) :: J_val     ! total angular momentum quantum number
    integer, intent(in) :: n_rrr_res ! how many RRR resonances to go back
    integer :: i_res ! index of the RRR resonance cnt_res resonances back

    integer :: cnt_res ! how many RRR resonances have we gone back

    cnt_res   = 0
    select case (this % LRF(this % i_urr - 1))
    case (SLBW)
      do i_res = size(this % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (this % bw_resonances(l_val + 1) % res(i_res) % AJ == J_val &
             .and. this % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < this % EL(this % i_urr)) cnt_res = cnt_res + 1
        if (cnt_res == n_rrr_res) exit
      end do

    case (MLBW)
      do i_res = size(this % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (this % bw_resonances(l_val + 1) % res(i_res) % AJ == J_val &
             .and. this % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < this % EL(this % i_urr)) cnt_res = cnt_res + 1
        if (cnt_res == n_rrr_res) exit
      end do
    
    case (REICH_MOORE)
      do i_res = size(this % rm_resonances(l_val + 1) % res(:)), 1, -1
        if (this % rm_resonances(l_val + 1) % res(i_res) % AJ == J_val &
             .and. this % rm_resonances(l_val + 1) % res(i_res) % E_lam&
             < this % EL(this % i_urr)) cnt_res = cnt_res + 1
        if (cnt_res == n_rrr_res) exit
      end do
    
    case default
      call exit_status(EXIT_FAILURE, 'Unrecognized/unsupported RRR formalism')
      return

    end select

    ! if there aren't enough contributing RRR resonances
    if (cnt_res < n_rrr_res) i_res = 0

  end function resolved_resonance_index


!> Add the contribution from an additional resonance
  subroutine resonance_contribution(this, res, i_E, i_T)

    class(Isotope),  intent(inout) :: this ! isotope object
    type(Resonance), intent(in)    :: res  ! resonance object
    integer, intent(in) :: i_E ! tabulated URR energy index
    integer, intent(in) :: i_T ! temperature index

    this % prob_tables(i_E, i_T) % avg_n % xs&
         = this % prob_tables(i_E, i_T) % avg_n % xs + res % xs_contribution % n
    this % prob_tables(i_E, i_T) % avg_g % xs&
         = this % prob_tables(i_E, i_T) % avg_g % xs + res % xs_contribution % g
    this % prob_tables(i_E, i_T) % avg_f % xs&
         = this % prob_tables(i_E, i_T) % avg_f % xs + res % xs_contribution % f
    this % prob_tables(i_E, i_T) % avg_x % xs&
         = this % prob_tables(i_E, i_T) % avg_x % xs + res % xs_contribution % x
    this % prob_tables(i_E, i_T) % avg_t % xs&
         = this % prob_tables(i_E, i_T) % avg_t % xs + res % xs_contribution % t

  end subroutine resonance_contribution


!> Zero out probability table batch accumulators
  subroutine flush_batches(this)

    class(Isotope), intent(inout) :: this ! isotope object

    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index
    integer :: i_b ! probability band index

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
        do i_b = 1, this % n_bands
          this % prob_tables(i_E, i_T) % t(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % n(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % g(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % f(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % x(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % n(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % g(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % f(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % x(i_b) % cnt_tmp_sum = 0
        end do
        this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_t % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_n % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_g % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_f % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_x % cnt_tmp_sum = 0
      end do
    end do

  end subroutine flush_batches


!> Add URR resonance parameters for a single URR resonance to the realization
  subroutine set_parameters(this, res, i_ens, i_l, i_J)

    class(Isotope),  intent(inout) :: this ! isotope object
    type(Resonance), intent(inout) :: res  ! resonance object
    integer, intent(in) :: i_ens ! resonance ensemble index
    integer, intent(in) :: i_l   ! orbital angular momentum index
    integer, intent(in) :: i_J   ! total angular momentum quantum number index

    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % E_lam = res % E_lam
    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % AJ = this % J
    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GN = res % Gam_n
    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GG = res % Gam_g
    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GF = res % Gam_f
    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GX = res % Gam_x
    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GT = res % Gam_t

    ! point to next resonance
    allocate(this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % next)
    this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
         => this % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % next

  end subroutine set_parameters


!> Get the URR resonance parameters for a single resonance
  subroutine get_parameters(this, res, i_res, i_l, i_J, i_ER)

    class(Isotope),  intent(inout) :: this ! isotope object
    type(Resonance), intent(inout) :: res  ! resonance object
    integer, intent(inout) :: i_res ! resonance counter
    integer, intent(in) :: i_l   ! orbital quantum number index
    integer, intent(in) :: i_J   ! total angular momentum quantum number
    integer, intent(in) :: i_ER  ! resonance energy region index

    select case (this % LRU(i_ER))

    ! resolved parameters
    case (1)
      select case (this % LRF(i_ER))
      case (SLBW)
        if (i_res > size(this % bw_resonances(i_l) % res(:))) then
          i_res = size(this % bw_resonances(i_l) % res(:))
          res % E_lam = res % E_lam&
               + (this % bw_resonances(i_l) % res(i_res) % E_lam&
               -  this % bw_resonances(i_l) % res(i_res - 1) % E_lam)
        else
          res % E_lam = this % bw_resonances(i_l) % res(i_res) % E_lam
        end if

        res % Gam_n = this % bw_resonances(i_l) % res(i_res) % GN
        res % Gam_g = this % bw_resonances(i_l) % res(i_res) % GG
        res % Gam_f = this % bw_resonances(i_l) % res(i_res) % GF
        res % Gam_t = this % bw_resonances(i_l) % res(i_res) % GT
        res % Gam_x = res % Gam_t - res % Gam_n - res % Gam_g - res % Gam_f

      case (MLBW)
        if (i_res > size(this % bw_resonances(i_l) % res(:))) then
          i_res = size(this % bw_resonances(i_l) % res(:))
          res % E_lam = res % E_lam&
               + (this % bw_resonances(i_l) % res(i_res) % E_lam&
               -  this % bw_resonances(i_l) % res(i_res - 1) % E_lam)
        else
          res % E_lam = this % bw_resonances(i_l) % res(i_res) % E_lam
        end if
        
        res % Gam_n = this % bw_resonances(i_l) % res(i_res) % GN
        res % Gam_g = this % bw_resonances(i_l) % res(i_res) % GG
        res % Gam_f = this % bw_resonances(i_l) % res(i_res) % GF
        res % Gam_t = this % bw_resonances(i_l) % res(i_res) % GT
        res % Gam_x = res % Gam_t - res % Gam_n - res % Gam_g - res % Gam_f

      case (REICH_MOORE)
        if (i_res > size(this % rm_resonances(i_l) % res(:))) then
          i_res = size(this % rm_resonances(i_l) % res(:))
          res % E_lam = res % E_lam&
               + (this % rm_resonances(i_l) % res(i_res) % E_lam&
               -  this % rm_resonances(i_l) % res(i_res - 1) % E_lam)
        else
          res % E_lam = this % rm_resonances(i_l) % res(i_res) % E_lam
        end if

        res % Gam_n = this % rm_resonances(i_l) % res(i_res) % GN
        res % Gam_g = this % rm_resonances(i_l) % res(i_res) % GG
        res % Gam_f = this % rm_resonances(i_l) % res(i_res) % GFA&
             + this % rm_resonances(i_l) % res(i_res) % GFB
        res % Gam_x = ZERO
        res % Gam_t = res % Gam_n + res % Gam_g + res % Gam_f + res % Gam_x

      case default
        call exit_status(EXIT_FAILURE, 'Unrecognized resolved resonance region formalism')
        return

      end select

    ! unresolved parameters
    case (2)
      res % E_lam&
           = this % urr_resonances(i_l, i_realization) % J(i_J) % res(i_res) % E_lam
      res % Gam_n&
           = this % urr_resonances(i_l, i_realization) % J(i_J) % res(i_res) % GN
      res % Gam_g&
           = this % urr_resonances(i_l, i_realization) % J(i_J) % res(i_res) % GG
      res % Gam_f&
           = this % urr_resonances(i_l, i_realization) % J(i_J) % res(i_res) % GF
      res % Gam_x&
           = this % urr_resonances(i_l, i_realization) % J(i_J) % res(i_res) % GX
      res % Gam_t&
           = this % urr_resonances(i_l, i_realization) % J(i_J) % res(i_res) % GT

    case default
      call exit_status(EXIT_FAILURE, 'Only 1 and 2 are supported ENDF-6 LRU values')
      return

    end select

    this % local_realization(i_l) % J(i_J) % res(res % i_res) % E_lam&
         = res % E_lam
    this % local_realization(i_l) % J(i_J) % res(res % i_res) % GN = res % Gam_n
    this % local_realization(i_l) % J(i_J) % res(res % i_res) % GG = res % Gam_g
    this % local_realization(i_l) % J(i_J) % res(res % i_res) % GF = res % Gam_f
    this % local_realization(i_l) % J(i_J) % res(res % i_res) % GX = res % Gam_x
    this % local_realization(i_l) % J(i_J) % res(res % i_res) % GT = res % Gam_t

  end subroutine get_parameters


!> Get the mean resonance parameters at an energy
  subroutine get_mean_parameters(this, E_res, i_l, i_J)

    class(Isotope), intent(inout) :: this ! isotope object
    real(8), intent(in) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    integer, intent(in) :: i_l   ! orbital quantum number index
    integer, intent(in) :: i_J   ! total angular momentum quantum number

    integer :: i_E ! tabulated URR parameters energy index
    real(8) :: fendf ! ENDF6 URR parameters energy interpolation factor

    ! compute interpolation factor
    if (E_res < this % ES(1)) then
      i_E = 1
    else if (E_res > this % ES(this % NE - 1)) then
      i_E = this % NE - 1
    else
      i_E = binary_search(this % ES, this % NE, E_res)
    end if
    fendf = interp_factor(E_res, this % ES(i_E), this % ES(i_E + 1),&
      this % INT)

    ! set current mean unresolved resonance parameters
    this % D   = interpolate(fendf, &
         this % D_mean(i_l) % dim2(i_J) % dim1(i_E), &
         this % D_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
    this % GN0 = interpolate(fendf, &
         this % GN0_mean(i_l) % dim2(i_J) % dim1(i_E), &
         this % GN0_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
    this % GG  = interpolate(fendf, &
         this % GG_mean(i_l) % dim2(i_J) % dim1(i_E), &
         this % GG_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
    if (this % INT == LINEAR_LINEAR .or. &
         this % GF_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
      this % GF  = interpolate(fendf, &
           this % GF_mean(i_l) % dim2(i_J) % dim1(i_E), &
           this % GF_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
    else
      this % GF = ZERO
    end if
    if (this % INT == LINEAR_LINEAR .or. &
         this % GX_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
      this % GX  = interpolate(fendf, &
           this % GX_mean(i_l) % dim2(i_J) % dim1(i_E), &
           this % GX_mean(i_l) % dim2(i_J) % dim1(i_E + 1), this % INT)
    else
      this % GX = ZERO
    end if

  end subroutine get_mean_parameters


!> Interpolate the averaged, infinite-dilute URR cross sections computed via
!! Monte Carlo from mean resonance parameters (i.e., not the evaluator-supplied
!! File 3 background cross sections)
  subroutine interpolate_avg_xs(this, f, i_avg, avg_xs)

    class(Isotope), intent(in) :: this ! isotope object
    real(8), intent(in) :: f ! interpolation factor
    integer, intent(in) :: i_avg ! averaged cross section index
    type(CrossSections), intent(out) :: avg_xs ! interpolated cross sections

    ! elastic xs
    avg_xs % n = interpolate(&
         f, this % avg_xs(i_avg) % n, this % avg_xs(i_avg + 1) % n, this % INT)

    ! fission xs
    if (this % INT == LINEAR_LINEAR&
         .or. (this % avg_xs(i_avg) % f > ZERO&
         .and. this % avg_xs(i_avg + 1) % f > ZERO)) then
      avg_xs % f = interpolate(&
           f, this % avg_xs(i_avg) % f, this % avg_xs(i_avg + 1) % f, this % INT)

    else
      avg_xs % f = ZERO

    end if

    ! capture xs
    if (this % INT == LINEAR_LINEAR&
         .or. (this % avg_xs(i_avg) % g > ZERO&
         .and. this % avg_xs(i_avg + 1) % g > ZERO)) then
      avg_xs % g = interpolate(&
           f, this % avg_xs(i_avg) % g, this % avg_xs(i_avg + 1) % g, this % INT)

    else
      avg_xs % g = ZERO

    end if

    ! competitive xs
    if (this % INT == LINEAR_LINEAR&
         .or. (this % avg_xs(i_avg) % x > ZERO&
         .and. this % avg_xs(i_avg + 1) % x > ZERO)) then
      avg_xs % x = interpolate(&
           f, this % avg_xs(i_avg) % x, this % avg_xs(i_avg + 1) % x, this % INT)

    else
      avg_xs % x = ZERO

    end if

    avg_xs % t = avg_xs % n + avg_xs % g + avg_xs % f + avg_xs % x

  end subroutine interpolate_avg_xs


!> Interpolate the evaluator-supplied ENDF-6 File 3 background and either add a
!! resonance component (LSSF = 0) or interpolate the computed average URR xs
!! value and then multiply the File 3 xs by a resonance component divided by
!! the average, which is a self-shielding factor (LSSF = 1)
  subroutine self_shielding(this, xs)

    class(Isotope), intent(inout) :: this ! isotope object
    type(CrossSections), intent(inout) :: xs ! partial cross sections object

    type(CrossSections) :: avg_xs ! computed average cross sections object
    integer :: i_mf3 ! index in a File 3 grid
    integer :: i_avg ! index in the averaged xs grid
    real(8) :: mf3_n ! MF3 elastic xs
    real(8) :: mf3_g ! MF3 capture xs
    real(8) :: mf3_f ! MF3 fission xs
    real(8) :: mf3_x ! MF3 competitive xs
    real(8) :: f_mf3 ! interpolation factor in a File 3 grid
    real(8) :: f_avg ! interpolation factor in the averaged xs grid

    ! elastic xs
    if (this % E < this % MF3_n_e(1)) then
      call exit_status(EXIT_FAILURE, 'Energy is below File 3 elastic energy grid')
      return
    else if (this % E > this % MF3_n_e(size(this % MF3_n_e))) then
      call exit_status(EXIT_FAILURE, 'Energy is above File 3 elastic energy grid')
      return
    else
      i_mf3 = binary_search(this % MF3_n_e, size(this % MF3_n_e), this % E)
      if (this % INT == LINEAR_LINEAR&
           .or. (this % MF3_n(i_mf3) > ZERO&
           .and. this % MF3_n(i_mf3 + 1) > ZERO)) then
        f_mf3 = interp_factor(this % E, this % MF3_n_e(i_mf3),&
             this % MF3_n_e(i_mf3 + 1), this % INT)
        mf3_n = interpolate(f_mf3, this % MF3_n(i_mf3),&
             this % MF3_n(i_mf3 + 1), this % INT)
      else
        mf3_n = ZERO
      end if
    end if
    if (mf3_n < ZERO) mf3_n = ZERO

    ! competitive reaction xs
    if (this % E < this % MF3_x_e(1)) then
      mf3_x = ZERO
    else if (this % E > this % MF3_x_e(size(this % MF3_x_e))) then
      call exit_status(EXIT_FAILURE, 'Energy is above File 3 competitive energy grid')
      return
    else
      i_mf3 = binary_search(this % MF3_x_e, size(this % MF3_x_e), this % E)
      if (this % INT == LINEAR_LINEAR&
           .or. (this % MF3_x(i_mf3) > ZERO&
           .and. this % MF3_x(i_mf3 + 1) > ZERO)) then
        f_mf3 = interp_factor(this % E, this % MF3_x_e(i_mf3),&
             this % MF3_x_e(i_mf3 + 1), this % INT)
        mf3_x = interpolate(f_mf3, this % MF3_x(i_mf3),&
             this % MF3_x(i_mf3 + 1), this % INT)
      else
        mf3_x = ZERO
      end if
    end if
    if (mf3_x < ZERO) mf3_x = ZERO

    ! capture xs
    if (this % E < this % MF3_g_e(1)) then
      call exit_status(EXIT_FAILURE, 'Energy is below File 3 capture energy grid')
      return
    else if (this % E > this % MF3_g_e(size(this % MF3_g_e))) then
      call exit_status(EXIT_FAILURE, 'Energy is above File 3 capture energy grid')
      return
    else
      i_mf3 = binary_search(this % MF3_g_e, size(this % MF3_g_e), this % E)
      if (this % INT == LINEAR_LINEAR&
           .or. (this % MF3_g(i_mf3) > ZERO&
           .and. this % MF3_g(i_mf3 + 1) > ZERO)) then
        f_mf3 = interp_factor(this % E, this % MF3_g_e(i_mf3),&
             this % MF3_g_e(i_mf3 + 1), this % INT)
        mf3_g = interpolate(f_mf3, this % MF3_g(i_mf3),&
             this % MF3_g(i_mf3 + 1), this % INT)
      else
        mf3_g = ZERO
      end if
    end if
    if (mf3_g < ZERO) mf3_g = ZERO

    ! fission xs
    if (.not. (allocated(this % MF3_f_e))) then
      mf3_f = ZERO
    else
      if (this % E < this % MF3_f_e(1)) then
        mf3_f = ZERO
      else if (this % E > this % MF3_f_e(size(this % MF3_f_e))) then
        call exit_status(EXIT_FAILURE, 'Energy is above File 3 fission energy grid')
        return
      else
        i_mf3 = binary_search(this % MF3_f_e, size(this % MF3_f_e), this % E)
        if (this % INT == LINEAR_LINEAR&
             .or. (this % MF3_f(i_mf3) > ZERO&
             .and. this % MF3_f(i_mf3 + 1) > ZERO)) then
          f_mf3 = interp_factor(this % E, this % MF3_f_e(i_mf3),&
               this % MF3_f_e(i_mf3 + 1), this % INT)
          mf3_f = interpolate(f_mf3, this % MF3_f(i_mf3),&
               this % MF3_f(i_mf3 + 1), this % INT)
        else
          mf3_f = ZERO
        end if
      end if
    end if
    if (mf3_f < ZERO) mf3_f = ZERO

    ! add resonance component to File 3 xs or use it as a self-shielding factor
    if (this % LSSF == 0) then
      xs % n = mf3_n + xs % n
      xs % g = mf3_g + xs % g
      xs % f = mf3_f + xs % f
      if (competitive_structure) xs % x = mf3_x + xs % x

    else if (this % LSSF == 1) then

      if (this % E < this % E_avg_xs(1)) then
        i_avg = 1
      else if (this % E > this % E_avg_xs(this % num_avg_xs_grid - 1)) then
        i_avg = this % num_avg_xs_grid - 1
      else
        i_avg = binary_search(this % E_avg_xs, this % num_avg_xs_grid, this % E)
      end if

      f_avg = interp_factor(this % E,&
           this % E_avg_xs(i_avg), this % E_avg_xs(i_avg + 1), this % INT)

      call this % interpolate_avg_xs(f_avg, i_avg, avg_xs)

      ! competitive xs
      if (avg_xs % x > ZERO&
           .and. this % E <= (ONE + ENDF_PRECISION) * this % E_ex2&
           .and. competitive_structure)&
           xs % x = mf3_x * xs % x / avg_xs % x

      ! elastic scattering xs
      xs % n = mf3_n * xs % n / avg_xs % n

      ! set negative elastic xs and competitive xs to zero
      if (xs % n < ZERO) xs % n = ZERO
      if (xs % x < ZERO) xs % x = ZERO

      ! radiative capture xs
      if (avg_xs % g > ZERO) xs % g = mf3_g * xs % g / avg_xs % g

      ! fission xs
      if (avg_xs % f > ZERO) xs % f = mf3_f * xs % f / avg_xs % f

    else
      call exit_status(EXIT_FAILURE, 'ENDF-6 LSSF not allowed - must be 0 or 1.')
      return

    end if

    xs % t = xs % n + xs % g + xs % f + xs % x

  end subroutine self_shielding


!> Add the resonance xs component to the evaluator-supplied background xs in MF3
  subroutine add_mf3(this, xs)

    class(Isotope), intent(in) :: this ! isotope object
    type(CrossSections), intent(inout)  :: xs ! partial cross sections object

    type(CrossSections) :: xs_tmp ! local partial cross sections object
    integer :: i_nuc  ! nuclide index
    integer :: i_grid ! background energy grid index
    real(8) :: fmf3   ! File 3 interpolation factor

    call xs_tmp % flush()

    ! elastic xs
    if (this % E < this % MF3_n_e(1)) then
      call exit_status(EXIT_FAILURE, 'Energy is below File 3 elastic energy grid')
      return

    else if (this % E > this % MF3_n_e(size(this % MF3_n_e))) then
      call exit_status(EXIT_FAILURE, 'Energy is above File 3 elastic energy grid')
      return

    else
      i_grid = binary_search(this % MF3_n_e, size(this % MF3_n_e),this % E)

      if (this % INT == LINEAR_LINEAR &
           .or. (this % MF3_n(i_grid) > ZERO &
           .and. this % MF3_n(i_grid + 1) > ZERO)) then
        fmf3 = interp_factor(this % E, this % MF3_n_e(i_grid), &
             this % MF3_n_e(i_grid + 1), this % INT)
        xs % n = interpolate(fmf3, this % MF3_n(i_grid),&
             this % MF3_n(i_grid + 1), this % INT) + xs % n

      else
        xs % n = xs % n

      end if
    end if

    if (xs % n < ZERO) xs % n = ZERO

    ! competitive reaction xs
    if (this % E < this % MF3_x_e(1)) then
      xs % x = ZERO

    else if (this % E > this % MF3_x_e(size(this % MF3_x_e))) then
      call exit_status(EXIT_FAILURE, 'Energy is above File 3 inelastic energy grid')
      return

    else
      i_grid = binary_search(this % MF3_x_e, size(this % MF3_x_e), this % E)

      if (this % INT == LINEAR_LINEAR &
           .or. (this % MF3_x(i_grid) > ZERO &
           .and. this % MF3_x(i_grid + 1) > ZERO)) then
        fmf3 = interp_factor(this % E, this % MF3_x_e(i_grid), &
             this % MF3_x_e(i_grid + 1), this % INT)
        xs_tmp % x = interpolate(fmf3, this % MF3_x(i_grid), &
             this % MF3_x(i_grid + 1), this % INT)

      else
        xs_tmp % x = ZERO

      end if

      if (competitive_structure) xs % x = xs_tmp % x + xs % x

    end if

    if (xs % x < ZERO) xs % x = ZERO

    ! capture xs
    if (this % E < this % MF3_g_e(1)) then
      call exit_status(EXIT_FAILURE, 'Energy is below File 3 capture energy grid')
      return

    else if (this % E > this % MF3_g_e(size(this % MF3_g_e))) then
      call exit_status(EXIT_FAILURE, 'Energy is above File 3 capture energy grid')
      return

    else
      i_grid = binary_search(this % MF3_g_e, size(this % MF3_g_e), this % E)

      if (this % INT == LINEAR_LINEAR &
           .or. (this % MF3_g(i_grid) > ZERO &
           .and. this % MF3_g(i_grid + 1) > ZERO)) then
        fmf3 = interp_factor(this % E, this % MF3_g_e(i_grid), &
             this % MF3_g_e(i_grid + 1), this % INT)
        xs % g = interpolate(fmf3, this % MF3_g(i_grid), &
             this % MF3_g(i_grid + 1), this % INT) + xs % g

      else
        xs % g = xs % g

      end if
    end if

    ! fission xs
    if (.not. (allocated(this % MF3_f_e))) then
      xs % f = xs % f

    else
      if (this % E < this % MF3_f_e(1)) then
        xs % f = xs % f

      else if (this % E > this % MF3_f_e(size(this % MF3_f_e))) then
        call exit_status(EXIT_FAILURE, 'Energy is above File 3 fission energy grid')
        return

      else
        i_grid = binary_search(this % MF3_f_e, size(this % MF3_f_e), this%E)

        if (this % INT == LINEAR_LINEAR &
             .or. (this % MF3_f(i_grid) > ZERO &
             .and. this % MF3_f(i_grid + 1) > ZERO)) then
          fmf3 = interp_factor(this % E, this % MF3_f_e(i_grid), &
               this % MF3_f_e(i_grid + 1), this % INT)
          xs % f = interpolate(fmf3, this%MF3_f(i_grid),&
               this % MF3_f(i_grid + 1), this % INT) + xs % f

        else
          xs % f = xs % f

        end if
      end if
    end if

    xs % t = xs % n + xs % g + xs % f + xs % x

  end subroutine add_mf3


!> Deallocate ENDF-6 File 3 evaluator-supplied background cross sections
  subroutine dealloc_mf3(this)

    class(Isotope), intent(inout) :: this ! isotope object

    if (allocated(this % MF3_n_e)) deallocate(this % MF3_n_e)
    if (allocated(this % MF3_f_e)) deallocate(this % MF3_f_e)
    if (allocated(this % MF3_g_e)) deallocate(this % MF3_g_e)
    if (allocated(this % MF3_x_e)) deallocate(this % MF3_x_e)
    if (allocated(this % MF3_n)) deallocate(this % MF3_n)
    if (allocated(this % MF3_f)) deallocate(this % MF3_f)
    if (allocated(this % MF3_g)) deallocate(this % MF3_g)
    if (allocated(this % MF3_x)) deallocate(this % MF3_x)

  end subroutine dealloc_mf3


!> Sample the energy spacing between adjacent resonances
  subroutine level_spacing(this, res, i_l, i_J)

    class(Isotope), intent(inout) :: this ! isotope object
    type(Resonance), intent(inout) :: res ! resonance object
    integer :: n_res ! number of resonances to include for a given l-wave
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index
    real(8) :: D_lJ ! sampled nuclear level spacing

    if (this % E == this % E_last) then
      ! Energy hasn't changed since last realization, so use the same one
      if (res % i_res == 0) then
        res % E_lam = this % local_realization(i_l) % J(i_J) % res(1) % E_lam
      else
        res % E_lam&
             = this % local_realization(i_l) % J(i_J) % res(res%i_res) % E_lam
      end if

    else
      ! sample a level spacing from the Wigner distribution
      D_lJ = wigner_level_spacing(this % D, prn())

      if (res % i_res == 0) then
        ! set lowest-energy resonance for this ladder well below the energy grid
        ! point such that the ladder spans a sufficient energy range
        n_res = num_contributing_resonances(this % L)
        res % E_lam = (this % E - n_res/2 * this % D) &
             + (ONE - TWO * prn()) * D_lJ
      
      else
        ! add subsequent resonance energies at the sampled spacing above the
        ! last resonance
        res % E_lam = res % E_lam + D_lJ
        this % local_realization(i_l) % J(i_J) % res(res % i_res) % E_lam&
             = res % E_lam

      end if
    end if

  end subroutine level_spacing


!> Sample the channel partial widths at E_lambda when generating
!! a full resonance ensemble or at E_n when generating localized parameters for
!! an on-the-fly cross section calculation
  subroutine channel_widths(this, res, i_l, i_J)

    class(Isotope), intent(inout) :: this ! isotope object
    class(Resonance), intent(inout) :: res ! resonance object
    integer :: i_l    ! orbital angular momentum quantum number index
    integer :: i_J    ! total angular momentum quantum number index
    integer :: i_tabn ! elastic chi-squared table index
    integer :: i_tabg ! capture chi-squared table index
    integer :: i_tabf ! fission chi-squared table index
    integer :: i_tabx ! competitivechi-squared table index
    real(8) :: rho    ! derived variable
    real(8) :: nu     ! derived variable

    if (this % E == this % E_last) then
      if (res % i_res == 0) return
      ! Energy hasn't changed since last realization, so use the same one
      res % Gam_n&
           = this % local_realization(i_l) % J(i_J) % res(res % i_res) % GN
      res % Gam_f&
           = this % local_realization(i_l) % J(i_J) % res(res % i_res) % GF
      res % Gam_g&
           = this % local_realization(i_l) % J(i_J) % res(res % i_res) % GG
      res % Gam_x&
           = this % local_realization(i_l) % J(i_J) % res(res % i_res) % GX
      res % Gam_t&
           = this % local_realization(i_l) % J(i_J) % res(res % i_res) % GT

    else
      ! indices into table of equiprobable chi-squared function values
      i_tabn = 1 + floor(prn() * 100.0_8)
      i_tabf = 1 + floor(prn() * 100.0_8)
      i_tabg = 1 + floor(prn() * 100.0_8)
      i_tabx = 1 + floor(prn() * 100.0_8)

      ! calculate widths from sampled values
      ! neutron width
      if (this % AMUN > 0) then
        ! compute factors needed to go from the mean reduced width that is
        ! provided by ENDF for elastic scattering to a partial width
        ! (use absolute value energies when handling  bound levels which have
        ! negative resonance energies - this is ENDF-6 convention, not theory)
        if (parameter_energy_dependence == E_NEUTRON) then
          rho = this % k_n * this % ac(this % i_urr)
          nu  = this % P_l_n / rho
          res % Gam_n = this % GN0 * sqrt(abs(this % E)) * nu &
               * chi2(i_tabn, this % AMUN)
        else if (parameter_energy_dependence == E_RESONANCE) then
          rho = this % k_lam * this % ac(this % i_urr)
          nu  = this % P_l_lam / rho
          res % Gam_n = this % GN0 * sqrt(abs(res % E_lam)) * nu &
               * chi2(i_tabn, this % AMUN)
        end if
      else
        call exit_status(EXIT_FAILURE, 'Non-positive neutron width sampled')
        return
      end if

      ! fission width
      if (this % AMUF > 0) then
        res % Gam_f = this % GF * chi2(i_tabf, this % AMUF) / dble(this % AMUF)
      else
        res % Gam_f = ZERO
      end if

      ! constant radiative width
      ! (many channels --> many degrees of freedom --> Dirac delta)
      res % Gam_g = this % GG

      ! competitive width
      if (this % AMUX > 0) then
        res % Gam_x = this % GX * chi2(i_tabx, this % AMUX) / dble(this % AMUX)
      else
        res % Gam_x = ZERO
      end if

      ! total width (sum of partials)
      res % Gam_t = res % Gam_n + res % Gam_f + res % Gam_g + res % Gam_x

      if (res % i_res /= 0) then
        this % local_realization(i_l) % J(i_J) % res(res % i_res) % GN&
             = res % Gam_n
        this % local_realization(i_l) % J(i_J) % res(res % i_res) % GF&
             = res % Gam_f
        this % local_realization(i_l) % J(i_J) % res(res % i_res) % GG&
             = res % Gam_g
        this % local_realization(i_l) % J(i_J) % res(res % i_res) % GX&
             = res % Gam_x
        this % local_realization(i_l) % J(i_J) % res(res % i_res) % GT&
             = res % Gam_t

      end if
    end if

  end subroutine channel_widths


!> Sample resonance parameters for the next resonance added to the ladder
  subroutine resonance_parameters(this, res, i_l, i_J)

    class(Isotope), intent(inout) :: this ! isotope object
    type(Resonance), intent(inout) :: res ! pseudo-resonance object
    integer :: i_l ! orbital angular momentum quantum number index
    integer :: i_J ! total angular momentum quantum number index

    ! sample unresolved resonance parameters for this resonance
    call this % level_spacing(res, i_l, i_J)
    call this % channel_widths(res, i_l, i_J)

  end subroutine resonance_parameters


!> Determine the number of resonances to include from the lth wave that
!! contribute to a URR cross section value
  function num_contributing_resonances(L) result(num_resonances)

    integer :: L              ! orbital quantum number
    integer :: num_resonances ! number of resonances to include for this l

    select case(L)
    case(0)
      num_resonances = num_l_waves(1)
    case(1)
      num_resonances = num_l_waves(2)
    case(2)
      num_resonances = num_l_waves(3)
    case(3)
      call exit_status(EXIT_FAILURE, 'Only s, p, and d wave resonances are supported &
           & in ENDF-6')
      return
      num_resonances = num_l_waves(4)
    case default
      call exit_status(EXIT_FAILURE, 'Only s, p, and d wave resonances are supported &
           & in ENDF-6')
      return
    end select

  end function num_contributing_resonances


!> Calculate single-level Breit-Wigner cross sections
  subroutine slbw_xs(this, res)

    class(Isotope), intent(inout) :: this ! isotope object
    type(Resonance), intent(inout) :: res ! resonance object

    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: theta   ! total width / Doppler width
    real(8) :: x       ! derived variable
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: sig_lam ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi ! derived variable
    real(8) :: S_l_lam ! resonance energy shift factor at E_lam

    ! if URR parameters have resonance energy dependence
    if (parameter_energy_dependence == E_RESONANCE) then

      this % k_lam = wavenumber(this % AWR, abs(res % E_lam))
      this % P_l_lam = penetration(this % L,&
           this % k_lam * this % ac(this % i_urr))
      S_l_lam = energy_shift(this % L, this % k_lam * this % ac(this % i_urr))
      E_shift = res % E_lam &
           + res % Gam_n * (S_l_lam - this % S_l_n) &
           / (TWO * this % P_l_lam)

      Gam_n_n = res % Gam_n * this % P_l_n / this % P_l_lam

      if (this % E > (ONE + ENDF_PRECISION) * this % E_ex2) then
        ! two competitive rxns possible, can't calculate an energy-dependent
        ! width because it depends on the two (unprovided) rxn partial widths
        Gam_x_n = res % Gam_x
      else if (this % E >= this % E_ex1) then
        ! can compute an energy-dependent width for the one competitive reaction
        k_n_x = wavenumber(this % AWR, this % E - this % E_ex1)
        k_lam_x = wavenumber(this % AWR, abs(res % E_lam - this % E_ex1))
        Gam_x_n = res % Gam_x &
             * penetration(this % L, k_n_x   * this % ac(this % i_urr)) &
             / penetration(this % L, k_lam_x * this % ac(this % i_urr))
      else
        Gam_x_n = ZERO
      end if

    else
      ! assume all URR parameters already have neutron energy dependence

      E_shift = res % E_lam
      Gam_n_n = res % Gam_n
      Gam_x_n = res % Gam_x

    end if

    Gam_t_n = res % Gam_t - res % Gam_n - res % Gam_x &
         + Gam_n_n + Gam_x_n

    theta = HALF * Gam_t_n &
         / sqrt(K_BOLTZMANN * 1.0E6_8 * this % T * this % E / this % AWR)

    x = (TWO * (this % E - E_shift)) / Gam_t_n

    sig_lam = FOUR * PI / (this % k_n * this % k_n) * this % g_J &
         * Gam_n_n / Gam_t_n

    ! this particular form comes from the NJOY2012 manual
    res % xs_contribution % n = sig_lam * &
         ((cos(TWO * this % phi_l_n) &
         - (ONE - Gam_n_n / Gam_t_n)) * psi(this % T, theta, x) &
         + sin(TWO * this % phi_l_n) * chi(this % T, theta, x))

    sig_lam_Gam_t_n_psi = sig_lam * psi(this%T, theta, x) / Gam_t_n

    if (res % Gam_g > ZERO) then
      res % xs_contribution % g = sig_lam_Gam_t_n_psi * res % Gam_g
    else
      res % xs_contribution % g = ZERO
    end if

    if (res % Gam_f > ZERO) then
      res % xs_contribution % f = sig_lam_Gam_t_n_psi * res % Gam_f
    else
      res % xs_contribution % f = ZERO
    end if

    if (Gam_x_n > ZERO .and. this % E >= this % E_ex1) then
      res % xs_contribution % x = sig_lam_Gam_t_n_psi * Gam_x_n
    else
      res % xs_contribution % x = ZERO
    end if

    res % xs_contribution % t&
         = res % xs_contribution % n&
         + res % xs_contribution % g&
         + res % xs_contribution % f&
         + res % xs_contribution % x

  end subroutine slbw_xs


!> Calculate a multi-level Breit-Wigner elastic scattering cross section and
!! single-level Breit-Wigner cross sections for other reactions at E_n
  subroutine mlbw_xs(this, res)

    class(Isotope), intent(inout) :: this ! isotope object
    type(Resonance), intent(inout) :: res ! resonance object

    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: theta   ! total width / Doppler width
    real(8) :: x       ! derived variable
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: sig_lam ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi ! derived variable
    real(8) :: S_l_lam ! resonance energy shift factor at E_lam

    if (this % J == this % SPI(this % i_urr)) then
      call exit_status(EXIT_FAILURE, 'Computing MLBW elastic scattering cross section&
           & for resonance with total orbital angular momentum quantum&
           & number, J, equal to the nuclear spin, I')
      return
    end if

    ! if URR parameters have resonance energy dependence
    if (parameter_energy_dependence == E_RESONANCE) then

      this % k_lam = wavenumber(this % AWR, abs(res % E_lam))
      this % P_l_lam = penetration(this % L,&
           this % k_lam * this % ac(this % i_urr))
      S_l_lam = energy_shift(this % L, this % k_lam * this % ac(this % i_urr))
      E_shift = res % E_lam &
           + res % Gam_n * (S_l_lam - this % S_l_n) &
           / (TWO * this % P_l_lam)

      Gam_n_n = res % Gam_n * this % P_l_n / this % P_l_lam

      if (this % E > (ONE + ENDF_PRECISION) * this % E_ex2) then
        ! two competitive rxns possible, can't calculate an energy-dependent
        ! width because it depends on the two (unprovided) rxn partial widths
        Gam_x_n = res % Gam_x
      else if (this % E >= this % E_ex1) then
        ! can compute an energy-dependent width for the one competitive reaction
        k_n_x = wavenumber(this % AWR, this % E - this % E_ex1)
        k_lam_x = wavenumber(this % AWR, abs(res % E_lam - this % E_ex1))
        Gam_x_n = res % Gam_x &
             * penetration(this % L, k_n_x   * this % ac(this % i_urr)) &
             / penetration(this % L, k_lam_x * this % ac(this % i_urr))
      else
        Gam_x_n = ZERO
      end if

    else
      ! assume all URR parameters already have neutron energy dependence
      E_shift = res % E_lam
      Gam_n_n = res % Gam_n
      Gam_x_n = res % Gam_x

    end if

    Gam_t_n = res % Gam_t - res % Gam_n - res % Gam_x &
         + Gam_n_n + Gam_x_n

    theta = HALF * Gam_t_n &
         / sqrt(K_BOLTZMANN * 1.0E6_8 * this % T * this % E / this % AWR)

    x = (TWO * (this % E - E_shift)) / Gam_t_n

    sig_lam = FOUR * PI / (this % k_n * this % k_n) * this % g_J &
         * Gam_n_n / Gam_t_n

    ! this particular form comes from the NJOY2012 manual
    res % xs_contribution % n = sig_lam *&
         ((cos(TWO * this % phi_l_n)&
         - (ONE - Gam_n_n / Gam_t_n)&
         + HALF * this % G_func(E_shift, Gam_n_n, Gam_t_n, res % i_res)&
         / Gam_n_n) * psi(this % T, theta, x)&
         + (sin(TWO * this % phi_l_n)&
         + this % H_func(E_shift, Gam_n_n, Gam_t_n, res % i_res)&
         / Gam_n_n) * chi(this % T, theta, x))

    sig_lam_Gam_t_n_psi = sig_lam * psi(this % T, theta, x) / Gam_t_n

    if (res % Gam_g > ZERO) then
      res % xs_contribution % g = sig_lam_Gam_t_n_psi * res % Gam_g
    else
      res % xs_contribution % g = ZERO
    end if

    if (res % Gam_f > ZERO) then
      res % xs_contribution % f = sig_lam_Gam_t_n_psi * res % Gam_f
    else
      res % xs_contribution % f = ZERO
    end if

    if (Gam_x_n > ZERO .and. this % E >= this % E_ex1) then
      res % xs_contribution % x = sig_lam_Gam_t_n_psi * Gam_x_n
    else
      res % xs_contribution % x = ZERO
    end if

    res % xs_contribution % t&
         = res % xs_contribution % n&
         + res % xs_contribution % g&
         + res % xs_contribution % f&
         + res % xs_contribution % x

  end subroutine mlbw_xs


!> Wrapper for the calculation of partial cross sections
  subroutine calc_xs(this, res)

    class(Isotope), intent(inout) :: this ! isotope object
    type(Resonance), intent(inout) :: res ! resonance object

    ! calculate cross section contributions from an additional resonance
    select case(formalism)
    
    case (SLBW)
      call this % slbw_xs(res)

    case (MLBW)
      call this % mlbw_xs(res)

    case (REICH_MOORE)
      call exit_status(EXIT_FAILURE, 'Reich-Moore formalism not yet supported for the URR')
      return

    case (MNBW)
      call exit_status(EXIT_FAILURE, 'MNBW formalism not yet supported for the URR')
      return

    case default
      call exit_status(EXIT_FAILURE, 'Unrecognized URR resonance formalism')
      return

    end select

  end subroutine calc_xs


end module URR_isotope
