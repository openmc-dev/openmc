module URR_settings

  implicit none
  private
  public :: use_urr,&
            competitive_structure,&
            num_isotopes,&
            xs_representation,&
            formalism,&
            background_xs_treatment,&
            parameter_energy_dependence,&
            num_l_waves,&
            faddeeva_method,&
            write_avg_xs,&
            write_prob_tables,&
            pregenerated_prob_tables,&
            num_energies_prob_tables,&
            num_temperatures_prob_tables,&
            num_bands_prob_tables,&
            num_histories_prob_tables,&
            min_num_batches_prob_tables,&
            max_num_batches_prob_tables,&
            E_grid_scheme_prob_tables,&
            temperature_interp_scheme,&
            rel_err_tolerance_avg_xs,&
            E_grid_prob_tables,&
            T_grid_prob_tables,&
            realization_frequency,&
            num_urr_realizations,&
            i_realization,&
            i_realization_user,&
            min_delta_E_pointwise,&
            rel_err_tolerance_pointwise,&
            path_avg_xs,&
            path_prob_tables,&
            path_endf_files,&
            endf_filenames

  ! general
  logical :: use_urr               ! use any of the URR methods?
  logical :: competitive_structure ! use competitve reaction xs resonance structure?
  integer :: num_isotopes                ! # URR isotopes being processed
  integer :: xs_representation           ! URR xs representation scheme
  integer :: formalism                   ! URR resonance formalism
  integer :: background_xs_treatment     ! where to get background cross sections
  integer :: parameter_energy_dependence ! URR resonance parameter interp scheme
  integer :: num_l_waves(4)              ! # contributing l-wave resonances
  integer :: faddeeva_method             ! Faddeeva function evaluation method

  ! probability tables/avg xs
  logical :: write_avg_xs             ! write averaged xs values to a file?
  logical :: write_prob_tables        ! write probability tables to a file?
  logical :: pregenerated_prob_tables ! load tables or generate new ones?
  integer :: num_energies_prob_tables     ! # prob table energies
  integer :: num_temperatures_prob_tables ! # prob table temperatures
  integer :: num_bands_prob_tables        ! # prob table xs bands
  integer :: num_histories_prob_tables    ! # histories/batch in prob table generation
  integer :: min_num_batches_prob_tables  ! min batches in prob table generation
  integer :: max_num_batches_prob_tables  ! max batches in prob table generation
  integer :: E_grid_scheme_prob_tables    ! prob table energy spacing scheme
  integer :: temperature_interp_scheme    ! URR temperature interpolation scheme
  real(8) :: rel_err_tolerance_avg_xs ! max rel err for avg xs calc termination
  real(8), allocatable :: E_grid_prob_tables(:) ! probability table energies
  real(8), allocatable :: T_grid_prob_tables(:) ! probability table temperatures

  ! fixed realizations
  integer :: realization_frequency ! frequency of URR realizations
  integer :: num_urr_realizations  ! # independent URR realizations
  integer :: i_realization         ! index of URR realization used for calc
  integer :: i_realization_user    ! user-specified realization index

  ! pointwise
  real(8) :: min_delta_E_pointwise       ! min diff between pointwise xs energies [eV]
  real(8) :: rel_err_tolerance_pointwise ! max rel err for pointwise xs reconstruction

  ! files
  character(:), allocatable :: path_avg_xs      ! path to averaged URR xs files
  character(:), allocatable :: path_prob_tables ! path to probability table files
  character(:), allocatable :: path_endf_files  ! path to ENDF-6 files
  character(255), allocatable :: endf_filenames(:) ! list of ENDF-6 filenames

end module URR_settings
