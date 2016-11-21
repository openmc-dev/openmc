module URR_settings

  use global, only: MAX_FILE_LEN

  implicit none
  private

  ! logicals
  public :: competitive_structure,&
            load_urr_prob_tables,&
            use_urr,&
            write_avg_urr_xs,&
            write_urr_prob_tables

  ! integers
  public :: background_xs_treatment,&
            E_grid_scheme_urr_prob_tables,&
            formalism_urr,&
            n_histories_urr_prob_tables,&
            INT_T,&
            i_realiz,&
            i_realiz_user,&
            n_l_waves,&
            max_n_batch_urr_prob_tables,&
            min_n_batch_urr_prob_tables,&
            n_bands_urr_prob_tables,&
            n_temperatures_urr_prob_tables,&
            n_isotopes_urr,&
            n_energies_urr,&
            n_realiz_urr,&
            realiz_frequency_urr,&
            parameter_interp_scheme,&
            xs_representation_urr,&
            faddeeva_method

  ! reals
  public :: min_dE_point_urr,&
            tol_avg_xs_urr,&
            tol_point_urr,&
            E_grid_urr_prob_tables,&
            T_grid_urr_prob_tables

  ! characters
  public :: path_avg_urr_xs,&
            path_endf_files,&
            path_urr_prob_tables,&
            endf_filenames

  logical :: competitive_structure ! use competitve reaction xs resonance structure?
  logical :: load_urr_prob_tables  ! load tables or generate new ones?
  logical :: use_urr               ! use any of the URR methods?
  logical :: write_avg_urr_xs      ! write averaged xs values to a file?
  logical :: write_urr_prob_tables ! write probability tables to a file?

  integer :: background_xs_treatment        ! where to get background cross sections
  integer :: E_grid_scheme_urr_prob_tables  ! prob table energy spacing scheme
  integer :: formalism_urr                  ! URR resonance formalism
  integer :: n_histories_urr_prob_tables    ! histories for prob table calc
  integer :: INT_T                          ! URR temperature interpolation scheme
  integer :: i_realiz                       ! index of URR realization used for calc
  integer :: i_realiz_user                  ! user-specified realization index
  integer :: n_l_waves(4)                   ! # contributing l-wave resonances
  integer :: max_n_batch_urr_prob_tables    ! max batches for prob table calc
  integer :: min_n_batch_urr_prob_tables    ! min batches for prob table calc
  integer :: n_bands_urr_prob_tables        ! # prob table xs bands
  integer :: n_temperatures_urr_prob_tables ! # prob table temperatures
  integer :: n_isotopes_urr                 ! # URR isotopes being processed
  integer :: n_energies_urr                 ! # prob table energies
  integer :: n_realiz_urr                   ! # independent URR realizations
  integer :: realiz_frequency_urr           ! frequency of URR realizations
  integer :: parameter_interp_scheme        ! URR resonance parameter interp scheme
  integer :: xs_representation_urr          ! URR xs representation scheme
  integer :: faddeeva_method                ! Faddeeva function evaluation method

  real(8) :: min_dE_point_urr    ! min diff between pointwise URR xs energies [eV]
  real(8) :: tol_avg_xs_urr      ! max rel err for average URR xs calc termination
  real(8) :: tol_point_urr       ! max rel err for pointwise URR xs reconstruction
  real(8), allocatable :: E_grid_urr_prob_tables(:) ! probability table energies
  real(8), allocatable :: T_grid_urr_prob_tables(:) ! probability table temperatures

  character(MAX_FILE_LEN) :: path_avg_urr_xs      ! path to averaged URR xs files
  character(MAX_FILE_LEN) :: path_endf_files      ! path to ENDF-6 files
  character(MAX_FILE_LEN) :: path_urr_prob_tables ! path to probability table files
  character(80), allocatable :: endf_filenames(:) ! list of ENDF-6 filenames

end module URR_settings
