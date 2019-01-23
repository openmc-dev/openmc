module output

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use constants
  use eigenvalue,      only: openmc_get_keff
  use endf,            only: reaction_name
  use error,           only: fatal_error, warning
  use geometry_header
  use math,            only: t_percentile
  use message_passing, only: master, n_procs
  use mgxs_interface
  use nuclide_header
  use particle_header, only: LocalCoord, Particle
  use settings
  use simulation_header
  use surface_header,  only: surfaces
  use string,          only: to_upper, to_str
  use tally_header
  use tally_derivative_header
  use tally_filter
  use tally_filter_mesh, only: MeshFilter
  use tally_filter_header, only: TallyFilterMatch
  use timer_header

  implicit none

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

  interface
    function entropy(i) result(h) bind(C, name='entropy_c')
      import C_INT, C_DOUBLE
      integer(C_INT), value :: i
      real(C_DOUBLE) :: h
    end function

    subroutine print_particle(p) bind(C)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine
  end interface

contains

!===============================================================================
! TITLE prints the main title banner as well as information about the program
! developers, version, and date/time which the problem was run.
!===============================================================================

  subroutine title() bind(C)

#ifdef _OPENMP
    use omp_lib
#endif

    write(UNIT=OUTPUT_UNIT, FMT='(/23(A/))') &
         '                               %%%%%%%%%%%%%%%', &
         '                          %%%%%%%%%%%%%%%%%%%%%%%%', &
         '                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', &
         '                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', &
         '                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', &
         '                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', &
         '                                   %%%%%%%%%%%%%%%%%%%%%%%%', &
         '                                    %%%%%%%%%%%%%%%%%%%%%%%%', &
         '                ###############      %%%%%%%%%%%%%%%%%%%%%%%%', &
         '               ##################     %%%%%%%%%%%%%%%%%%%%%%%', &
         '               ###################     %%%%%%%%%%%%%%%%%%%%%%%', &
         '               ####################     %%%%%%%%%%%%%%%%%%%%%%', &
         '               #####################     %%%%%%%%%%%%%%%%%%%%%', &
         '               ######################     %%%%%%%%%%%%%%%%%%%%', &
         '               #######################     %%%%%%%%%%%%%%%%%%', &
         '                #######################     %%%%%%%%%%%%%%%%%', &
         '                ######################     %%%%%%%%%%%%%%%%%', &
         '                 ####################     %%%%%%%%%%%%%%%%%', &
         '                   #################     %%%%%%%%%%%%%%%%%', &
         '                    ###############     %%%%%%%%%%%%%%%%', &
         '                      ############     %%%%%%%%%%%%%%%', &
         '                         ########     %%%%%%%%%%%%%%', &
         '                                     %%%%%%%%%%%'

    ! Write version information
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '                  | The OpenMC Monte Carlo Code'
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '        Copyright | 2011-2018 MIT and OpenMC contributors'
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '          License | http://openmc.readthedocs.io/en/latest/license.html'
    write(UNIT=OUTPUT_UNIT, FMT='(11X,"Version | ",I1,".",I2,".",I1)') &
         VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#ifdef GIT_SHA1
    write(UNIT=OUTPUT_UNIT, FMT='(10X,"Git SHA1 | ",A)') GIT_SHA1
#endif

    ! Write the date and time
    write(UNIT=OUTPUT_UNIT, FMT='(9X,"Date/Time | ",A)') time_stamp()

#ifdef OPENMC_MPI
    ! Write number of processors
    write(UNIT=OUTPUT_UNIT, FMT='(5X,"MPI Processes | ",A)') &
         trim(to_str(n_procs))
#endif

#ifdef _OPENMP
    ! Write number of OpenMP threads
    write(UNIT=OUTPUT_UNIT, FMT='(4X,"OpenMP Threads | ",A)') &
         trim(to_str(omp_get_max_threads()))
#endif
    write(UNIT=OUTPUT_UNIT, FMT=*)

  end subroutine title

!===============================================================================
! TIME_STAMP returns the current date and time in a formatted string
!===============================================================================

  function time_stamp() result(current_time)

    character(19) :: current_time ! ccyy-mm-dd hh:mm:ss
    character(8)  :: date_        ! ccyymmdd
    character(10) :: time_        ! hhmmss.sss

    call date_and_time(DATE=date_, TIME=time_)
    current_time = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8) // &
         " " // time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  end function time_stamp

!===============================================================================
! HEADER displays a header block according to a specified level. If no level is
! specified, it is assumed to be a minor header block.
!===============================================================================

  subroutine header(msg, level, unit)
    character(*), intent(in)      :: msg   ! header message
    integer, intent(in)           :: level
    integer, intent(in), optional :: unit  ! unit to write to

    integer :: n            ! number of = signs on left
    integer :: m            ! number of = signs on right
    integer :: unit_        ! unit to write to
    character(MAX_LINE_LEN) :: line

    ! set default unit
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! determine how many times to repeat '=' character
    n = (63 - len_trim(msg))/2
    m = n
    if (mod(len_trim(msg),2) == 0) m = m + 1

    ! convert line to upper case
    line = to_upper(msg)

    ! print header based on level
    if (verbosity >= level) then
      write(UNIT=unit_, FMT='(/1X,A/)') repeat('=', n) // '>     ' // &
           trim(line) // '     <' // repeat('=', m)
    end if

  end subroutine header

!===============================================================================
! PRINT_VERSION shows the current version as well as copright and license
! information
!===============================================================================

  subroutine print_version() bind(C)

    if (master) then
      write(UNIT=OUTPUT_UNIT, FMT='(1X,A,1X,I1,".",I2,".",I1)') &
           "OpenMC version", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#ifdef GIT_SHA1
      write(UNIT=OUTPUT_UNIT, FMT='(1X,A,A)') "Git SHA1: ", GIT_SHA1
#endif
      write(UNIT=OUTPUT_UNIT, FMT=*) "Copyright (c) 2011-2018 &
           &Massachusetts Institute of Technology and OpenMC contributors"
      write(UNIT=OUTPUT_UNIT, FMT=*) "MIT/X license at &
           &<http://openmc.readthedocs.io/en/latest/license.html>"
    end if

  end subroutine print_version

!===============================================================================
! PRINT_USAGE displays information about command line usage of OpenMC
!===============================================================================

  subroutine print_usage() bind(C)

    if (master) then
      write(OUTPUT_UNIT,*) 'Usage: openmc [options] [directory]'
      write(OUTPUT_UNIT,*)
      write(OUTPUT_UNIT,*) 'Options:'
      write(OUTPUT_UNIT,*) '  -c, --volume           Run in stochastic volume calculation mode'
      write(OUTPUT_UNIT,*) '  -g, --geometry-debug   Run with geometry debugging on'
      write(OUTPUT_UNIT,*) '  -n, --particles        Number of particles per generation'
      write(OUTPUT_UNIT,*) '  -p, --plot             Run in plotting mode'
      write(OUTPUT_UNIT,*) '  -r, --restart          Restart a previous run from a state point'
      write(OUTPUT_UNIT,*) '                         or a particle restart file'
      write(OUTPUT_UNIT,*) '  -s, --threads          Number of OpenMP threads'
      write(OUTPUT_UNIT,*) '  -t, --track            Write tracks for all particles'
      write(OUTPUT_UNIT,*) '  -v, --version          Show version information'
      write(OUTPUT_UNIT,*) '  -h, --help             Show this message'
    end if

  end subroutine print_usage

!===============================================================================
! PRINT_COLUMNS displays a header listing what physical values will displayed
! below them
!===============================================================================

  subroutine print_columns() bind(C)

    write(UNIT=ou, FMT='(2X,A9,3X)', ADVANCE='NO') "Bat./Gen."
    write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "   k    "
    if (entropy_on) write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "Entropy "
    write(UNIT=ou, FMT='(A20,3X)', ADVANCE='NO') "     Average k      "
    write(UNIT=ou, FMT=*)

    write(UNIT=ou, FMT='(2X,A9,3X)', ADVANCE='NO') "========="
    write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "========"
    if (entropy_on) write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "========"
    write(UNIT=ou, FMT='(A20,3X)', ADVANCE='NO') "===================="
    write(UNIT=ou, FMT=*)

  end subroutine print_columns

!===============================================================================
! PRINT_GENERATION displays information for a generation of neutrons.
!===============================================================================

  subroutine print_generation() bind(C)

    integer :: i  ! overall generation
    integer :: n  ! number of active generations

    ! Determine overall generation and number of active generations
    i = overall_generation()
    if (current_batch > n_inactive) then
      n = gen_per_batch*n_realizations + current_gen
    else
      n = 0
    end if

    ! write out information about batch and generation
    write(UNIT=OUTPUT_UNIT, FMT='(2X,A9)', ADVANCE='NO') &
         trim(to_str(current_batch)) // "/" // trim(to_str(current_gen))
    write(UNIT=OUTPUT_UNIT, FMT='(3X,F8.5)', ADVANCE='NO') k_generation(i)

    ! write out entropy info
    if (entropy_on) write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
         entropy(i)

    if (n > 1) then
      write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5," +/-",F8.5)', ADVANCE='NO') &
           keff, keff_std
    end if

    ! next line
    write(UNIT=OUTPUT_UNIT, FMT=*)

  end subroutine print_generation

!===============================================================================
! PRINT_BATCH_KEFF displays the last batch's tallied value of the neutron
! multiplication factor as well as the average value if we're in active batches
!===============================================================================

  subroutine print_batch_keff() bind(C)

    integer :: i  ! overall generation
    integer :: n  ! number of active generations

    ! Determine overall generation and number of active generations
    i = current_batch*gen_per_batch
    n = n_realizations*gen_per_batch

    ! write out information batch and option independent output
    write(UNIT=OUTPUT_UNIT, FMT='(2X,A9)', ADVANCE='NO') &
         trim(to_str(current_batch)) // "/" // trim(to_str(gen_per_batch))
    write(UNIT=OUTPUT_UNIT, FMT='(3X,F8.5)', ADVANCE='NO') &
         k_generation(i)

    ! write out entropy info
    if (entropy_on) write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
         entropy(i)

    ! write out accumulated k-effective if after first active batch
    if (n > 1) then
      write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5," +/-",F8.5)', ADVANCE='NO') &
           keff, keff_std
    else
      write(UNIT=OUTPUT_UNIT, FMT='(23X)', ADVANCE='NO')
    end if

    ! next line
    write(UNIT=OUTPUT_UNIT, FMT=*)

  end subroutine print_batch_keff

!===============================================================================
! PRINT_RUNTIME displays the total time elapsed for the entire run, for
! initialization, for computation, and for intergeneration synchronization.
!===============================================================================

  subroutine print_runtime() bind(C)

    integer       :: n_active
    real(8)       :: speed_inactive  ! # of neutrons/second in inactive batches
    real(8)       :: speed_active    ! # of neutrons/second in active batches
    character(15) :: string

    ! display header block
    call header("Timing Statistics", 6)

    ! display time elapsed for various sections
    write(ou,100) "Total time for initialization", time_initialize_elapsed()
    write(ou,100) "  Reading cross sections", time_read_xs % elapsed
    write(ou,100) "Total time in simulation", time_inactive_elapsed() + &
         time_active_elapsed()
    write(ou,100) "  Time in transport only", time_transport_elapsed()
    if (run_mode == MODE_EIGENVALUE) then
      write(ou,100) "  Time in inactive batches", time_inactive_elapsed()
    end if
    write(ou,100) "  Time in active batches", time_active_elapsed()
    if (run_mode == MODE_EIGENVALUE) then
      write(ou,100) "  Time synchronizing fission bank", time_bank_elapsed()
      write(ou,100) "    Sampling source sites", time_bank_sample_elapsed()
      write(ou,100) "    SEND/RECV source sites", time_bank_sendrecv_elapsed()
    end if
    write(ou,100) "  Time accumulating tallies", time_tallies_elapsed()
    write(ou,100) "Total time for finalization", time_finalize_elapsed()
    write(ou,100) "Total time elapsed", time_total_elapsed()

    ! Calculate particle rate in active/inactive batches
    n_active = current_batch - n_inactive
    if (restart_run) then
      if (restart_batch < n_inactive) then
        speed_inactive = real(n_particles * (n_inactive - restart_batch) * &
             gen_per_batch) / time_inactive_elapsed()
        speed_active = real(n_particles * n_active * gen_per_batch) / &
             time_active_elapsed()
      else
        speed_inactive = ZERO
        speed_active = real(n_particles * (n_batches - restart_batch) * &
             gen_per_batch) / time_active_elapsed()
      end if
    else
      if (n_inactive > 0) then
        speed_inactive = real(n_particles * n_inactive * gen_per_batch) / &
             time_inactive_elapsed()
      end if
      speed_active = real(n_particles * n_active * gen_per_batch) / &
           time_active_elapsed()
    end if

    ! display calculation rate
    if (.not. (restart_run .and. (restart_batch >= n_inactive)) &
         .and. n_inactive > 0) then
      string = to_str(speed_inactive)
      write(ou,101) "Calculation Rate (inactive)", trim(string)
    end if
    string = to_str(speed_active)
    write(ou,101) "Calculation Rate (active)", trim(string)

    ! format for write statements
100 format (1X,A,T36,"= ",ES11.4," seconds")
101 format (1X,A,T36,"=  ",A," particles/second")

  end subroutine print_runtime

!===============================================================================
! PRINT_RESULTS displays various estimates of k-effective as well as the global
! leakage rate.
!===============================================================================

  subroutine print_results() bind(C)

    integer :: n       ! number of realizations
    real(8) :: alpha   ! significance level for CI
    real(8) :: t_n1    ! t-value with N-1 degrees of freedom
    real(8) :: t_n3    ! t-value with N-3 degrees of freedom
    real(8) :: x(2)    ! mean and standard deviation
    real(C_DOUBLE) :: k_combined(2)
    integer(C_INT) :: err

    ! display header block for results
    call header("Results", 4)

    n = n_realizations

    if (confidence_intervals) then
      ! Calculate t-value for confidence intervals
      alpha = ONE - CONFIDENCE_LEVEL
      t_n1 = t_percentile(ONE - alpha/TWO, n - 1)
      t_n3 = t_percentile(ONE - alpha/TWO, n - 3)
    else
      t_n1 = ONE
      t_n3 = ONE
    end if

    ! write global tallies
    if (n > 1) then
      associate (r => global_tallies(RESULT_SUM:RESULT_SUM_SQ, :))
        if (run_mode == MODE_EIGENVALUE) then
          x(:) = mean_stdev(r(:, K_COLLISION), n)
          write(ou,102) "k-effective (Collision)", x(1), t_n1 * x(2)
          x(:) = mean_stdev(r(:, K_TRACKLENGTH), n)
          write(ou,102) "k-effective (Track-length)", x(1), t_n1 * x(2)
          x(:) = mean_stdev(r(:, K_ABSORPTION), n)
          write(ou,102) "k-effective (Absorption)", x(1), t_n1 * x(2)
          if (n > 3) then
            err = openmc_get_keff(k_combined)
            write(ou,102) "Combined k-effective", k_combined(1), &
                 t_n3 * k_combined(2)
          end if
        end if
        x(:) = mean_stdev(r(:, LEAKAGE), n)
        write(ou,102) "Leakage Fraction", x(1), t_n1 * x(2)
      end associate
    else
      if (master) call warning("Could not compute uncertainties -- only one &
           &active batch simulated!")

      if (run_mode == MODE_EIGENVALUE) then
        write(ou,103) "k-effective (Collision)", global_tallies(RESULT_SUM, K_COLLISION) / n
        write(ou,103) "k-effective (Track-length)", global_tallies(RESULT_SUM, K_TRACKLENGTH) / n
        write(ou,103) "k-effective (Absorption)", global_tallies(RESULT_SUM, K_ABSORPTION) / n
      end if
      write(ou,103) "Leakage Fraction", global_tallies(RESULT_SUM, LEAKAGE) / n
    end if
    write(ou,*)

102 format (1X,A,T30,"= ",F8.5," +/- ",F8.5)
103 format (1X,A,T30,"= ",F8.5)

  end subroutine print_results

!===============================================================================
! WRITE_TALLIES creates an output file and writes out the mean values of all
! tallies and their standard deviations
!===============================================================================

  subroutine write_tallies() bind(C)

    integer :: i            ! index in tallies array
    integer :: j            ! level in tally hierarchy
    integer :: k            ! loop index for scoring bins
    integer :: n            ! loop index for nuclides
    integer :: h            ! loop index for tally filters
    integer :: indent       ! number of spaces to preceed output
    integer :: filter_index ! index in results array for filters
    integer :: score_index  ! scoring bin index
    integer :: i_nuclide    ! index in nuclides array
    integer :: i_filt
    integer :: unit_tally   ! tallies.out file unit
    integer :: nr           ! number of realizations
    real(8) :: t_value      ! t-values for confidence intervals
    real(8) :: alpha        ! significance level for CI
    real(8) :: x(2)         ! mean and standard deviation
    character(MAX_FILE_LEN) :: filename                    ! name of output file
    character(36)           :: score_names(N_SCORE_TYPES)  ! names of scoring function
    character(36)           :: score_name                  ! names of scoring function
                                                           ! to be applied at write-time
    integer, allocatable :: filter_bins(:)
    character(MAX_WORD_LEN) :: temp_name

    ! Skip if there are no tallies
    if (n_tallies == 0) return

    allocate(filter_bins(n_filters))

    ! Initialize names for scores
    score_names(abs(SCORE_FLUX))               = "Flux"
    score_names(abs(SCORE_TOTAL))              = "Total Reaction Rate"
    score_names(abs(SCORE_SCATTER))            = "Scattering Rate"
    score_names(abs(SCORE_NU_SCATTER))         = "Scattering Production Rate"
    score_names(abs(SCORE_ABSORPTION))         = "Absorption Rate"
    score_names(abs(SCORE_FISSION))            = "Fission Rate"
    score_names(abs(SCORE_NU_FISSION))         = "Nu-Fission Rate"
    score_names(abs(SCORE_KAPPA_FISSION))      = "Kappa-Fission Rate"
    score_names(abs(SCORE_EVENTS))             = "Events"
    score_names(abs(SCORE_DECAY_RATE))         = "Decay Rate"
    score_names(abs(SCORE_DELAYED_NU_FISSION)) = "Delayed-Nu-Fission Rate"
    score_names(abs(SCORE_PROMPT_NU_FISSION))  = "Prompt-Nu-Fission Rate"
    score_names(abs(SCORE_INVERSE_VELOCITY))   = "Flux-Weighted Inverse Velocity"
    score_names(abs(SCORE_FISS_Q_PROMPT))      = "Prompt fission power"
    score_names(abs(SCORE_FISS_Q_RECOV))       = "Recoverable fission power"
    score_names(abs(SCORE_CURRENT))            = "Current"

    ! Create filename for tally output
    filename = trim(path_output) // "tallies.out"

    ! Open tally file for writing
    open(FILE=filename, NEWUNIT=unit_tally, STATUS='replace', ACTION='write')

    ! Calculate t-value for confidence intervals
    if (confidence_intervals) then
      alpha = ONE - CONFIDENCE_LEVEL
      t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)
    else
      t_value = ONE
    end if

    TALLY_LOOP: do i = 1, n_tallies
      associate (t => tallies(i) % obj)
      nr = t % n_realizations

      if (confidence_intervals) then
        ! Calculate t-value for confidence intervals
        alpha = ONE - CONFIDENCE_LEVEL
        t_value = t_percentile(ONE - alpha/TWO, nr - 1)
      else
        t_value = ONE
      end if

      ! Write header block
      if (t % name == "") then
        call header("TALLY " // trim(to_str(t % id)), 1, unit=unit_tally)
      else
        call header("TALLY " // trim(to_str(t % id)) // ": " &
             // trim(t % name), 1, unit=unit_tally)
      endif

      ! Write derivative information.
      if (t % deriv /= NONE) then
        associate(deriv => tally_derivs(t % deriv))
          select case (deriv % variable)
          case (DIFF_DENSITY)
            write(unit=unit_tally, fmt="(' Density derivative  Material ',A)") &
                 to_str(deriv % diff_material)
          case (DIFF_NUCLIDE_DENSITY)
            write(unit=unit_tally, fmt="(' Nuclide density derivative  &
                 &Material ',A,'  Nuclide ',A)") &
                 trim(to_str(deriv % diff_material)), &
                 trim(nuclides(deriv % diff_nuclide) % name)
          case (DIFF_TEMPERATURE)
            write(unit=unit_tally, fmt="(' Temperature derivative  Material ',&
                 &A)") to_str(deriv % diff_material)
          case default
            call fatal_error("Differential tally dependent variable for tally "&
                 // trim(to_str(t % id)) // " not defined in output.F90.")
          end select
        end associate
      end if

      ! WARNING: Admittedly, the logic for moving for printing results is
      ! extremely confusing and took quite a bit of time to get correct. The
      ! logic is structured this way since it is not practical to have a do
      ! loop for each filter variable (given that only a few filters are likely
      ! to be used for a given tally.

      ! Initialize bins, filter level, and indentation
      do h = 1, t % n_filters()
        filter_bins(t % filter(h)) = 0
      end do
      j = 1
      indent = 0

      print_bin: do
        find_bin: do
          ! Check for no filters
          if (t % n_filters() == 0) exit find_bin

          ! Increment bin combination
          filter_bins(t % filter(j)) = filter_bins(t % filter(j)) + 1

          ! =================================================================
          ! REACHED END OF BINS FOR THIS FILTER, MOVE TO NEXT FILTER

          if (filter_bins(t % filter(j)) > &
               filters(t % filter(j)) % obj % n_bins) then
            ! If this is the first filter, then exit
            if (j == 1) exit print_bin

            filter_bins(t % filter(j)) = 0
            j = j - 1
            indent = indent - 2

            ! =================================================================
            ! VALID BIN -- WRITE FILTER INFORMATION OR EXIT TO WRITE RESULTS

          else
            ! Check if this is last filter
            if (j == t % n_filters()) exit find_bin

            ! Print current filter information
            i_filt = t % filter(j)
            write(UNIT=unit_tally, FMT='(1X,2A)') repeat(" ", indent), &
                 trim(filters(i_filt) % obj % &
                 text_label(filter_bins(i_filt)))
            indent = indent + 2
            j = j + 1
          end if

        end do find_bin

        ! Print filter information
        if (t % n_filters() > 0) then
          i_filt = t % filter(j)
          write(UNIT=unit_tally, FMT='(1X,2A)') repeat(" ", indent), &
               trim(filters(i_filt) % obj % &
               text_label(filter_bins(i_filt)))
        end if

        ! Determine scoring index for this bin combination -- note that unlike
        ! in the score_tally subroutine, we have to use max(bins,1) since all
        ! bins below the lowest filter level will be zeros

        filter_index = 1
        do h = 1, t % n_filters()
          filter_index = filter_index &
               + (max(filter_bins(t % filter(h)) ,1) - 1) * t % stride(h)
        end do

        ! Write results for this filter bin combination
        score_index = 0
        if (t % n_filters() > 0) indent = indent + 2
        do n = 1, t % n_nuclide_bins
          ! Write label for nuclide
          i_nuclide = t % nuclide_bins(n)
          if (i_nuclide == -1) then
            write(UNIT=unit_tally, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                 "Total Material"
          else
            if (run_CE) then
              write(UNIT=unit_tally, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                   trim(nuclides(i_nuclide) % name)
            else
              call get_name_c(i_nuclide, len(temp_name), temp_name)
              write(UNIT=unit_tally, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                   trim(temp_name)
            end if
          end if

          indent = indent + 2
          do k = 1, t % n_score_bins
            score_index = score_index + 1

            associate(r => t % results(RESULT_SUM:RESULT_SUM_SQ, :, :))

            if (t % score_bins(k) > 0) then
              score_name = reaction_name(t % score_bins(k))
            else
              score_name = score_names(abs(t % score_bins(k)))
            end if
            x(:) = mean_stdev(r(:, score_index, filter_index), nr)
            write(UNIT=unit_tally, FMT='(1X,2A,1X,A,"+/- ",A)') &
                 repeat(" ", indent), score_name, &
                 to_str(x(1)), trim(to_str(t_value * x(2)))
            end associate
          end do
          indent = indent - 2

        end do
        indent = indent - 2

        if (t % n_filters() == 0) exit print_bin

      end do print_bin

      end associate
    end do TALLY_LOOP

    close(UNIT=unit_tally)

  end subroutine write_tallies

!===============================================================================
! MEAN_STDEV computes the sample mean and standard deviation of the mean of a
! single tally score
!===============================================================================

  pure function mean_stdev(result_, n) result(x)
    real(8), intent(in) :: result_(2) ! sum and sum-of-squares
    integer, intent(in) :: n          ! number of realizations
    real(8)  :: x(2)                  ! mean and standard deviation

    x(1) = result_(1) / n
    if (n > 1) then
      x(2) = sqrt((result_(2) / n - x(1)*x(1))/(n - 1))
    else
      x(2) = ZERO
    end if
  end function mean_stdev

end module output
