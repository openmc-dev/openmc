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

  implicit none

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

  interface
    subroutine print_particle(p) bind(C)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine
  end interface

contains

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
      do h = 1, size(t % filter)
        filter_bins(t % filter(h)) = 0
      end do
      j = 1
      indent = 0

      print_bin: do
        find_bin: do
          ! Check for no filters
          if (size(t % filter) == 0) exit find_bin

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
            if (j == size(t % filter)) exit find_bin

            ! Print current filter information
            write(UNIT=unit_tally, FMT='(1X,2A)') repeat(" ", indent), &
                 trim(filters(t % filter(j)) % obj % &
                 text_label(filter_bins(t % filter(j))))
            indent = indent + 2
            j = j + 1
          end if

        end do find_bin

        ! Print filter information
        if (size(t % filter) > 0) then
          write(UNIT=unit_tally, FMT='(1X,2A)') repeat(" ", indent), &
               trim(filters(t % filter(j)) % obj % &
               text_label(filter_bins(t % filter(j))))
        end if

        ! Determine scoring index for this bin combination -- note that unlike
        ! in the score_tally subroutine, we have to use max(bins,1) since all
        ! bins below the lowest filter level will be zeros

        filter_index = 1
        do h = 1, size(t % filter)
          filter_index = filter_index &
               + (max(filter_bins(t % filter(h)) ,1) - 1) * t % stride(h)
        end do

        ! Write results for this filter bin combination
        score_index = 0
        if (size(t % filter) > 0) indent = indent + 2
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

        if (size(t % filter) == 0) exit print_bin

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
