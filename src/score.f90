module score

  use global
  use output, only: message, error
  use ace,    only: get_macro_xs
  use search, only: binary_search

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!=====================================================================
! CALCULATE_KEFF
!=====================================================================

  subroutine calculate_keff(i_cycle)

    integer, intent(in) :: i_cycle ! index of current cycle

    integer(8) :: total_bank
    integer :: n
    integer :: ierr
    real(8) :: kcoll   ! keff collision estimator         
    real(8) :: ktemp   ! MPI-reduced keff and stdev
    real(8), save :: k1 = 0. ! accumulated keff
    real(8), save :: k2 = 0. ! accumulated keff**2
    real(8) :: std     ! stdev of keff over active cycles
    character(250) :: msg

    msg = "Calculate cycle keff..."
    call message(msg, 8)

    ! set k1 and k2 at beginning of run
    if (i_cycle == 1) then
       k1 = ZERO
       k2 = ZERO
    end if

#ifdef MPI
    ! Collect number bank sites onto master process
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         & MPI_COMM_WORLD, ierr)
#else
    total_bank = n_bank
#endif

    ! Collect statistics and print output
    if (master) then
       kcoll = real(total_bank)/real(n_particles)*keff
       if (i_cycle > n_inactive) then
          n = i_cycle - n_inactive
          k1 = k1 + kcoll
          k2 = k2 + kcoll**2
          keff = k1/n
          std  = sqrt((k2/n-keff**2)/n)
          if (i_cycle > n_inactive+1) then
             write(6,101) i_cycle, kcoll, keff, std
          else
             write(6,100) i_cycle, kcoll
          end if
       else
          write(6,100) i_cycle, kcoll
          keff = kcoll
       end if
    end if

#ifdef MPI
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
#endif

100 format (2X,I4,3X,F7.5)
101 format (2X,I4,3X,F7.5,10X,F7.5,2X,F7.5)

  end subroutine calculate_keff

!=====================================================================
! SCORE_TALLY
!=====================================================================

  subroutine score_tally(p, flux)

    type(Particle), pointer :: p     ! particle
    real(8), intent(in)     :: flux  ! estimator of flux

    integer :: i          ! index for tallies in cell
    integer :: j          ! index for reactions in tally
    integer :: e_bin      ! energy bin
    integer :: c_bin      ! cell bin
    integer :: r_bin      ! reaction bin
    integer :: n_energy   ! number of energy bins
    integer :: n_reaction ! number of reactions
    integer :: type       ! cell or reaction type
    integer :: MT         ! reaction MT value
    real(8) :: E          ! energy of particle
    real(8) :: val        ! value to score
    real(8) :: Sigma      ! macroscopic cross section of reaction
    character(250) :: msg ! output/error message
    type(Cell),  pointer :: c => null()
    type(Tally), pointer :: t => null()

    ! ================================================================
    ! HANDLE LOCAL TALLIES

    c => cells(p % cell)
    if (.not. allocated(c % tallies)) return

    ! if so, loop over each tally
    do i = 1, size(c % tallies)
       t => tallies(c % tallies(i))

       ! =============================================================
       ! DETERMINE ENERGY BIN
       if (allocated(t % energies)) then
          E = p % E
          n_energy = size(t % energies)
          if (E < t % energies(1) .or. E > t % energies(n_energy)) then
             ! Energy outside of specified range
             cycle
          else
             e_bin = binary_search(t % energies, n_energy, E)
          end if
       end if
 
       ! =============================================================
       ! DETERMINE CELL BIN
       type = t % cell_type
       if (type == TALLY_SUM) then
          ! Sum tallies from separate cells into one bin
          c_bin = 1
       elseif (type == TALLY_BINS) then
          ! Need to determine cell bin
          do j = 1, size(t % cells)
             if (p % cell == t % cells(j)) then
                c_bin = j
                exit
             end if
          end do
       else
          msg = "Invalid type for cell bins in tally " // int_to_str(t % uid)
          call error(msg)
       end if

       ! =============================================================
       ! DETERMINE REACTION BIN AND ADD VALUE TO SCORE
       type = t % reaction_type
       if (type == TALLY_FLUX) then
          ! Tally flux only
          val = flux
          call add_to_score(t % score(r_bin, c_bin, e_bin), val)

       elseif (type == TALLY_ALL) then
          ! Tally total reaction rate
          val = ONE
          r_bin = 1
          call add_to_score(t % score(r_bin, c_bin, e_bin), val)

       elseif (type == TALLY_BINS) then
          ! Individually bin reactions
          n_reaction = t % reaction_type
          
          r_bin = 0
          do j = 1, n_reaction
             MT = t % reactions(j)
             Sigma = get_macro_xs(p, cMaterial, MT)
             val = Sigma * flux
             r_bin = r_bin + 1
             call add_to_score(t % score(r_bin, c_bin, e_bin), &
                  & val)
          end do
       elseif (type == TALLY_SUM) then
          ! Tally reactions in one bin
          n_reaction = t % reaction_type
          
          r_bin = 1
          do j = 1, n_reaction
             MT = t % reactions(j)
             Sigma = get_macro_xs(p, cMaterial, MT)
             val = Sigma * flux
             call add_to_score(t % score(r_bin, c_bin, e_bin), &
                  & val)
          end do
       end if
    end do

    ! ================================================================
    ! HANDLE GLOBAL TALLIES

    do i = 1, n_tallies_global
       t => tallies_global(i)

       ! =============================================================
       ! DETERMINE ENERGY BIN
       if (allocated(t % energies)) then
          E = p % E
          n_energy = size(t % energies)
          if (E < t % energies(1) .or. E > t % energies(n_energy)) then
             ! Energy outside of specified range
             cycle
          else
             e_bin = binary_search(t % energies, n_energy, E)
          end if
       end if

       ! Since it's a global tally, the cell bin is unity
       c_bin = 1
       
       ! =============================================================
       ! DETERMINE REACTION BIN AND ADD VALUE TO SCORE
       type = t % reaction_type
       if (type == TALLY_FLUX) then
          ! Tally flux only
          val = flux
          call add_to_score(t % score(r_bin, c_bin, e_bin), val)

       elseif (type == TALLY_ALL) then
          ! Tally total reaction rate
          val = ONE
          r_bin = 1
          call add_to_score(t % score(r_bin, c_bin, e_bin), val)

       elseif (type == TALLY_BINS) then
          ! Individually bin reactions
          n_reaction = t % reaction_type
          
          r_bin = 0
          do j = 1, n_reaction
             MT = t % reactions(j)
             Sigma = get_macro_xs(p, cMaterial, MT)
             val = Sigma * flux
             r_bin = r_bin + 1
             call add_to_score(t % score(r_bin, c_bin, e_bin), &
                  & val)
          end do
       elseif (type == TALLY_SUM) then
          ! Tally reactions in one bin
          n_reaction = t % reaction_type
          
          r_bin = 1
          do j = 1, n_reaction
             MT = t % reactions(j)
             Sigma = get_macro_xs(p, cMaterial, MT)
             val = Sigma * flux
             call add_to_score(t % score(r_bin, c_bin, e_bin), &
                  & val)
          end do
       end if

    end do

    ! ================================================================
    ! TODO: Add lattice tallies

    ! ================================================================
    ! TODO: Add mesh tallies

  end subroutine score_tally

!=====================================================================
! ADD_TO_SCORE
!=====================================================================

  subroutine add_to_score(score, val)

    type(TallyScore), intent(inout) :: score
    real(8),          intent(in)    :: val
    
    score % n_events = score % n_events + 1
    score % val      = score % val + val
    score % val_sq   = score % val_sq + val*val
    
  end subroutine add_to_score

end module score
