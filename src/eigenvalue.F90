module eigenvalue

  use, intrinsic :: ISO_C_BINDING

  use constants,   only: ZERO
  use message_passing
  use settings
  use simulation_header
  use tally_header

  implicit none

  interface
    subroutine calculate_generation_keff() bind(C)
    end subroutine

    subroutine calculate_average_keff() bind(C)
    end subroutine
  end interface

contains

!===============================================================================
! OPENMC_GET_KEFF calculates a minimum variance estimate of k-effective based on
! a linear combination of the collision, absorption, and tracklength
! estimates. The theory behind this can be found in M. Halperin, "Almost
! linearly-optimum combination of unbiased estimates," J. Am. Stat. Assoc., 56,
! 36-43 (1961), doi:10.1080/01621459.1961.10482088. The implementation here
! follows that described in T. Urbatsch et al., "Estimation and interpretation
! of keff confidence intervals in MCNP," Nucl. Technol., 111, 169-182 (1995).
!===============================================================================

  function openmc_get_keff(k_combined) result(err) bind(C)
    real(C_DOUBLE), intent(out) :: k_combined(2)
    integer(C_INT) :: err

    integer :: l        ! loop index
    integer :: i, j, k  ! indices referring to collision, absorption, or track
    real(8) :: n        ! number of realizations
    real(8) :: kv(3)    ! vector of k-effective estimates
    real(8) :: cov(3,3) ! sample covariance matrix
    real(8) :: f        ! weighting factor
    real(8) :: g        ! sum of weighting factors
    real(8) :: S(3)     ! sums used for variance calculation

    k_combined = ZERO

    ! Make sure we have at least four realizations. Notice that at the end,
    ! there is a N-3 term in a denominator.
    if (n_realizations <= 3) then
      err = -1
      return
    end if

    ! Initialize variables
    n = real(n_realizations, 8)

    ! Copy estimates of k-effective and its variance (not variance of the mean)
    kv(1) = global_tallies(RESULT_SUM, K_COLLISION) / n
    kv(2) = global_tallies(RESULT_SUM, K_ABSORPTION) / n
    kv(3) = global_tallies(RESULT_SUM, K_TRACKLENGTH) / n
    cov(1, 1) = (global_tallies(RESULT_SUM_SQ, K_COLLISION) - &
         n * kv(1) * kv(1)) / (n - ONE)
    cov(2, 2) = (global_tallies(RESULT_SUM_SQ, K_ABSORPTION) - &
         n * kv(2) * kv(2)) / (n - ONE)
    cov(3, 3) = (global_tallies(RESULT_SUM_SQ, K_TRACKLENGTH) - &
         n * kv(3) * kv(3)) / (n - ONE)

    ! Calculate covariances based on sums with Bessel's correction
    cov(1, 2) = (k_col_abs - n * kv(1) * kv(2)) / (n - ONE)
    cov(1, 3) = (k_col_tra - n * kv(1) * kv(3)) / (n - ONE)
    cov(2, 3) = (k_abs_tra - n * kv(2) * kv(3)) / (n - ONE)
    cov(2, 1) = cov(1, 2)
    cov(3, 1) = cov(1, 3)
    cov(3, 2) = cov(2, 3)

    ! Check to see if two estimators are the same; this is guaranteed to happen
    ! in MG-mode with survival biasing when the collision and absorption
    ! estimators are the same, but can theoretically happen at anytime.
    ! If it does, the standard estimators will produce floating-point
    ! exceptions and an expression specifically derived for the combination of
    ! two estimators (vice three) should be used instead.

    ! First we will identify if there are any matching estimators
    if ((abs(kv(1) - kv(2)) / kv(1) < FP_REL_PRECISION) .and. &
         (abs(cov(1, 1) - cov(2, 2)) / cov(1, 1) < FP_REL_PRECISION)) then
      ! 1 and 2 match, so only use 1 and 3 in our comparisons
      i = 1
      j = 3

    else if ((abs(kv(1) - kv(3)) / kv(1) < FP_REL_PRECISION) .and. &
         (abs(cov(1, 1) - cov(3, 3)) / cov(1, 1) < FP_REL_PRECISION)) then
      ! 1 and 3 match, so only use 1 and 2 in our comparisons
      i = 1
      j = 2

    else if ((abs(kv(2) - kv(3)) / kv(2) < FP_REL_PRECISION) .and. &
         (abs(cov(2, 2) - cov(3, 3)) / cov(2, 2) < FP_REL_PRECISION)) then
      ! 2 and 3 match, so only use 1 and 2 in our comparisons
      i = 1
      j = 2

    else
      ! No two estimators match, so set i to 0 and this will be the indicator
      ! to use all three estimators.
      i = 0
    end if

    if (i == 0) then
      ! Use three estimators as derived in the paper by Urbatsch

      ! Initialize variables
      g = ZERO
      S = ZERO

      do l = 1, 3
        ! Permutations of estimates
        if (l == 1) then
          ! i = collision, j = absorption, k = tracklength
          i = 1
          j = 2
          k = 3
        elseif (l == 2) then
          ! i = absortion, j = tracklength, k = collision
          i = 2
          j = 3
          k = 1
        elseif (l == 3) then
          ! i = tracklength, j = collision, k = absorption
          i = 3
          j = 1
          k = 2
        end if

        ! Calculate weighting
        f = cov(j, j) * (cov(k, k) - cov(i, k)) - cov(k, k) * cov(i, j) + &
             cov(j, k) * (cov(i, j) + cov(i, k) - cov(j, k))

        ! Add to S sums for variance of combined estimate
        S(1) = S(1) + f * cov(1, l)
        S(2) = S(2) + (cov(j, j) + cov(k, k) - TWO * cov(j, k)) * kv(l) * kv(l)
        S(3) = S(3) + (cov(k, k) + cov(i, j) - cov(j, k) - &
             cov(i, k)) * kv(l) * kv(j)

        ! Add to sum for combined k-effective
        k_combined(1) = k_combined(1) + f * kv(l)
        g = g + f
      end do

      ! Complete calculations of S sums
      S = (n - ONE) * S
      S(1) = (n - ONE)**2 * S(1)

      ! Calculate combined estimate of k-effective
      k_combined(1) = k_combined(1) / g

      ! Calculate standard deviation of combined estimate
      g = (n - ONE)**2 * g
      k_combined(2) = sqrt(S(1) / &
           (g * n * (n - THREE)) * (ONE + n * ((S(2) - TWO * S(3)) / g)))

    else
      ! Use only two estimators
      ! These equations are derived analogously to that done in the paper by
      ! Urbatsch, but are simpler than for the three estimators case since the
      ! block matrices of the three estimator equations reduces to scalars here

      ! Store the commonly used term
      f = kv(i) - kv(j)
      g = cov(i, i) + cov(j, j) - TWO * cov(i, j)

      ! Calculate combined estimate of k-effective
      k_combined(1) = kv(i) - (cov(i, i) - cov(i, j)) / g * f

      ! Calculate standard deviation of combined estimate
      k_combined(2) = (cov(i, i) * cov(j, j) - cov(i, j) * cov(i, j)) * &
           (g + n * f * f) / (n * (n - TWO) * g * g)
      k_combined(2) = sqrt(k_combined(2))

    end if
    err = 0

  end function openmc_get_keff

#ifdef _OPENMP
!===============================================================================
! JOIN_BANK_FROM_THREADS joins threadprivate fission banks into a single fission
! bank that can be sampled. Note that this operation is necessarily sequential
! to preserve the order of the bank when using varying numbers of threads.
!===============================================================================

  subroutine join_bank_from_threads()

    integer(8) :: total ! total number of fission bank sites
    integer    :: i     ! loop index for threads

    ! Initialize the total number of fission bank sites
    total = 0

!$omp parallel

    ! Copy thread fission bank sites to one shared copy
!$omp do ordered schedule(static)
    do i = 1, n_threads
!$omp ordered
      master_fission_bank(total+1:total+n_bank) = fission_bank(1:n_bank)
      total = total + n_bank
!$omp end ordered
    end do
!$omp end do

    ! Make sure all threads have made it to this point
!$omp barrier

    ! Now copy the shared fission bank sites back to the master thread's copy.
    if (thread_id == 0) then
      n_bank = total
      fission_bank(1:n_bank) = master_fission_bank(1:n_bank)
    else
      n_bank = 0
    end if

!$omp end parallel

  end subroutine join_bank_from_threads
#endif

end module eigenvalue
