module cmfd_solver

! This module contains routines to execute the power iteration solver

  use constants,          only: MAX_LINE_LEN
  use cmfd_loss_operator, only: init_loss_matrix, build_loss_matrix
  use cmfd_prod_operator, only: init_prod_matrix, build_prod_matrix
  use matrix_header,      only: Matrix
  use vector_header,      only: Vector

  implicit none
  private
  public :: cmfd_solver_execute

  real(8) :: k_n          ! New k-eigenvalue
  real(8) :: k_o          ! Old k-eigenvalue
  real(8) :: k_s          ! Shift of eigenvalue
  real(8) :: k_ln         ! New shifted eigenvalue
  real(8) :: k_lo         ! Old shifted eigenvalue
  real(8) :: norm_n       ! Current norm of source vector
  real(8) :: norm_o       ! Old norm of source vector
  real(8) :: kerr         ! Error in keff
  real(8) :: serr         ! Error in source
  real(8) :: ktol         ! Tolerance on keff
  real(8) :: stol         ! Tolerance on source
  logical :: adjoint_calc ! Run an adjoint calculation
  type(Matrix) :: loss    ! Cmfd loss matrix
  type(Matrix) :: prod    ! Cmfd prod matrix
  type(Vector) :: phi_n   ! New flux vector
  type(Vector) :: phi_o   ! Old flux vector
  type(Vector) :: s_n     ! New source vector
  type(Vector) :: s_o     ! Old flux vector
  type(Vector) :: serr_v  ! Error in source

  ! CMFD linear solver interface
  procedure(linsolve), pointer :: cmfd_linsolver => null()
  abstract interface
    subroutine linsolve(A, b, x, tol, i)
      import :: Matrix
      import :: Vector
      type(Matrix), intent(inout) :: A
      type(Vector), intent(inout) :: b
      type(Vector), intent(inout) :: x
      real(8), intent(in) :: tol
      integer, intent(out) :: i
    end subroutine linsolve
  end interface

contains

!===============================================================================
! CMFD_SOLVER_EXECUTE sets up and runs power iteration solver for CMFD
!===============================================================================

  subroutine cmfd_solver_execute(adjoint)

    use global,  only: cmfd_adjoint_type, time_cmfdbuild, time_cmfdsolve

    logical, optional, intent(in) :: adjoint  ! adjoint calc

    logical :: physical_adjoint = .false.

    ! Check for adjoint execution
    adjoint_calc = .false.
    if (present(adjoint)) adjoint_calc = adjoint

    ! Check for physical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical') &
        physical_adjoint = .true.

    ! Start timer for build
    call time_cmfdbuild % start()

    ! Initialize matrices and vectors
    call init_data(physical_adjoint)

    ! Check for mathematical adjoint calculation
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math') &
        call compute_adjoint()

    ! Stop timer for build
    call time_cmfdbuild % stop()

    ! Begin power iteration
    call time_cmfdsolve % start()
    call execute_power_iter()
    call time_cmfdsolve % stop()

    ! Extract results
    call extract_results()

    ! Deallocate data
    call finalize()

  end subroutine cmfd_solver_execute

!===============================================================================
! INIT_DATA allocates matrices and vectors for CMFD solution
!===============================================================================

  subroutine init_data(adjoint)

    use constants, only: ONE, ZERO
    use global,    only: cmfd_shift, keff, cmfd_ktol, cmfd_stol, &
                         cmfd_write_matrices

    logical, intent(in) :: adjoint

    integer :: n      ! problem size
    real(8) :: guess  ! initial guess
    real(8) :: dw     ! eigenvalue shift

    ! Set up matrices
    call init_loss_matrix(loss)
    call init_prod_matrix(prod)

    ! Get problem size
    n = loss % n

    ! Set up flux vectors
    call phi_n % create(n)
    call phi_o % create(n)

    ! Set up source vectors
    call s_n % create(n)
    call s_o % create(n)
    call serr_v % create(n)

    ! Set initial guess
    guess = ONE
    phi_n % val = guess
    phi_o % val = guess
    k_n = keff
    k_o = k_n
    dw = cmfd_shift
    k_s = k_o + dw
    k_ln = ONE/(ONE/k_n - ONE/k_s)
    k_lo = k_ln

    ! Fill in loss matrix
    call build_loss_matrix(loss, adjoint=adjoint)

    ! Fill in production matrix
    call build_prod_matrix(prod, adjoint=adjoint)

    ! Finalize setup of CSR matrices
    call loss % assemble()
    call prod % assemble()
    if (cmfd_write_matrices) then
      call loss % write('loss.dat')
      call prod % write('prod.dat')
    end if

    ! Set norms to 0
    norm_n = ZERO
    norm_o = ZERO

    ! Set tolerances
    ktol = cmfd_ktol
    stol = cmfd_stol

  end subroutine init_data

!===============================================================================
! COMPUTE_ADJOINT computes a mathematical adjoint of CMFD problem
!===============================================================================

  subroutine compute_adjoint()

    use error,   only: fatal_error
#ifdef PETSC
    use global,  only: cmfd_write_matrices
#endif

#ifdef PETSC
    ! Transpose matrices
    call loss % transpose()
    call prod % transpose()

    ! Write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call loss % write_petsc_binary('adj_lossmat.bin')
      call prod % write_petsc_binary('adj_prodmat.bin')
    end if
#else
    call fatal_error('Adjoint calculations only allowed with PETSc')
#endif

  end subroutine compute_adjoint

!===============================================================================
! EXECUTE_POWER_ITER  is the main power iteration routine
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter()

    use constants,  only: ONE
    use error,      only: fatal_error
    use global,     only: cmfd_atoli, cmfd_rtoli

    integer :: i ! iteration counter
    integer :: innerits ! # of inner iterations
    integer :: totalits ! total number of inners
    logical :: iconv ! did the problem converged
    real(8) :: atoli ! absolute minimum tolerance
    real(8) :: rtoli ! relative tolerance based on source conv
    real(8) :: toli ! the current tolerance of inners

    ! Reset convergence flag
    iconv = .false.

    ! Set up tolerances
    atoli = cmfd_atoli
    rtoli = cmfd_rtoli
    toli = rtoli*100._8

    ! Perform shift
    call wielandt_shift()
    totalits = 0

    ! Begin power iteration
    do i = 1, 10000

      ! Check if reached iteration 10000
      if (i == 10000) then
        call fatal_error('Reached maximum iterations in CMFD power iteration &
             &solver.')
      end if

      ! Compute source vector
      call prod % vector_multiply(phi_o, s_o)

      ! Normalize source vector
      s_o % val = s_o % val / k_lo

      ! Compute new flux vector
      call cmfd_linsolver(loss, s_o, phi_n, toli, innerits)

      ! Compute new source vector
      call prod % vector_multiply(phi_n, s_n)

      ! Compute new shifted eigenvalue
      k_ln = sum(s_n % val) / sum(s_o % val)

      ! Compute new eigenvalue
      k_n = ONE/(ONE/k_ln + ONE/k_s)

      ! Renormalize the old source
      s_o % val = s_o % val * k_lo

      ! Check convergence
      call convergence(i, innerits, iconv)
      totalits = totalits + innerits

      ! Break loop if converged
      if (iconv) exit

      ! Record old values
      phi_o % val = phi_n % val
      k_o = k_n
      k_lo = k_ln
      norm_o = norm_n

      ! Get new tolerance for inners
      toli = max(atoli, rtoli*serr)

    end do

  end subroutine execute_power_iter

!===============================================================================
! WIELANDT SHIFT
!===============================================================================

  subroutine wielandt_shift()

    use constants,  only: ONE

    integer :: irow ! row counter
    integer :: icol ! col counter
    integer :: jcol ! current col index in prod matrix

    ! perform subtraction
    jcol = 1
    ROWS: do irow = 1, loss % n
      COLS: do icol = loss % get_row(irow), loss % get_row(irow + 1) - 1
        if (loss % get_col(icol) == prod % get_col(jcol) .and. &
            jcol < prod % get_row(irow + 1)) then
          loss % val(icol) = loss % val(icol) - ONE/k_s*prod % val(jcol)
          jcol = jcol + 1
        end if
      end do COLS
    end do ROWS

  end subroutine wielandt_shift

!===============================================================================
! CONVERGENCE checks the convergence of the CMFD problem
!===============================================================================

  subroutine convergence(iter, innerits, iconv)

    use constants,  only: ONE, ZERO
    use global,     only: cmfd_power_monitor, master
    use, intrinsic :: ISO_FORTRAN_ENV

    integer, intent(in) :: iter     ! outer iteration number
    integer, intent(in) :: innerits ! inner iteration nubmer
    logical, intent(out) :: iconv   ! convergence logical

    ! Reset convergence flag
    iconv = .false.

    ! Calculate error in keff
    kerr = abs(k_o - k_n)/k_n

    ! Calculate max error in source
    where (s_n % val > ZERO)
      serr_v % val = ((s_n % val - s_o % val)/s_n % val)**2
    end where
    serr = sqrt(ONE/dble(s_n % n) * sum(serr_v % val))

    ! Check for convergence
    if(kerr < ktol .and. serr < stol) iconv = .true.

    ! Save the L2 norm of the source
    norm_n = serr

    ! Print out to user
    if (cmfd_power_monitor .and. master) then
      write(OUTPUT_UNIT,FMT='(I0,":",T10,"k-eff: ",F0.8,T30,"k-error: ", &
           &1PE12.5,T55, "src-error: ",1PE12.5,T80,I0)') iter, k_n, kerr, &
            serr, innerits
    end if

  end subroutine convergence

!===============================================================================
! EXTRACT_RESULTS takes results and puts them in CMFD global data object
!===============================================================================

  subroutine extract_results()

    use global, only: cmfd, cmfd_write_matrices, current_batch

    character(len=25)    :: filename  ! name of file to write data
    integer              :: n         ! problem size

    ! Get problem size
    n = loss % n

    ! Allocate in cmfd object if not already allocated
    if (adjoint_calc) then
      if (.not. allocated(cmfd%adj_phi)) allocate(cmfd%adj_phi(n))
    else
      if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))
    end if

    ! Save values
    if (adjoint_calc) then
      cmfd % adj_phi = phi_n % val
    else
      cmfd % phi = phi_n % val
    end if

    ! Save eigenvalue
    if(adjoint_calc) then
      cmfd%adj_keff = k_n
    else
      cmfd%keff = k_n
    end if

    ! Normalize phi to 1
    if (adjoint_calc) then
      cmfd%adj_phi = cmfd%adj_phi/sqrt(sum(cmfd%adj_phi*cmfd%adj_phi))
    else
      cmfd%phi = cmfd%phi/sqrt(sum(cmfd%phi*cmfd%phi))
    end if

    ! Save dominance ratio
    cmfd % dom(current_batch) = norm_n/norm_o

    ! Write out results
    if (cmfd_write_matrices) then
      if (adjoint_calc) then
        filename = 'adj_fluxvec.bin'
      else
        filename = 'fluxvec.bin'
      end if
#ifdef PETSC
      call phi_n % write_petsc_binary(filename)
#endif
    end if

  end subroutine extract_results

!===============================================================================
! FINALIZE frees all memory associated with power iteration
!===============================================================================

  subroutine finalize()

    ! Destroy all objects
    call loss   % destroy()
    call prod   % destroy()
    call phi_n  % destroy()
    call phi_o  % destroy()
    call s_n    % destroy()
    call s_o    % destroy()
    call serr_v % destroy

  end subroutine finalize

end module cmfd_solver
