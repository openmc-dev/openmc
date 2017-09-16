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

    use cmfd_header,  only: cmfd_adjoint_type, time_cmfdbuild, time_cmfdsolve

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
    use cmfd_header, only: cmfd_shift, cmfd_ktol, cmfd_stol, cmfd_write_matrices
    use simulation_header, only: keff

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
    use cmfd_header, only: cmfd_write_matrices

    ! Transpose matrices
    loss = loss % transpose()
    prod = prod % transpose()

    ! Write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call loss % write('adj_loss.dat')
      call prod % write('adj_prod.dat')
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXECUTE_POWER_ITER  is the main power iteration routine
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter()

    use constants,  only: ONE
    use error,      only: fatal_error
    use cmfd_header, only: cmfd, cmfd_atoli, cmfd_rtoli

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
      select case(cmfd % indices(4))
      case(1)
        call cmfd_linsolver_1g(loss, s_o, phi_n, toli, innerits)
      case(2)
        call cmfd_linsolver_2g(loss, s_o, phi_n, toli, innerits)
      case default
        call cmfd_linsolver_ng(loss, s_o, phi_n, toli, innerits)
      end select

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

    use, intrinsic :: ISO_FORTRAN_ENV

    use constants,  only: ONE, ZERO
    use cmfd_header, only: cmfd_power_monitor
    use message_passing, only: master

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
! CMFD_LINSOLVER_1g solves the CMFD linear system
!===============================================================================

  subroutine cmfd_linsolver_1g(A, b, x, tol, its)

    use constants,  only: ONE, ZERO
    use error,      only: fatal_error
    use cmfd_header, only: cmfd, cmfd_spectral

    type(Matrix), intent(inout) :: A ! coefficient matrix
    type(Vector), intent(inout) :: b ! right hand side vector
    type(Vector), intent(inout) :: x ! unknown vector
    real(8), intent(in) :: tol ! tolerance on final error
    integer, intent(out) :: its ! number of inner iterations

    integer :: g ! group index
    integer :: i ! loop counter for x
    integer :: j ! loop counter for y
    integer :: k ! loop counter for z
    integer :: n  ! total size of vector
    integer :: nx ! maximum dimension in x direction
    integer :: ny ! maximum dimension in y direction
    integer :: nz ! maximum dimension in z direction
    integer :: ng ! number of energy groups
    integer :: igs   ! Gauss-Seidel iteration counter
    integer :: irb   ! Red/Black iteration switch
    integer :: irow  ! row iteration
    integer :: icol  ! iteration counter over columns
    integer :: didx  ! index for diagonal component
    logical :: found ! did we find col
    real(8) :: tmp1 ! temporary sum g1
    real(8) :: x1 ! new g1 value of x
    real(8) :: err ! error in convergence of solution
    real(8) :: w ! overrelaxation parameter
    type(Vector) :: tmpx ! temporary solution vector

    ! Set overrelaxation parameter
    w = ONE

    ! Dimensions
    ng = 1
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    n = A % n

    ! Perform Gauss Seidel iterations
    GS: do igs = 1, 10000

      ! Check for max iterations met
      if (igs == 10000) then
        call fatal_error('Maximum Gauss-Seidel iterations encountered.')
      endif

      ! Copy over x vector
      call  tmpx % copy(x)

      ! Perform red/black gs iterations
      REDBLACK: do irb = 0,1

        ! Begin loop around matrix rows
        ROWS: do irow = 1, n

          ! Get spatial location
          call  matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

          ! Filter out black cells (even)
          if (mod(i+j+k,2) == irb) cycle

          ! Get the index of the diagonals for both rows
          call A % search_indices(irow, irow, didx, found)

          ! Perform temporary sums, first do left of diag block, then right of diag block
          tmp1 = ZERO
          do icol = A % get_row(irow), didx - 1
            tmp1 = tmp1 + A % val(icol)*x % val(A % get_col(icol))
          end do
          do icol = didx + 1, A % get_row(irow + 1) - 1
            tmp1 = tmp1 + A % val(icol)*x % val(A % get_col(icol))
          end do

          ! Solve for new x
          x1 = (b % val(irow) - tmp1)/A % val(didx)

          ! Perform overrelaxation
          x % val(irow) = (ONE - w)*x % val(irow) + w*x1

        end do ROWS

      end do REDBLACK

      ! Check convergence
      err = sqrt(sum(((tmpx % val - x % val)/tmpx % val)**2)/n)
      its = igs
      if (err < tol) exit

      ! Calculation new overrelaxation parameter
      w = ONE/(ONE - 0.25_8*cmfd_spectral*w)

    end do GS

    call tmpx % destroy()

  end subroutine cmfd_linsolver_1g

!===============================================================================
! CMFD_LINSOLVER_2G solves the CMFD linear system
!===============================================================================

  subroutine cmfd_linsolver_2g(A, b, x, tol, its)

    use constants,  only: ONE, ZERO
    use error,      only: fatal_error
    use cmfd_header, only: cmfd, cmfd_spectral

    type(Matrix), intent(inout) :: A ! coefficient matrix
    type(Vector), intent(inout) :: b ! right hand side vector
    type(Vector), intent(inout) :: x ! unknown vector
    real(8), intent(in) :: tol ! tolerance on final error
    integer, intent(out) :: its ! number of inner iterations

    integer :: g ! group index
    integer :: i ! loop counter for x
    integer :: j ! loop counter for y
    integer :: k ! loop counter for z
    integer :: n  ! total size of vector
    integer :: nx ! maximum dimension in x direction
    integer :: ny ! maximum dimension in y direction
    integer :: nz ! maximum dimension in z direction
    integer :: ng ! number of energy groups
    integer :: d1idx ! index of row "1" diagonal
    integer :: d2idx ! index of row "2" diagonal
    integer :: igs   ! Gauss-Seidel iteration counter
    integer :: irb   ! Red/Black iteration switch
    integer :: irow  ! row iteration
    integer :: icol  ! iteration counter over columns
    logical :: found ! did we find col
    real(8) :: m11 ! block diagonal component 1,1
    real(8) :: m12 ! block diagonal component 1,2
    real(8) :: m21 ! block diagonal component 2,1
    real(8) :: m22 ! block diagonal component 2,2
    real(8) :: dm  ! determinant of block diagonal
    real(8) :: d11 ! inverse component 1,1
    real(8) :: d12 ! inverse component 1,2
    real(8) :: d21 ! inverse component 2,1
    real(8) :: d22 ! inverse component 2,2
    real(8) :: tmp1 ! temporary sum g1
    real(8) :: tmp2 ! temporary sum g2
    real(8) :: x1 ! new g1 value of x
    real(8) :: x2 ! new g2 value of x
    real(8) :: err ! error in convergence of solution
    real(8) :: w ! overrelaxation parameter
    type(Vector) :: tmpx ! temporary solution vector

    ! Set tolerance and overrelaxation parameter
    w = ONE

    ! Dimensions
    ng = 2
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    n = A % n

    ! Perform Gauss Seidel iterations
    GS: do igs = 1, 10000

      ! Check for max iterations met
      if (igs == 10000) then
        call fatal_error('Maximum Gauss-Seidel iterations encountered.')
      endif

      ! Copy over x vector
      call  tmpx % copy(x)

      ! Perform red/black gs iterations
      REDBLACK: do irb = 0,1

        ! Begin loop around matrix rows
        ROWS: do irow = 1, n, 2

          ! Get spatial location
          call  matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

          ! Filter out black cells (even)
          if (mod(i+j+k,2) == irb) cycle

          ! Get the index of the diagonals for both rows
          call A % search_indices(irow, irow, d1idx, found)
          call A % search_indices(irow + 1, irow + 1, d2idx, found)

          ! Get block diagonal
          m11 = A % val(d1idx) ! group 1 diagonal
          m12 = A % val(d1idx + 1) ! group 1 right of diagonal (sorted by col)
          m21 = A % val(d2idx - 1) ! group 2 left of diagonal (sorted by col)
          m22 = A % val(d2idx) ! group 2 diagonal

          ! Analytically invert the diagonal
          dm = m11*m22 - m12*m21
          d11 = m22/dm
          d12 = -m12/dm
          d21 = -m21/dm
          d22 = m11/dm

          ! Perform temporary sums, first do left of diag block, then right of diag block
          tmp1 = ZERO
          tmp2 = ZERO
          do icol = A % get_row(irow), d1idx - 1
            tmp1 = tmp1 + A % val(icol)*x % val(A % get_col(icol))
          end do
          do icol = A % get_row(irow + 1), d2idx - 2
            tmp2 = tmp2 + A % val(icol)*x % val(A % get_col(icol))
          end do
          do icol = d1idx + 2, A % get_row(irow + 1) - 1
            tmp1 = tmp1 + A % val(icol)*x % val(A % get_col(icol))
          end do
          do icol = d2idx + 1, A % get_row(irow + 2) - 1
            tmp2 = tmp2 + A % val(icol)*x % val(A % get_col(icol))
          end do

          ! Adjust with RHS vector
          tmp1 = b % val(irow) - tmp1
          tmp2 = b % val(irow + 1) - tmp2

          ! Solve for new x
          x1 = d11*tmp1 + d12*tmp2
          x2 = d21*tmp1 + d22*tmp2

          ! Perform overrelaxation
          x % val(irow) = (ONE - w)*x % val(irow) + w*x1
          x % val(irow + 1) = (ONE - w)*x % val(irow + 1) + w*x2

        end do ROWS

      end do REDBLACK

      ! Check convergence
      err = sqrt(sum(((tmpx % val - x % val)/tmpx % val)**2)/n)
      its = igs
      if (err < tol) exit

      ! Calculation new overrelaxation parameter
      w = ONE/(ONE - 0.25_8*cmfd_spectral*w)

    end do GS

    call tmpx % destroy()

  end subroutine cmfd_linsolver_2g

!===============================================================================
! CMFD_LINSOLVER_ng solves the CMFD linear system
!===============================================================================

  subroutine cmfd_linsolver_ng(A, b, x, tol, its)

    use constants,  only: ONE, ZERO
    use error,      only: fatal_error
    use cmfd_header, only: cmfd, cmfd_spectral

    type(Matrix), intent(inout) :: A ! coefficient matrix
    type(Vector), intent(inout) :: b ! right hand side vector
    type(Vector), intent(inout) :: x ! unknown vector
    real(8), intent(in) :: tol ! tolerance on final error
    integer, intent(out) :: its ! number of inner iterations

    integer :: g ! group index
    integer :: i ! loop counter for x
    integer :: j ! loop counter for y
    integer :: k ! loop counter for z
    integer :: n  ! total size of vector
    integer :: nx ! maximum dimension in x direction
    integer :: ny ! maximum dimension in y direction
    integer :: nz ! maximum dimension in z direction
    integer :: ng ! number of energy groups
    integer :: igs   ! Gauss-Seidel iteration counter
    integer :: irow  ! row iteration
    integer :: icol  ! iteration counter over columns
    integer :: didx  ! index for diagonal component
    logical :: found ! did we find col
    real(8) :: tmp1 ! temporary sum g1
    real(8) :: x1 ! new g1 value of x
    real(8) :: err ! error in convergence of solution
    real(8) :: w ! overrelaxation parameter
    type(Vector) :: tmpx ! temporary solution vector

    ! Set overrelaxation parameter
    w = ONE

    ! Dimensions
    ng = 1
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    n = A % n

    ! Perform Gauss Seidel iterations
    GS: do igs = 1, 10000

      ! Check for max iterations met
      if (igs == 10000) then
        call fatal_error('Maximum Gauss-Seidel iterations encountered.')
      endif

      ! Copy over x vector
      call  tmpx % copy(x)

      ! Begin loop around matrix rows
      ROWS: do irow = 1, n

        ! Get spatial location
        call  matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

        ! Get the index of the diagonals for both rows
        call A % search_indices(irow, irow, didx, found)

        ! Perform temporary sums, first do left of diag block, then right of diag block
        tmp1 = ZERO
        do icol = A % get_row(irow), didx - 1
          tmp1 = tmp1 + A % val(icol)*x % val(A % get_col(icol))
        end do
        do icol = didx + 1, A % get_row(irow + 1) - 1
          tmp1 = tmp1 + A % val(icol)*x % val(A % get_col(icol))
        end do

        ! Solve for new x
        x1 = (b % val(irow) - tmp1)/A % val(didx)

        ! Perform overrelaxation
        x % val(irow) = (ONE - w)*x % val(irow) + w*x1

      end do ROWS

      ! Check convergence
      err = sqrt(sum(((tmpx % val - x % val)/tmpx % val)**2)/n)
      its = igs

      if (err < tol) exit

      ! Calculation new overrelaxation parameter
      w = ONE/(ONE - 0.25_8*cmfd_spectral*w)

    end do GS

    call tmpx % destroy()

  end subroutine cmfd_linsolver_ng

!===============================================================================
! EXTRACT_RESULTS takes results and puts them in CMFD global data object
!===============================================================================

  subroutine extract_results()

    use cmfd_header, only: cmfd, cmfd_write_matrices
    use simulation_header, only: current_batch

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
        filename = 'adj_fluxvec.dat'
      else
        filename = 'fluxvec.dat'
      end if
      ! TODO: call phi_n % write(filename)
    end if

  end subroutine extract_results

!===============================================================================
! MATRIX_TO_INDICES converts a matrix index to spatial and group indicies
!===============================================================================

  subroutine matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

    use cmfd_header, only: cmfd, cmfd_coremap

    integer, intent(out) :: i    ! iteration counter for x
    integer, intent(out) :: j    ! iteration counter for y
    integer, intent(out) :: k    ! iteration counter for z
    integer, intent(out) :: g    ! iteration counter for groups
    integer, intent(in)  :: irow ! iteration counter over row (0 reference)
    integer, intent(in)  :: nx   ! maximum number of x cells
    integer, intent(in)  :: ny   ! maximum number of y cells
    integer, intent(in)  :: nz   ! maximum number of z cells
    integer, intent(in)  :: ng   ! maximum number of groups

    ! Check for core map
    if (cmfd_coremap) then

      ! Get indices from indexmap
      g = mod(irow-1, ng) + 1
      i = cmfd % indexmap((irow-1)/ng+1,1)
      j = cmfd % indexmap((irow-1)/ng+1,2)
      k = cmfd % indexmap((irow-1)/ng+1,3)

    else

      ! Compute indices
      g = mod(irow-1, ng) + 1
      i = mod(irow-1, ng*nx)/ng + 1
      j = mod(irow-1, ng*nx*ny)/(ng*nx)+ 1
      k = mod(irow-1, ng*nx*ny*nz)/(ng*nx*ny) + 1

    end if

  end subroutine matrix_to_indices

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
