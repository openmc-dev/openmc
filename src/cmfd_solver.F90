module cmfd_solver

! This module contains routines to execute the power iteration solver

  use cmfd_loss_operator, only: init_loss_matrix, build_loss_matrix
  use cmfd_prod_operator, only: init_prod_matrix, build_prod_matrix
  use matrix_header,      only: Matrix
  use vector_header,      only: Vector

  implicit none
  private
  public :: cmfd_solver_execute 

  logical :: iconv                ! did the problem converged
  real(8) :: k_n                  ! new k-eigenvalue
  real(8) :: k_o                  ! old k-eigenvalue
  real(8) :: ktol = 1.e-8_8       ! tolerance on keff
  real(8) :: stol = 1.e-8_8       ! tolerance on source
  real(8) :: norm_n               ! current norm of source vector
  real(8) :: norm_o               ! old norm of source vector
  real(8) :: kerr                 ! error in keff
  real(8) :: serr                 ! error in source
  logical :: adjoint_calc         ! run an adjoint calculation
  type(Matrix) :: loss            ! cmfd loss matrix
  type(Matrix) :: prod            ! cmfd prod matrix
  type(Vector) :: phi_n           ! new flux vector
  type(Vector) :: phi_o           ! old flux vector
  type(Vector) :: s_n             ! new source vector
  type(Vector) :: s_o             ! old flux vector
  type(Vector) :: serr_v          ! error in source

contains

!===============================================================================
! CMFD_SOLVER_EXECUTE sets up and runs power iteration solver for CMFD
!===============================================================================

  subroutine cmfd_solver_execute(k_tol, s_tol, adjoint)

    use global,  only: cmfd_adjoint_type, time_cmfdbuild, time_cmfdsolve

    real(8), optional :: k_tol    ! tolerance on keff
    real(8), optional :: s_tol    ! tolerance on source
    logical, optional :: adjoint  ! adjoint calc

    logical :: physical_adjoint = .false.

    ! Set tolerances if present
    if (present(k_tol)) ktol = k_tol
    if (present(s_tol)) stol = s_tol

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
    use global,    only: cmfd_write_matrices

    logical :: adjoint

    integer :: n      ! problem size
    real(8) :: guess  ! initial guess

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
    k_n = guess
    k_o = guess

    ! Fill in loss matrix
    call build_loss_matrix(loss, adjoint=adjoint) 

    ! Fill in production matrix
    call build_prod_matrix(prod, adjoint=adjoint)

    ! Setup petsc for everything
    call loss % assemble()
    call prod % assemble()
    call loss % write('loss.dat')
    call prod % write('prod.dat')

    ! Set norms to 0
    norm_n = ZERO
    norm_o = ZERO

  end subroutine init_data

!===============================================================================
! COMPUTE_ADJOINT computes a mathematical adjoint of CMFD problem 
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    ! Transpose matrices
    call loss % transpose()
    call prod % transpose()

    ! Write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call loss % write_petsc_binary('adj_lossmat.bin')
      call prod % write_petsc_binary('adj_prodmat.bin')
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXECUTE_POWER_ITER  is the main power iteration routine
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter()

    integer :: i ! iteration counter

    ! Reset convergence flag
    iconv = .false.

    ! Begin power iteration
    do i = 1, 10000

      ! Compute source vector
      call prod % vector_multiply(phi_o, s_o)

      ! Normalize source vector
      s_o % val = s_o % val / k_o

      ! Compute new flux vector
      call cmfd_linsolver(s_o, phi_n)

      ! Compute new source vector
      call prod % vector_multiply(phi_n, s_n)

      ! Compute new k-eigenvalue
      k_n = sum(s_n % val) / sum(s_o % val)

      ! Renormalize the old source
      s_o % val = s_o % val * k_o

      ! Check convergence
      call convergence(i)

      ! Break loop if converged
      if (iconv) exit

      ! Record old values
      phi_o % val = phi_n % val
      k_o = k_n
      norm_o = norm_n

    end do

  end subroutine execute_power_iter 

!===============================================================================
! CONVERGENCE checks the convergence of the CMFD problem
!===============================================================================

  subroutine convergence(iter)

    use constants,  only: ONE, TINY_BIT
    use global,     only: cmfd_power_monitor, master
    use, intrinsic :: ISO_FORTRAN_ENV

    integer     :: iter           ! iteration number

    ! Reset convergence flag
    iconv = .false.

    ! Calculate error in keff
    kerr = abs(k_o - k_n)/k_n

    ! Calculate max error in source
    where (s_n % val > TINY_BIT)
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
           &1PE12.5,T55, "src-error: ",1PE12.5)') iter, k_n, kerr, serr
    end if

  end subroutine convergence

!===============================================================================
! CMFD_LINSOLVER solves the CMFD linear system
!===============================================================================

  subroutine cmfd_linsolver(b, x)

    use constants,  only: ONE, ZERO
    use global,     only: cmfd

    type(Vector) :: b ! right hand side vector
    type(Vector) :: x ! unknown vector

    integer :: i ! loop counter for x 
    integer :: j ! loop counter for y
    integer :: k ! loop counter for z
    integer :: n  ! total size of vector
    integer :: nx ! maximum dimension in x direction
    integer :: ny ! maximum dimension in y direction
    integer :: nz ! maximum dimension in z direction
    integer :: ng ! number of energy groups
    integer :: matidx ! matrix index of row
    integer :: d1idx ! index of row "1" diagonal
    integer :: d2idx ! index of row "2" diagonal
    integer :: igs   ! Gauss-Seidel iteration counter
    integer :: irb   ! Red/Black iteration switch
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
    real(8) :: err ! error in convergence of solution
    real(8) :: tol ! tolerance on final error
    type(Vector) :: tmpx ! temporary solution vector

    ! Set tolerance
    tol = 1.e-10_8

    ! Dimensions
    ng = 2
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    n = ng*nx*ny*nz

    ! Perform Gauss Seidel iterations
    GS: do igs = 1, 10000

      ! Copy over x vector
      call  tmpx % copy(x) 

      ! Perform red/black gs iterations
      REDBLACK: do irb = 0,1

        ZLOOP: do k = 1, nz

          YLOOP: do j = 1, ny

            XLOOP: do i = 1, nx

              ! Filter out black cells (even)
              if (mod(i+j+k,2) == irb) cycle

              ! Get starting row in matrix for this block
              call indices_to_matrix(1, i, j, k, matidx, ng, nx, ny)

              ! Get the index of the diagonals for both rows
              call loss % search_indices(matidx, matidx, d1idx, found)

              call loss % search_indices(matidx + 1, matidx + 1, d2idx, found)

              ! Get block diagonal
              m11 = loss % val(d1idx) ! group 1 diagonal
              m12 = loss % val(d1idx + 1) ! group 1 right of diagonal (sorted by col)
              m21 = loss % val(d2idx - 1) ! group 2 left of diagonal (sorted by col)
              m22 = loss % val(d2idx) ! group 2 diagonal

              ! Analytically invert the diagonal
              dm = m11*m22 - m12*m21
              d11 = m22/dm
              d12 = -m12/dm
              d21 = -m21/dm
              d22 = m11/dm

              ! Perform temporary sums, first do left of diag block, then right of diag block
              tmp1 = ZERO
              tmp2 = ZERO
              do icol = loss % get_row(matidx), d1idx - 1
                tmp1 = tmp1 + loss % val(icol)*x % val(loss % get_col(icol))
              end do
              do icol = loss % get_row(matidx + 1), d2idx - 2
                tmp2 = tmp2 + loss % val(icol)*x % val(loss % get_col(icol))
              end do
              do icol = d1idx + 2, loss % get_row(matidx + 1) - 1
                tmp1 = tmp1 + loss % val(icol)*x % val(loss % get_col(icol))
              end do
              do icol = d2idx + 1, loss % get_row(matidx + 2) - 1
                tmp2 = tmp2 + loss % val(icol)*x % val(loss % get_col(icol))
              end do

              ! Adjust with RHS vector
              tmp1 = b % val(matidx) - tmp1
              tmp2 = b % val(matidx + 1) - tmp2

              ! Solve for new x
              x % val(matidx) = d11*tmp1 + d12*tmp2
              x % val(matidx + 1) = d21*tmp1 + d22*tmp2

            end do XLOOP

          end do YLOOP

        end do ZLOOP

      end do REDBLACK

      ! Check convergence
      err = sqrt(sum(((tmpx % val - x % val)/tmpx % val)**2)/n)
      if (err < tol) exit

    end do GS

    call tmpx % destroy()

  end subroutine cmfd_linsolver

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
      call phi_n % write_petsc_binary(filename)
    end if

  end subroutine extract_results

!===============================================================================
! INDICES_TO_MATRIX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  subroutine indices_to_matrix(g, i, j, k, matidx, ng, nx, ny)

    use global,  only: cmfd, cmfd_coremap

    integer :: matidx ! the index location in matrix
    integer :: i      ! current x index
    integer :: j      ! current y index
    integer :: k      ! current z index
    integer :: g      ! current group index
    integer :: nx     ! maximum number of x cells
    integer :: ny     ! maximum number of y cells
    integer :: ng     ! maximum number of groups

    ! Check if coremap is used
    if (cmfd_coremap) then

      ! Get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! Compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

  end subroutine indices_to_matrix

!===============================================================================
! FINALIZE frees all memory associated with power iteration
!===============================================================================

  subroutine finalize()

    ! Destroy all objects 
#ifdef PETSC
    call gmres  % destroy()
#endif
    call loss   % destroy() 
    call prod   % destroy()
    call phi_n  % destroy()
    call phi_o  % destroy()
    call s_n    % destroy()
    call s_o    % destroy()
    call serr_v % destroy

  end subroutine finalize

end module cmfd_solver
