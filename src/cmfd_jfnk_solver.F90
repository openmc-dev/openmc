module cmfd_jfnk_solver

  use cmfd_loss_operator,  only: init_loss_matrix, build_loss_matrix
  use cmfd_power_solver,   only: cmfd_power_execute
  use cmfd_prod_operator,  only: init_prod_matrix, build_prod_matrix
  use matrix_header,       only: Matrix
  use petsc_solver,        only: Petsc_jfnk, Jfnk_ctx
  use vector_header,       only: Vector

  implicit none
  private
  public :: cmfd_jfnk_execute

  logical          :: adjoint_calc
  type(Jfnk_ctx)   :: jfnk_data
  type(Matrix), target     :: jac_prec
  type(Matrix)     :: loss
  type(Matrix)     :: prod
  type(Petsc_jfnk) :: jfnk 
  type(Vector), target     :: resvec
  type(Vector), target     :: xvec

contains

!===============================================================================
! CMFD_JFNK_EXECUTE
!===============================================================================

  subroutine cmfd_jfnk_execute(adjoint)

    logical, optional :: adjoint ! adjoint calculation

    ! check for adjoint
    adjoint_calc = .false.
    if (present(adjoint)) adjoint_calc = adjoint

    ! seed with power iteration
    call cmfd_power_execute(k_tol=1.E-3_8, s_tol=1.E-3_8, adjoint=adjoint_calc)

    ! initialize data
    call init_data()

    ! initialize solver
    call jfnk % create()

    ! set up residual and jacobian routines
    jfnk_data % res_proc_ptr => compute_nonlinear_residual
    jfnk_data % jac_proc_ptr => build_jacobian_matrix
    call jfnk % set_functions(jfnk_data, resvec, jac_prec)

    ! solve the system
    call jfnk % solve(xvec)

    ! extracts results to cmfd object
    call extract_results()

    ! deallocate all slepc data
    call finalize()

  end subroutine cmfd_jfnk_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

    use constants, only: ONE
    use global,    only: cmfd, cmfd_adjoint_type

    logical :: physical_adjoint
    integer :: n
    type(Vector), pointer :: x
    type(Vector), pointer :: res
    type(Matrix), pointer :: jac
    x => xvec
    res => resvec
    jac => jac_prec
    ! Set up all matrices
    call init_loss_matrix(loss)
    call init_prod_matrix(prod)
    call init_jacobian_matrix()

    ! Check for physical adjoint
    physical_adjoint = .false.
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical') &
         physical_adjoint = .true.

    ! Create matrix operators
    call build_loss_matrix(loss, adjoint = physical_adjoint)
    call build_prod_matrix(prod, adjoint = physical_adjoint)

    ! Assembly matrices and use petsc
    call loss % assemble()
    call prod % assemble()
    call loss % setup_petsc()
    call prod % setup_petsc()

    ! Check for mathematical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math') &
         call compute_adjoint()

    ! Get problem size
    n = loss % n

    ! Create problem vectors
    call resvec % create(n + 1)
    call xvec % create(n + 1)

    ! set flux in guess
    if (adjoint_calc) then
      xvec % val(1:n) = cmfd % adj_phi
    else
      xvec % val(1:n) = cmfd % phi
    end if

    ! set keff in guess
    if (adjoint_calc) then
      xvec % val(n + 1) = ONE/cmfd % adj_keff
    else
      xvec % val(n + 1) = ONE/cmfd % keff
    end if

    ! Set up vectors for petsc
    call resvec   % setup_petsc()
    call xvec     % setup_petsc()

    ! Build jacobian from initial guess
    call build_jacobian_matrix(xvec)

    ! Set up Jacobian for Petsc
    call jac_prec % setup_petsc()

    ! write all out
    call loss % write_petsc_binary('loss.bin')
    call prod % write_petsc_binary('prod.bin')
    call jac_prec % write_petsc_binary('jac.bin')
!   call xvec % write_petsc_binary('x.bin')

  end subroutine init_data

!==============================================================================
! INIT_JACOBIAN_MATRIX preallocates jacobian matrix and initializes it
!==============================================================================

  subroutine init_jacobian_matrix()

    integer :: nnz
    integer :: n

    ! Get length of matrix and number of nonzeros total in loss matrix
    nnz = loss % nnz
    n   = loss % n

    ! There is one more nonzero for each row and and additional row of nonzeros
    nnz = nnz + 2*n + 1

    ! We need 1 more row for the jacobian
    n = n + 1

    ! Configure Jacobian matrix
    call jac_prec % create(n, nnz)

  end subroutine init_jacobian_matrix

!===============================================================================
! BUILD_JACOBIAN_MATRIX creates the Jacobian of nonlinear eigenvalue problem
!===============================================================================

  subroutine build_jacobian_matrix(x)

    use constants,  only: ONE

    type(Vector) :: x

    integer      :: i      ! loop counter for jacobian rows
    integer      :: jjac   ! loop counter for jacobian cols
    integer      :: jloss  ! loop counter for loss matrix cols
    integer      :: jprod  ! loop counter for prod matrix cols
    integer      :: n      ! problem size
    integer      :: shift  ! integer to account for index shift for Petsc
    real(8)      :: lambda ! eigenvalue
    real(8)      :: val    ! temporary real scalar
    type(Vector) :: flux   ! flux vector
    type(Vector) :: fsrc   ! fission source vector
print *, 'IN JACOBIAN'
    ! get the problem size
    n = loss % n

    ! Check to use shift
    shift = 0
    if (loss % petsc_active) shift = 1

    ! get flux and eigenvalue
    flux % val => xvec % val(1:n)
    lambda = xvec % val(n + 1)

    ! create fission source vector
    call fsrc % create(n)
    call prod % vector_multiply(flux, fsrc)

    ! reset counters in Jacobian matrix
    jac_prec % n_kount  = 1
    jac_prec % nz_kount = 1

    ! begin loop around rows of Jacobian matrix
    ROWS: do i = 1, n

      ! add a new row in Jacobian
      call jac_prec % new_row()

      ! begin loop around columns of loss matrix
      COLS_LOSS: do jloss = loss % row(i) + shift, &
                            loss % row(i + 1) - 1 + shift

        ! start with the value in the loss matrix
        val = loss % val(jloss)

        ! loop around columns of prod matrix
        COLS_PROD: do jprod = prod % row(i) + shift, &
                              prod % row(i + 1) - 1 + shift

          ! see if columns agree with loss matrix
          if (prod % col(jprod) == loss % col(jloss)) then
            val = val - lambda*prod % val(jprod)
            exit
          end if

        end do COLS_PROD 

        ! record value in jacobian
        call jac_prec % add_value(loss % col(jloss) + shift, val)

      end do COLS_LOSS

      ! add fission source value
      val = -fsrc % val(i)
      call jac_prec % add_value(n + 1, val)

    end do ROWS

    ! need to add negative transpose of flux vector on last row
    call jac_prec % new_row()
    do jjac = 1, n
      val = -flux % val(jjac)
      call jac_prec % add_value(jjac, val)
    end do

    ! add unity on bottom right corner of matrix
    call jac_prec % add_value(n + 1, ONE)

    ! CRS requires a final value in row
    call jac_prec % new_row()

    ! If the Jacobian is already associated with PETSc need to reshift row/col
    if (jac_prec % petsc_active) then
      jac_prec % row = jac_prec % row - 1
      jac_prec % col = jac_prec % col - 1
    end if

    ! free all allocated memory
    flux % val => null()
    call fsrc % destroy()

  end subroutine build_jacobian_matrix

!===============================================================================
! COMPUTE_NONLINEAR_RESIDUAL
!===============================================================================

  subroutine compute_nonlinear_residual(x, res)

    use global,       only: cmfd_write_matrices

    type(Vector) :: x
    type(Vector) :: res

    character(len=25) :: filename
    integer           :: n
    real(8)           :: lambda
    type(Vector)      :: res_loss
    type(Vector)      :: res_prod
    type(Vector)      :: flux

    type(Vector), pointer :: res_ptr
    type(Vector), pointer :: x_ptr

    res_ptr => resvec
    x_ptr => xvec

print *, 'IN RESIDUAL'

    ! Get problem size
    n = loss % n

    ! Set up temporary vectors
    call res_loss % create(n)
    call res_prod % create(n)

    ! Extract flux
    flux % n = n
    flux % val => x % val(1:n)

    ! Extract eigenvalue
    lambda = x % val(n + 1)

    ! Calculate M*flux then F*flux
    call loss % vector_multiply(flux, res_loss)
    call prod % vector_multiply(flux, res_prod)

    ! Put flux component values in residual vector
    res % val(1:n) = res_loss % val - lambda*res_prod % val

    ! Put eigenvalue component in residual vector
    res % val(n+1) = 0.5_8 - 0.5_8*sum(flux % val * flux % val)

    ! write out data in binary files (debugging)
    if (cmfd_write_matrices) then

      ! write out residual vector 
      if (adjoint_calc) then
        filename = 'adj_res.bin'
      else
        filename = 'res.bin'
      end if
      call res % write_petsc_binary(filename)

      ! write out solution vector
      if (adjoint_calc) then
        filename = 'adj_x.bin'
      else
        filename = 'x.bin'
      end if
      call x % write_petsc_binary(filename)

    end if
!   print *, x % val
!   print *, res % val
  end subroutine compute_nonlinear_residual

!===============================================================================
! COMPUTE_ADJOINT
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    ! transpose matrices
    call loss % transpose()
    call prod % transpose()

    ! write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call loss % write_petsc_binary('adj_lossmat.bin')
      call loss % write_petsc_binary('adj_prodmat.bin')
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXTRACT_RESULTS 
!===============================================================================

  subroutine extract_results()

    use constants, only: ZERO, ONE
    use global,    only: cmfd

    integer :: n ! problem size
call loss % write_petsc_binary('loss.mat')
call prod % write_petsc_binary('prod.mat')
call xvec % write_petsc_binary('x.bin')
call resvec % write_petsc_binary('res.bin')
    ! get problem size
    n = loss % n

    ! also allocate in cmfd object
    if (adjoint_calc) then
      if (.not. allocated(cmfd % adj_phi)) allocate(cmfd % adj_phi(n))
    else
      if (.not. allocated(cmfd % phi)) allocate(cmfd % phi(n))
    end if

    ! Get flux and eigenvalue
    if (adjoint_calc) then
      cmfd % adj_phi  = xvec % val(1:n)
      cmfd % adj_keff = ONE / xvec % val(n+1)
    else
      cmfd % phi  = xvec % val(1:n)
      cmfd % keff = ONE / xvec % val(n+1)
    end if
stop
  end subroutine extract_results

!===============================================================================
! FINALIZE frees all memory from a JFNK calculation
!===============================================================================

  subroutine finalize()

    call loss % destroy()
    call prod % destroy()
    call jac_prec % destroy()
    call xvec % destroy()
    call resvec % destroy()
    call jfnk % destroy() 
    nullify(jfnk_data % res_proc_ptr)
    nullify(jfnk_data % jac_proc_ptr)

  end subroutine finalize

end module cmfd_jfnk_solver
