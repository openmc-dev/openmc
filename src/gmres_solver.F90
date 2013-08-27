module gmres_solver

  implicit none

contains

subroutine pgmres ( n, im, rhs, sol, vv, eps, maxits, iout, &
  aa, ja, ia, alu, jlu, ju, ierr )

!*****************************************************************************80
!
!! PGMRES is an ILUT - Preconditioned GMRES solver.
!                                                                      
!  Discussion:
!
!    This is a simple version of the ILUT preconditioned GMRES algorithm. 
!    The ILUT preconditioner uses a dual strategy for dropping elements   
!    instead  of the usual level of-fill-in approach. See details in ILUT 
!    subroutine documentation. PGMRES uses the L and U matrices generated 
!    from the subroutine ILUT to precondition the GMRES algorithm.        
!    The preconditioning is applied to the right. The stopping criterion  
!    utilized is based simply on reducing the residual norm by epsilon.   
!    This preconditioning is more reliable than ilu0 but requires more    
!    storage. It seems to be much less prone to difficulties related to   
!    strong nonsymmetries in the matrix. We recommend using a nonzero tol 
!    (tol=.005 or .001 usually give good results) in ILUT. Use a large    
!    lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    
!    more reliable the code is. Efficiency may also be much improved.     
!    Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as 
!    Gaussian elimination without pivoting.                               
!                                                                      
!    ILU(0) and MILU(0) are also provided for comparison purposes         
!    USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and 
!    then call pgmres.                                                    
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad              
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IM, the size of the Krylov subspace.  IM 
!    should not exceed 50 in this version.  This restriction can be reset by 
!    changing the parameter command for KMAX below.
!                                                
!    Input/output, real RHS(N), on input, the right hand side vector.
!    On output, the information in this vector has been destroyed.
!
! sol   == real vector of length n containing an initial guess to the  
!          solution on input. approximate solution on output           
!
! eps   == tolerance for stopping criterion. process is stopped        
!          as soon as ( ||.|| is the euclidean norm):                  
!          || current residual||/||initial residual|| <= eps           
!
! maxits== maximum number of iterations allowed                        
!
! iout  == output unit number number for printing intermediate results 
!          if (iout <= 0) nothing is printed out.                    
!                                                                      
!    Input, real AA(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!                                                                      
! alu,jlu== A matrix stored in Modified Sparse Row format containing   
!           the L and U factors, as computed by routine ilut.       
!                                                                      
! ju     == integer ( kind = 4 ) array of length n containing the pointers to       
!           the beginning of each row of U in alu, jlu as computed     
!           by routine ILUT.                                        
!                                                                      
! on return:                                                           
!                                                          
! sol   == contains an approximate solution (upon successful return).  
! ierr  == integer ( kind = 4 ). Error message with the following meaning.          
!          ierr = 0 --> successful return.                            
!          ierr = 1 --> convergence not achieved in itmax iterations. 
!          ierr =-1 --> the initial guess seems to be the exact        
!                       solution (initial residual computed was zero) 
!                                                                      
! work arrays:                                                        
!                                                       
! vv    == work array of length  n x (im+1) (used to store the Arnoli  
!          basis)                                                      
!
  integer ( kind = 4 ), parameter :: kmax = 50
  integer ( kind = 4 ) n

  real ( kind = 8 ) aa(*)
  real ( kind = 8 ) alu(*)
  real ( kind = 8 ) c(kmax)
  real ( kind = 8 ) eps
  real ( kind = 8 ) eps1
  real ( kind = 8 ), parameter :: epsmac = 1.0D-16
  real ( kind = 8 ) gam
  real ( kind = 8 ) hh(kmax+1,kmax)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) im
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) its
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) maxits
  integer ( kind = 4 ) n1
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) ro
  real ( kind = 8 ) rs(kmax+1)
  real ( kind = 8 ) s(kmax)
  real ( kind = 8 ) sol(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) vv(n,*)
!
!  Arnoldi size should not exceed KMAX=50 in this version.
!  To reset modify parameter KMAX accordingly.
!
  n1 = n + 1
  its = 0
!
!  Outer loop starts here.
!  Compute initial residual vector.
!
  call ope ( n, sol, vv, aa, ja, ia )

  vv(1:n,1) = rhs(1:n) - vv(1:n,1)

  do

    ro = sqrt ( ddot ( n, vv, 1, vv, 1 ) )

    if ( 0 < iout .and. its == 0 ) then
      write(iout, 199) its, ro
    end if

    if ( ro == 0.0D+00 ) then
      ierr = -1
      exit
    end if

    t = 1.0D+00 / ro
    vv(1:n,1) = vv(1:n,1) * t

    if ( its == 0 ) then
      eps1 = eps * ro
    end if
!
!  Initialize first term of RHS of Hessenberg system.
!
     rs(1) = ro
     i = 0

 4   continue

     i = i + 1
     its = its + 1
     i1 = i + 1
     call lusol0 ( n, vv(1,i), rhs, alu, jlu, ju )
     call ope ( n, rhs, vv(1,i1), aa, ja, ia )
!
!  Modified Gram - Schmidt.
!
     do j = 1, i
       t = ddot ( n, vv(1,j), 1, vv(1,i1), 1 )
       hh(j,i) = t
       call daxpy ( n, -t, vv(1,j), 1, vv(1,i1), 1 )
     end do

     t = sqrt ( ddot ( n, vv(1,i1), 1, vv(1,i1), 1 ) )
     hh(i1,i) = t

     if ( t /= 0.0D+00 ) then
       t = 1.0D+00 / t
       vv(1:n,i1) = vv(1:n,i1) * t
     end if
!
!  Update factorization of HH.
!
    if ( i == 1 ) then
      go to 121
    end if
!
!  Perform previous transformations on I-th column of H.
!
    do k = 2, i
       k1 = k-1
       t = hh(k1,i)
       hh(k1,i) = c(k1) * t + s(k1) * hh(k,i)
       hh(k,i) = -s(k1) * t + c(k1) * hh(k,i)
    end do

121 continue

    gam = sqrt ( hh(i,i)**2 + hh(i1,i)**2 )
!
!  If GAMMA is zero then any small value will do.
!  It will affect only residual estimate.
!
    if ( gam == 0.0D+00 ) then
      gam = epsmac
    end if
!
!  Get the next plane rotation.
!
    c(i) = hh(i,i) / gam
    s(i) = hh(i1,i) / gam
    rs(i1) = -s(i) * rs(i)
    rs(i) = c(i) * rs(i)
!
!  Determine residual norm and test for convergence.
!
    hh(i,i) = c(i) * hh(i,i) + s(i) * hh(i1,i)
    ro = abs ( rs(i1) )
131 format(1h ,2e14.4)

    if ( 0 < iout ) then
      write(iout, 199) its, ro
    end if

    if ( i < im .and. eps1 < ro ) then
      go to 4
    end if
!
!  Now compute solution.  First solve upper triangular system.
!
    rs(i) = rs(i) / hh(i,i)

    do ii = 2, i
      k = i - ii + 1
      k1 = k + 1
      t = rs(k)
      do j = k1, i
        t = t - hh(k,j) * rs(j)
      end do
      rs(k) = t / hh(k,k)
    end do
!
!  Form linear combination of V(*,i)'s to get solution.
!
    t = rs(1)
    rhs(1:n) = vv(1:n,1) * t

    do j = 2, i
      t = rs(j)
      rhs(1:n) = rhs(1:n) + t * vv(1:n,j)
    end do
!
!  Call preconditioner.
!
    call lusol0 ( n, rhs, rhs, alu, jlu, ju )

    sol(1:n) = sol(1:n) + rhs(1:n)
!
!  Restart outer loop when necessary.
!
    if ( ro <= eps1 ) then
      ierr = 0
      exit
    end if

    if ( maxits < its ) then
      ierr = 1
      exit
    end if
!
!  Else compute residual vector and continue.
!
    do j = 1, i
      jj = i1 - j + 1
      rs(jj-1) = -s(jj-1) * rs(jj)
      rs(jj) = c(jj-1) * rs(jj)
    end do

    do j = 1, i1
      t = rs(j)
      if ( j == 1 ) then
        t = t - 1.0D+00
      end if
      call daxpy ( n, t, vv(1,j), 1,  vv, 1 )
    end do

199 format(' its =', i4, ' res. norm =', G14.6)

  end do

  return
end

subroutine ilut ( n, a, ja, ia, lfil, tol, alu, jlu, ju, iwk, wu, wl, jr, &
  jwl, jwu, ierr )

!*****************************************************************************80
!
!! ILUT is an ILUT preconditioner.
!
!  Discussion:
!
!    This routine carries ouot incomplete LU factorization with dual
!    truncation mechanism.  Sorting is done for both L and U. 
!
!    The dual drop-off strategy works as follows:
!
!    1) Theresholding in L and U as set by TOL.  Any element whose size
!       is less than some tolerance (relative to the norm of current
!       row in u) is dropped.
!         
!    2) Keeping only the largest lenl0+lfil elements in L and the
!       largest lenu0+lfil elements in U, where lenl0=initial number 
!       of nonzero elements in a given row of lower part of A 
!       and lenlu0 is similarly defined.
!  
!    Flexibility: one can use tol=0 to get a strategy based on keeping the
!    largest elements in each row of L and U. Taking tol /= 0 but lfil=n
!    will give the usual threshold strategy (however, fill-in is then     
!    unpredictible).                                                      
!
!    A must have all nonzero diagonal elements.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!                                                                      
!  Parameters:
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! lfil    = integer ( kind = 4 ). The fill-in parameter. Each row of L and
!           each row of U will have a maximum of lfil elements
!           in addition to the original number of nonzero elements.
!           Thus storage can be determined beforehand.
!           lfil must be >= 0.
!
! iwk     = integer ( kind = 4 ). The minimum length of arrays alu and jlu
!
! On return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer ( kind = 4 ) array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! ierr    = integer ( kind = 4 ). Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr > 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is >  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero pivot encountered.
!
! work arrays:
!
! jr,jwu,jwl, integer ( kind = 4 ) work arrays of length n.
! wu, wl, real work arrays of length n+1, and n resp.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) alu(*)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) idiag
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iwk
  integer ( kind = 4 ) j 
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) jpos
  integer ( kind = 4 ) jr(*)
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) ju0
  integer ( kind = 4 ) jwl(n)
  integer ( kind = 4 ) jwu(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) len
  integer ( kind = 4 ) lenl
  integer ( kind = 4 ) lenl0
  integer ( kind = 4 ) lenu
  integer ( kind = 4 ) lenu0
  integer ( kind = 4 ) lfil
  integer ( kind = 4 ) nl
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) tnorm
  real ( kind = 8 ) tol
  real ( kind = 8 ) wl(n)
  real ( kind = 8 ) wu(n+1)

  if ( lfil < 0 ) then
    ierr = -4
    return
  end if
!
!  Initialize JU0 (points to next element to be added to ALU, JLU)
!  and pointer.
!
  ju0 = n + 2
  jlu(1) = ju0
!
!  integer ( kind = 4 ) double pointer array.
!
  jr(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    j1 = ia(ii)
    j2 = ia(ii+1) - 1
    lenu = 0
    lenl = 0

    tnorm = 0.0D+00
    do k = j1, j2
      tnorm = tnorm + abs ( a(k) )
    end do
    tnorm = tnorm / real ( j2-j1+1, kind = 8 )
!
!  Unpack L-part and U-part of row of A in arrays WL, WU.
!
    do j = j1, j2

      k = ja(j)
      t = a(j)

      if ( tol * tnorm <= abs ( t ) ) then

        if ( k < ii ) then
          lenl = lenl + 1
          jwl(lenl) = k
          wl(lenl) = t
          jr(k) = lenl
        else
          lenu = lenu+1
          jwu(lenu) = k
          wu(lenu) = t
          jr(k) = lenu
        end if

      end if

    end do

    lenl0 = lenl
    lenu0 = lenu
    jj = 0
    nl = 0
!
!  Eliminate previous rows.
!
150 continue

    jj = jj + 1

    if ( lenl < jj ) then
      go to 160
    end if
!
!  In order to do the elimination in the correct order we need to
!  exchange the current row number with the one that has
!  smallest column number, among JJ, JJ+1, ..., LENL.
!
    jrow = jwl(jj)
    k = jj
!
!  Determine the smallest column index.
!
    do j = jj+1, lenl
       if ( jwl(j) < jrow ) then
          jrow = jwl(j)
          k = j
       end if
    end do
!
!  Exchange in JWL.
!
    j = jwl(jj)
    jwl(jj) = jrow
    jwl(k) = j
!
!  Exchange in JR.
!
    jr(jrow) = jj
    jr(j) = k
!
!  Exchange in WL.
!
    s = wl(k)
    wl(k) = wl(jj)
    wl(jj) = s

    if ( ii <= jrow ) then
      go to 160
    end if
!
!  Get the multiplier for row to be eliminated: JROW.
!
    fact = wl(jj) * alu(jrow)
    jr(jrow) = 0

    if ( abs ( fact ) * wu(n+2-jrow) <= tol * tnorm ) then
      go to 150
    end if
!
!  Combine current row and row JROW.
!
    do k = ju(jrow), jlu(jrow+1)-1
       s = fact * alu(k)
       j = jlu(k)
       jpos = jr(j)
!
!  If fill-in element and small disregard.
!
       if ( abs ( s ) < tol * tnorm .and. jpos == 0 ) then
         cycle
       end if

       if ( ii <= j ) then
!
!  Dealing with upper part.
!
          if ( jpos == 0 ) then
!
!  This is a fill-in element.
!
             lenu = lenu + 1

             if ( n < lenu ) then
               go to 995
             end if

             jwu(lenu) = j
             jr(j) = lenu
             wu(lenu) = - s
          else
!
!  No fill-in element.
!
             wu(jpos) = wu(jpos) - s
          end if
       else
!
!  Dealing with lower part.
!
          if ( jpos == 0 ) then
!
!  This is a fill-in element.
!
             lenl = lenl + 1

             if ( n < lenl ) then
               go to 995
             end if

             jwl(lenl) = j
             jr(j) = lenl
             wl(lenl) = -s
          else
!
!  No fill-in element.
!
             wl(jpos) = wl(jpos) - s
          end if
       end if

  end do

    nl = nl + 1
    wl(nl) = fact
    jwl(nl) = jrow
  go to 150
!
!  Update the L matrix.
!
 160 continue

    len = min ( nl, lenl0 + lfil )

    call bsort2 ( wl, jwl, nl, len )

    do k = 1, len

       if ( iwk < ju0 ) then
         ierr = -2
         return
       end if

       alu(ju0) =  wl(k)
       jlu(ju0) =  jwl(k)
       ju0 = ju0 + 1

    end do
!
!  Save pointer to beginning of row II of U.
!
    ju(ii) = ju0
!
!  Reset double pointer JR to zero (L-part - except first
!  JJ-1 elements which have already been reset).
!
  do k = jj, lenl
    jr(jwl(k)) = 0
  end do
!
!  Be sure that the diagonal element is first in W and JW.
!
    idiag = jr(ii)

    if ( idiag == 0 ) then
      go to 900
    end if

    if ( idiag /= 1 ) then

       s = wu(1)
       wu(j) = wu(idiag)
       wu(idiag) = s

       j = jwu(1)
       jwu(1) = jwu(idiag)
       jwu(idiag) = j

    end if

    len = min ( lenu, lenu0 + lfil )

    call bsort2 ( wu(2), jwu(2), lenu-1, len )
!
! Update the U-matrix.
!
    t = 0.0D+00

    do k = 2, len

       if ( iwk < ju0 ) then
         ierr = -3
         return
       end if

       jlu(ju0) = jwu(k)
       alu(ju0) = wu(k)
       t = t + abs ( wu(k) )
       ju0 = ju0 + 1

    end do
!
!  Save norm in WU (backwards). Norm is in fact average absolute value.
!
    wu(n+2-ii) = t / real ( len + 1, kind = 8 )
!
!  Store inverse of diagonal element of U.
!
    if ( wu(1) == 0.0D+00 ) then
      ierr = -5
      return
    end if

    alu(ii) = 1.0D+00 / wu(1)
!
!  Update pointer to beginning of next row of U.
!
  jlu(ii+1) = ju0
!
!  Reset double pointer JR to zero (U-part).
!
  do k = 1, lenu
    jr(jwu(k)) = 0
  end do

  end do

  ierr = 0

  return
!
!  Zero pivot :
!
 900    ierr = ii
    return
!
!  Incomprehensible error. Matrix must be wrong.
!
 995    ierr = -1
    return
end

subroutine amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do

    y(i) = t

  end do

  return
end
subroutine amuxd ( n, x, y, diag, ndiag, idiag, ioff )

!*****************************************************************************80
!
!! AMUXD multiplies a DIA matrix times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector when the original matrix 
!    is stored in the DIA diagonal storage format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), array of length equal to the column dimension of
!    the A matrix.
!
!    Output, real Y(N), the product A * X.
!
!    Input, real DIAG(NDIAG,IDIAG), the diagonals.
!
!    Input, integer ( kind = 4 ) NDIAG, the first dimension of array adiag as 
!    declared in the calling program.
!
!    Input, integer ( kind = 4 ) IDIAG, the number of diagonals in the matrix.
!
!    Input, integer ( kind = 4 ) IOFF(IDIAG), the offsets of the diagonals of 
!    the matrix: diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
  integer ( kind = 4 ) idiag
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndiag

  real ( kind = 8 ) diag(ndiag,idiag)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) io
  integer ( kind = 4 ) ioff(idiag)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do j = 1, idiag
    io = ioff(j)
    i1 = max ( 1, 1 - io )
    i2 = min ( n, n - io )
    do k = i1, i2
      y(k) = y(k) + diag(k,j) * x(k+io)
    end do
  end do

  return
end

subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    Uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Dongarra, Moler, Bunch, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979.
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DA, the multiplier of DX.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of DX.
!
!    Input/output, real ( kind = 8 ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive entries of DY.
!
  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    return
  end if

  if ( da  == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    do i = 1, m
      dy(i) = dy(i) + da * dx(i)
    end do

    do i = m+1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end

subroutine bsort2 ( w, ind, n, ncut )

!*****************************************************************************80
!
!! BSORT2 returns the NCUT largest elements of an array, using bubble sort.
!
!  Discussion:
!
!    This routine carries out a simple bubble sort for getting the NCUT largest
!    elements in modulus, in array W.  IND is sorted accordingly.
!    (Ought to be replaced by a more efficient sort especially
!    if NCUT is not that small).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(*)
  integer ( kind = 4 ) iswp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncut
  logical test
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) wswp

  i = 1

  do

    test = .false.

    do j = n-1, i, -1

      if ( abs ( w(j) ) < abs ( w(j+1) ) ) then
!
!  Swap.
!
        wswp = w(j)
        w(j) = w(j+1)
        w(j+1) = wswp
!
!  Reorder the original ind array accordingly.
!
        iswp = ind(j)
        ind(j) = ind(j+1)
        ind(j+1) = iswp
!
!  Set indicator that sequence is still unsorted.
!
        test = .true.

      end if

    end do

    i = i + 1

    if ( .not. test .or. ncut < i ) then
      exit
    end if

  end do

  return
end

subroutine ope ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! OPE sparse matrix * vector multiplication
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real Y(N), the product A * X.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) -1
    y(i) = 0.0D+00
    do k = k1, k2
      y(i) = y(i) + a(k) * x(ja(k))
    end do
  end do

  return
end

subroutine lusol0 ( n, y, x, alu, jlu, ju )

!*****************************************************************************80
!
!! LUSOL0 performs a forward followed by a backward solve
! for LU matrix as produced by  ILUT
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution.
!
!    ALU, JLU, JU, ...
!
  integer ( kind = 4 ) n

  real ( kind = 8 ) alu(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Forward solve
!
  do i = 1, n
    x(i) = y(i)
    do k = jlu(i), ju(i)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
  end do
!
!  Backward solve.
!
  do i = n, 1, -1
    do k = ju(i), jlu(i+1)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
    x(i) = alu(i) * x(i)
  end do

  return
end

function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Dongarra, Moler, Bunch, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries in X.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive entries in Y.
!
!    Output, real DDOT, the sum of the product of the corresponding
!    entries of X and Y.
!
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m+1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end

end module gmres_solver
