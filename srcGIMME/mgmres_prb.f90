subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real  A(NZ_NUM), the matrix values.
!
!    Input, real  X(N), the vector to be multiplied by A'.
!
!    Output, real  W(N), the value of A'*X.
!
  implicit none

  integer  n
  integer  nz_num

  real  a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  ja(nz_num)
  integer  k
  integer  k1
  integer  k2
  real  w(n)
  real  x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
  end do

  return
end
subroutine atx_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_ST computes A'*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real  A(NZ_NUM), the matrix values.
!
!    Input, real  X(N), the vector to be multiplied by A'.
!
!    Output, real  W(N), the value of A'*X.
!
  implicit none

  integer  n
  integer  nz_num

  real  a(nz_num)
  integer  i
  integer  ia(nz_num)
  integer  j
  integer  ja(nz_num)
  integer  k
  real  w(n)
  real  x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(j) = w(j) + a(k) * x(i)
  end do

  return
end
subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real  A(NZ_NUM), the matrix values.
!
!    Input, real  X(N), the vector to be multiplied by A.
!
!    Output, real  W(N), the value of A*X.
!
  implicit none

  integer  n
  integer  nz_num

  real  a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  ja(nz_num)
  integer  k
  integer  k1
  integer  k2
  real  w(n)
  real  x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end
subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_ST computes A*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real  A(NZ_NUM), the matrix values.
!
!    Input, real  X(N), the vector to be multiplied by A.
!
!    Output, real  W(N), the value of A*X.
!
  implicit none

  integer  n
  integer  nz_num

  real  a(nz_num)
  integer  i
  integer  ia(nz_num)
  integer  j
  integer  ja(nz_num)
  integer  k
  real  w(n)
  real  x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(i) = w(i) + a(k) * x(j)
  end do

  return
end
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

!*****************************************************************************80
!
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer  UA(N), the index of the diagonal element
!    of each row.
!
  implicit none

  integer  n
  integer  nz_num

  integer  i
  integer  ia(n+1)
  integer  k
  integer  ja(nz_num)
  integer  ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end
subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

!*****************************************************************************80
!
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real  A(NZ_NUM), the matrix values.
!
!    Input, integer  UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real  L(NZ_NUM), the ILU factorization of A.
!
  implicit none

  integer  n
  integer  nz_num

  real  a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  iw(n)
  integer  j
  integer  ja(nz_num)
  integer  jj
  integer  jrow
  integer  jw
  integer  k
  real  l(nz_num)
  real  tl
  integer  ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end
subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

!*****************************************************************************80
!
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Discussion:
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real  L(NZ_NUM), the matrix values.
!
!    Input, integer  UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real  R(N), the right hand side.
!
!    Output, real  Z(N), the solution of the system M * Z = R.
!
  implicit none

  integer  n
  integer  nz_num

  integer  i
  integer  ia(n+1)
  integer  j
  integer  ja(nz_num)
  real  l(nz_num)
  real  r(n)
  integer  ua(n)
  real  w(n)
  real  z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end
subroutine mgmres_st ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, tol_abs, &
  tol_rel )

!*****************************************************************************80
!
!! MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
!
!  Discussion:
!
!    The linear system A*X=B is solved iteratively.
!
!    The matrix A is assumed to be stored in sparse triplet form.  Only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real  A(NZ_NUM), the matrix values.
!
!    Input/output, real  X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real  RHS(N), the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations
!    to take.  0 < MR <= N.
!
!    Input, real  TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real  TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n
  integer  nz_num

  real  a(nz_num)
  real  av
  real  c(1:mr)
  real , parameter :: delta = 1.0D-03
  real  g(1:mr+1)
  real  h(1:mr+1,1:mr)
  real  htmp
  integer  i
  integer  ia(nz_num)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  real  mu
  real  r(1:n)
  real  rho
  real  rho_tol
  real  rhs(1:n)
  real  s(1:mr)
  real  tol_abs
  real  tol_rel
  real  v(1:n,1:mr+1)
  !logical, parameter :: verbose = .true.
  logical, parameter :: verbose = .false.
  real  x(1:n)
  real  y(1:mr+1)

  itr_used = 0

  if ( n < mr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
    write ( *, '(a)' ) '  N < MR.'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  MR = ', mr
    stop
  end if

  do itr = 1, itr_max

    call ax_st ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )

    if ( verbose ) then
      write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_st ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - h(j,k) * v(1:n,j)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( av + delta * h(k+1,k) == av ) then

        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do

        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then

        y(1:k+1) = h(1:k+1,k)

        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y(1:k+1) )
        end do

        h(1:k+1,k) = y(1:k+1)

      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g(1:k+1) )
      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)'       ) ' '
    write ( *, '(a)'       ) 'MGMRES_ST:'
    write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end
subroutine mult_givens ( c, s, k, g )

!*****************************************************************************80
!
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, real  C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer  K, indicates the location of the first
!    vector entry.
!
!    Input/output, real  G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer  k

  real  c
  real  g(1:k+1)
  real  g1
  real  g2
  real  s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end
subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel )

!*****************************************************************************80
!
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the linear system.
!
!    Input, integer  NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer  IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real  A(NZ_NUM), the matrix values.
!
!    Input/output, real  X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real  RHS(N), the right hand side of the linear system.
!
!    Input, integer  ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer  MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real  TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real  TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer  mr
  integer  n
  integer  nz_num

  real  a(nz_num)
  real  av
  real  c(mr+1)
  real , parameter :: delta = 1.0D-03
  real  g(mr+1)
  real  h(mr+1,mr)
  real  htmp
  integer  i
  integer  ia(n+1)
  integer  itr
  integer  itr_max
  integer  itr_used
  integer  j
  integer  ja(nz_num)
  integer  k
  integer  k_copy
  real  l(ia(n+1)+1)
  real  mu
  real  r(n)
  real  rho
  real  rho_tol
  real  rhs(n)
  real  s(mr+1)
  real  tol_abs
  real  tol_rel
  integer  ua(n)
  real  v(n,mr+1);
  !logical, parameter :: verbose = .true.
  logical, parameter :: verbose = .false.
  real  x(n)
  real  y(mr+1)

  itr_used = 0

  call rearrange_cr ( n, nz_num, ia, ja, a )

  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end subroutine pmgmres_ilu_cr
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real  values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer  N, the number of entries in the vector.
!
!    Input/output, integer  SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real  R(N), the vector of pseudorandom values.
!
  implicit none

  integer  n

  integer  i
  integer  k
  integer  seed
  real  r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed) * 4.656612875D-10

  end do

  return
end
subroutine rearrange_cr ( n, nz_num, ia, ja, a )

!*****************************************************************************80
!
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Discussion:
!
!    This routine guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer  N, the order of the system.
!
!    Input, integer  NZ_NUM, the number of nonzeros.
!
!    Input, integer  IA(N+1), the compressed row indices.
!
!    Input/output, integer  JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real  A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer  n
  integer  nz_num

  real  a(nz_num)
  integer  i
  integer  ia(n+1)
  integer  i4temp
  integer  ja(nz_num)
  integer  k
  integer  l
  real  r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer  d
  integer  h
  integer  m
  integer  mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer  n
  integer  s
  integer  values(8)
  integer  y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
