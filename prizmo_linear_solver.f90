module prizmo_linear_solver
contains

  ! *******************************
  ! solve Ax=B linear system
  function linear_solver(a, b) result(c)
    implicit none
    real*8,intent(in)::a(:, :), b(:)
    integer::n
    real*8::c(size(b))

    n = size(b)

    if(n==0) then
      print *, "ERROR: linear solver called with zero-order matix!"
      stop
    elseif(n==2) then
       c(:) = linear_solver_n2(a(:, :), b(:))
    elseif(n==3) then
       c(:) = linear_solver_n3(a(:, :), b(:))
    elseif(n==4) then
        c(:) = linear_solver_n4(a(:, :), b(:))
    elseif(n==5) then
        c(:) = linear_solver_n5(a(:, :), b(:))
    else
       c(:) = linear_solver_nn(a(:, :), b(:), n)
    end if

  end function linear_solver

  ! ***********************
  ! solve Ax=B analytically for a 2-levels system with first row of ones
  ! (  1,   1) * (c1) = (b1)
  ! (a21, a22)   (c2)   (0 )
  function linear_solver_n2(a, b) result(c)
    implicit none
    integer,parameter::n=2
    real*8,intent(in)::a(n, n), b(n)
    real*8::c(n), iab

    c(1) = a(2, 2) * b(1) / (a(2, 2) - a(2, 1))
    c(2) = b(1) - c(1)

  end function linear_solver_n2

  ! ************************
  ! solve Ax=B analytically for a 3-levels system with first row of ones
  ! (  1,   1,   1) * (c1) = (b1)
  ! (a21, a22, a23)   (c2)   (0 )
  ! (a31, a32, a33)   (c3)   (0 )
  function linear_solver_n3(a, b) result(c)
    implicit none
    integer,parameter::n=3
    real*8,intent(in)::a(n, n), b(n)
    real*8::iab, c(n)

    iab = b(1) / (a(2, 1) * (a(3, 3) - a(3, 2)) &
      + a(2, 2) * (a(3, 1) - a(3, 3)) &
      + a(2, 3) * (a(3, 2) - a(3, 1)))

    c(1) = (a(2, 3) * a(3, 2) - a(2, 2) * a(3, 3)) * iab
    c(2) = (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) * iab
    c(3) = b(1) - c(1) - c(2)

  end function linear_solver_n3

  ! ********************************
  function linear_solver_n4(a, b) result(c)
  implicit none
    integer,parameter::n=4
    real*8,intent(in)::a(n, n), b(n)
    real*8::det, c(n)
    real*8::a0, a1, a2, a3
    real*8::b0, b1, b2, b3
    real*8::c0, c1, c2, c3
    real*8::d01, d02, d03, d12, d13, d23

    a0 = a(2,1)
    a1 = a(2,2)
    a2 = a(2,3)
    a3 = a(2,4)
    b0 = a(3,1)
    b1 = a(3,2)
    b2 = a(3,3)
    b3 = a(3,4)
    c0 = a(4,1)
    c1 = a(4,2)
    c2 = a(4,3)
    c3 = a(4,4)

    d01 = c0 - c1
    d02 = c0 - c2
    d03 = c0 - c3
    d12 = c1 - c2
    d13 = c1 - c3
    d23 = c2 - c3

    det = -a0 * (b3 * d12 - b2 * d13 + b1 * d23) &
         + a1 * (b3 * d02 - b2 * d03 + b0 * d23) &
         - a2 * (b3 * d01 - b1 * d03 + b0 * d13) &
         + a3 * (b2 * d01 - b1 * d02 + b0 * d12)

    c(2) = a3 * (b2*c0 - b0*c2) - a2 * (b3*c0 - b0*c3) + a0 * (b3*c2 - b2*c3)
    c(3) = a3 * (b0*c1 - b1*c0) - a1 * (b0*c3 - b3*c0) + a0 * (b1*c3 - b3*c1)
    c(4) = a2 * (b1*c0 - b0*c1) - a1 * (b2*c0 - b0*c2) + a0 * (b2*c1 - b1*c2)
    c(1) = det - sum(c(2:n))
    c = b(1) * c / det

  end function linear_solver_n4

  ! ********************************
  function linear_solver_n5(a, b) result(c)
    implicit none
    integer,parameter::n=5
    real*8,intent(in)::a(n, n), b(n)
    real*8::det, c(n)
    real*8::a0, a1, a2, a3, a4
    real*8::b0, b1, b2, b3, b4
    real*8::c0, c1, c2, c3, c4
    real*8::d0, d1, d2, d3, d4
    real*8::h01, h02, h03, h04, h12, h13, h14, h23, h24, h34
    real*8::p012, p013, p014, p023, p024, p034, p102, p103, p104
    real*8::p123, p124, p134, p203, p204, p210, p213, p214, p304
    real*8::p310, p314, p320, p321, p410, p420, p421, p430, p431
    real*8::p432, p324, p234
    real*8::q210, q310, q320, q321, q410, q420, q421, q430, q431, q432
    real*8::w01, w02, w03, w04, w12, w13, w14, w23, w24, w34

    a0 = a(2,1)
    a1 = a(2,2)
    a2 = a(2,3)
    a3 = a(2,4)
    a4 = a(2,5)
    b0 = a(3,1)
    b1 = a(3,2)
    b2 = a(3,3)
    b3 = a(3,4)
    b4 = a(3,5)
    c0 = a(4,1)
    c1 = a(4,2)
    c2 = a(4,3)
    c3 = a(4,4)
    c4 = a(4,5)
    d0 = a(5,1)
    d1 = a(5,2)
    d2 = a(5,3)
    d3 = a(5,4)
    d4 = a(5,5)

    h01 = d0 - d1
    h02 = d0 - d2
    h03 = d0 - d3
    h04 = d0 - d4
    h12 = d1 - d2
    h13 = d1 - d3
    h14 = d1 - d4
    h23 = d2 - d3
    h24 = d2 - d4
    h34 = d3 - d4

    p012 = + c0 * h12
    p013 = + c0 * h13
    p014 = + c0 * h14
    p023 = + c0 * h23
    p024 = + c0 * h24
    p034 = + c0 * h34
    p102 = + c1 * h02
    p103 = + c1 * h03
    p104 = + c1 * h04
    p123 = + c1 * h23
    p124 = + c1 * h24
    p134 = + c1 * h34
    p203 = + c2 * h03
    p204 = + c2 * h04
    p210 = - c2 * h01
    p213 = + c2 * h13
    p214 = + c2 * h14
    p234 = + c2 * h34
    p304 = + c3 * h04
    p310 = - c3 * h01
    p314 = + c3 * h14
    p320 = - c3 * h02
    p321 = - c3 * h12
    p324 = + c3 * h24
    p410 = - c4 * h01
    p420 = - c4 * h02
    p421 = - c4 * h12
    p430 = - c4 * h03
    p431 = - c4 * h13
    p432 = - c4 * h23

    q210 = p210 + p102 - p012
    q310 = p310 + p103 - p013
    q320 = p320 + p203 - p023
    q321 = p321 + p213 - p123
    q410 = p410 + p104 - p014
    q420 = p420 + p204 - p024
    q421 = p421 + p214 - p124
    q430 = p430 + p304 - p034
    q431 = p431 + p314 - p134
    q432 = p432 + p324 - p234

    w01 = b0*d1 - b1*d0
    w02 = b0*d2 - b2*d0
    w03 = b0*d3 - b3*d0
    w04 = b0*d4 - b4*d0
    w12 = b1*d2 - b2*d1
    w13 = b1*d3 - b3*d1
    w14 = b1*d4 - b4*d1
    w23 = b2*d3 - b3*d2
    w24 = b2*d4 - b4*d2
    w34 = b3*d4 - b4*d3

    det =  a4 * (b3 * q210 - b2 * q310 + b1 * q320 - b0 * q321) &
         - a3 * (b4 * q210 - b2 * q410 + b1 * q420 - b0 * q421) &
         + a2 * (b4 * q310 - b3 * q410 + b1 * q430 - b0 * q431) &
         - a1 * (b4 * q320 - b3 * q420 + b2 * q430 - b0 * q432) &
         + a0 * (b4 * q321 - b3 * q421 + b2 * q431 - b1 * q432)

    c(2) = c4 * (a3*w02 - a2*w03 + a0*w23) &
         - c3 * (a4*w02 - a2*w04 + a0*w24) &
         + c2 * (a4*w03 - a3*w04 + a0*w34) &
         - c0 * (a4*w23 - a3*w24 + a2*w34)

    c(3) = -c4 * (a3*w01 - a1*w03 + a0*w13) &
         +  c3 * (a4*w01 - a1*w04 + a0*w14) &
         -  c1 * (a4*w03 - a3*w04 + a0*w34) &
         +  c0 * (a4*w13 - a3*w14 + a1*w34)

    c(4) = c4 * (a2*w01 - a1*w02 + a0*w12) &
         - c2 * (a4*w01 - a1*w04 + a0*w14) &
         + c1 * (a4*w02 - a2*w04 + a0*w24) &
         - c0 * (a4*w12 - a2*w14 + a1*w24)

    c(5) = -c3 * (a2*w01 - a1*w02 + a0*w12) &
         +  c2 * (a3*w01 - a1*w03 + a0*w13) &
         -  c1 * (a3*w02 - a2*w03 + a0*w23) &
         +  c0 * (a3*w12 - a2*w13 + a1*w23)

    c(1) = det - sum(c(2:n))

    c = b(1) * c / det

  end function linear_solver_n5


  ! ************************
  ! solve Ax=B for a n-levels system
  function linear_solver_nn(a, b, n) result(c)
    implicit none
    integer,intent(in)::n
    real*8,intent(in)::a(n, n), b(n)
    real*8::c(n)
    integer::ipiv(n), ierr, i, j

    c(:) = b(:)
    print *, "LAPACK DGESV not implemented!"
    stop

    ! call dgesv(n, 1, a(:, :), n, ipiv(:), c(:), n, ierr)
    !
    ! if(ierr /= 0) then
    !    print *, "ERROR: dgesv error", ierr
    !    if(ierr < 0) then
    !      print *, "the ith argument had an illegal value, ith=", ierr
    !    else
    !      print *, "U(i, i) is exactly zero (singular matrix), i=", ierr
    !      do i=1,n
    !        do j=1,n
    !          if(a(i,j) /= 0d0) then
    !            print '(a10,2I5,E17.8e3)', "A(i, j):", i, j, a(i, j)
    !          end if
    !        end do
    !      end do
    !    end if
    !    print *, "NOTE: the next floating point exception is triggered to get backtrace"
    !    c = -1d0
    !    c = log10(c)
    !    stop
    ! end if

  end function linear_solver_nn

  ! *************************
  ! solve a linear system as
  ! (a11, -a12,    0,    0)   (c1)   (0)
  ! (  0,  a22, -a23,    0)   (c2)   (0)
  ! (  0,    0,  a33, -a34)   (c3)   (0)
  ! (  1,    1,    1,    1) * (c4) = (b)
  function bidiagonal_solver(a, b, n) result(c)
    implicit none
    integer,intent(in)::n
    real*8,intent(in)::a(n, n), b
    real*8::c(n), w(n)
    integer::i

    ! compute coefficients
    w(2) = a(1, 1) / a(1, 2)
    do i=3,n
      w(i) = w(i-1) * a(i-1, i-1) / a(i-1, i)
    end do

    ! find unknowns
    c(1) = b / (1d0 + sum(w(2:n)))
    c(2:n) = w(2:n) * c(1)

  end function bidiagonal_solver

end module prizmo_linear_solver
