module prizmo_spline

  type, public :: spline1d_data_vec(nv, n1)
    integer,len::nv, n1
    real*8::fdata(nv, n1), xmin, invdx, dx, xmax
    real*8::b(nv, n1), c(nv, n1), d(nv, n1)
  end type

contains

  ! *************************
  subroutine load_spline1d_vec(fnames, nvec, nx, data)
    implicit none
    integer,intent(in)::nvec, nx
    character(len=*),intent(in)::fnames(nvec)
    real*8::xdata(nx), fdata(nvec, nx)
    real*8::b(nvec, nx), c(nvec, nx), d(nvec, nx)
    type(spline1d_data_vec(nv=nvec, n1=nx)),intent(inout)::data
    integer::i

    do i=1,nvec
      call load_1d(fnames(i), nx, xdata, fdata(i, :))
      call compute_coefficients(xdata, fdata(i, :), nx, b(i, :), c(i, :), d(i, :))
    end do

    data%xmin = minval(xdata)
    data%xmax = maxval(xdata)
    data%dx = xdata(2) - xdata(1)
    data%fdata = fdata
    data%invdx = 1d0 / data%dx
    data%b = b
    data%c = c
    data%d = d

  end subroutine load_spline1d_vec

  ! *************************
  subroutine load_1d(fname, nx, xdata, fdata)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nx
    real*8,intent(out)::xdata(nx), fdata(nx)
    integer::unit, i

    open(newunit=unit, file=trim(fname), status="old")
    do i=1,nx
      read(unit, *) xdata(i), fdata(i)
    end do
    close(unit)

  end subroutine load_1d

  ! ************************
  function interp_spline1d_vec(x, nvec, nx, data) result(f)
    implicit none
    integer,intent(in)::nvec, nx
    type(spline1d_data_vec(nv=nvec, n1=nx)),intent(inout)::data
    real*8,intent(in)::x
    real*8::f(nvec), dx
    integer::idx

    idx = (x - data%xmin) * data%invdx + 1
    dx = x - (data%xmin + (idx - 1) * data%dx)
    f = data%fdata(:, idx) + dx * (data%b(:, idx) + dx * (data%c(:, idx) + dx * data%d(:, idx)))

  end function interp_spline1d_vec

  ! *************************
  subroutine compute_coefficients(x, y, nx, b, c, d)
    implicit none
    integer,intent(in)::nx
    real*8,intent(in)::x(nx), y(nx)
    real*8,intent(out)::b(nx), c(nx), d(nx)
    real*8::h
    integer::i, j, gap

    gap = nx - 1

    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1)) / d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2.0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i)) / d(i)
      c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions
    b(1) = -d(1)
    b(nx) = -d(nx-1)
    c(1) = 0.0
    c(nx) = 0.0
    if(nx /= 3) then
      c(1) = c(3) / (x(4) - x(2)) - c(2) / (x(3) - x(1))
      c(nx) = c(nx-1) / (x(nx) - x(nx-2)) - c(nx-2) / (x(nx-1) - x(nx-3))
      c(1) = c(1) * d(1)**2 / (x(4) - x(1))
      c(nx) = -c(nx) * d(nx-1)**2 / (x(nx) - x(nx-3))
    end if

    ! step 3: forward elimination
    do i = 2, nx
      h = d(i-1) / b(i-1)
      b(i) = b(i) - h * d(i-1)
      c(i) = c(i) - h * c(i-1)
    end do

    ! step 4: back substitution
    c(nx) = c(nx) / b(nx)
    do j = 1, gap
      i = nx - j
      c(i) = (c(i) - d(i) * c(i+1)) / b(i)
    end do

    ! step 5: compute spline coefficients
    b(nx) = (y(nx) - y(gap)) / d(gap) + d(gap) * (c(gap) + 2.0*c(nx))
    do i = 1, gap
      b(i) = (y(i+1) - y(i)) / d(i) - d(i) * (c(i+1) + 2.0*c(i))
      d(i) = (c(i+1) - c(i)) / d(i)
      c(i) = 3. * c(i)
    end do
    c(nx) = 3.0*c(nx)
    d(nx) = d(nx-1)
  end subroutine compute_coefficients

  ! *************************
  subroutine spline_test()
    implicit none
    type(spline1d_data_vec(nv=2, n1=30))::spline_data
    real*8::x, f(2)
    integer::i

    call load_spline1d_vec((/"spline_cos.dat", "spline_sin.dat"/), &
                           2, 30, spline_data)

    do i=1,20
      x = (i - 1) / 19.
      f = interp_spline1d_vec(x, 2, 30, spline_data)
      print *, x, f
    end do

  end subroutine spline_test

end module prizmo_spline
