module prizmo_fit
  type, public :: fit1d_data(n1)
    integer,len::n1
    real*8::fdata(n1), xmin, invdx, dx, xmax, xfact
  end type

  type, public :: fit2d_data(n1, n2)
    integer,len::n1, n2
    real*8::fdata(n1, n2), xmin, invdx, dx, xmax, xfact
    real*8::ymin, invdy, dy, ymax, yfact
  end type

  type, public :: fit3d_data(n1, n2, n3)
    integer,len::n1, n2, n3
    real*8::fdata(n1, n2, n3), xmin, invdx, dx, xmax, xfact
    real*8::ymin, invdy, dy, ymax, yfact
    real*8::zmin, invdz, dz, zmax, zfact
  end type

  type, public :: fit4d_data(n1, n2, n3, n4)
    integer,len::n1, n2, n3, n4
    real*8::fdata(n1, n2, n3, n4), xmin, invdx, dx, xmax, xfact
    real*8::ymin, invdy, dy, ymax, yfact
    real*8::zmin, invdz, dz, zmax, zfact
    real*8::umin, invdu, du, umax, ufact
  end type

  type, public :: fit1d_data_vec(nv, n1)
    integer,len::nv, n1
    real*8::fdata(nv, n1), xmin, invdx, dx, xmax, xfact
  end type

  type, public :: fit2d_data_vec(nv, n1, n2)
    integer,len::nv, n1, n2
    real*8::fdata(nv, n1, n2), xmin, invdx, dx, xmax, xfact
    real*8::ymin, invdy, dy, ymax, yfact
  end type

  type, public :: fit3d_data_vec(nv, n1, n2, n3)
    integer,len::nv, n1, n2, n3
    real*8::fdata(nv, n1, n2, n3), xmin, invdx, dx, xmax, xfact
    real*8::ymin, invdy, dy, ymax, yfact
    real*8::zmin, invdz, dz, zmax, zfact
  end type

  type, public :: fit4d_data_vec(nv, n1, n2, n3, n4)
    integer,len::nv, n1, n2, n3, n4
    real*8::fdata(nv, n1, n2, n3, n4), xmin, invdx, dx, xmax, xfact
    real*8::ymin, invdy, dy, ymax, yfact
    real*8::zmin, invdz, dz, zmax, zfact
    real*8::umin, invdu, du, umax, ufact
  end type

contains

  ! ***************************
  subroutine load_1d_fit_vec(fname, nvec, nx, data, do_log, skip_lines)
    implicit none
    integer,intent(in)::nvec, nx
    character(len=*),intent(in)::fname(nvec)
    type(fit1d_data_vec(nv=nvec, n1=nx)),intent(inout)::data
    integer::i
    real*8::xdata(nx), fdata(nvec, nx), xmin, invdx, xfact
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    do i=1,nvec
      call load_data_1d(fname(i), fdata(i, :), xdata, nx, xmin, xfact, invdx, dolog, skipl)
    end do

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

  end subroutine load_1d_fit_vec

  ! ***************************
  subroutine load_1d_fit(fname, nx, data, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nx
    type(fit1d_data(n1=nx)),intent(inout)::data
    real*8::xdata(nx), fdata(nx), xmin, invdx, xfact
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    call load_data_1d(fname, fdata, xdata, nx, xmin, xfact, invdx, dolog, skipl)

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

  end subroutine load_1d_fit

  ! ********************
  subroutine load_data_1d(fname, fdata, xdata, nx, xmin, xfact, invdx, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nx
    integer::unit, i
    real*8::dx
    real*8,intent(out)::fdata(nx), xmin, xfact, invdx
    real*8,intent(out)::xdata(nx)
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    open(newunit=unit, file=trim(fname), status="old")
    do i=1,skipl
      read(unit, *)  ! read lines to be skipped
    end do
    do i=1,nx
      read(unit, *) xdata(i), fdata(i)
    end do
    close(unit)

    if(dolog) then
      xdata(:) = log10(xdata(:))
      fdata(:) = log10(fdata(:))
    end if

    call check_regular_grid(xdata(:), fname)

    xmin = xdata(1)
    dx = xdata(2) - xdata(1)
    xfact = (nx - 1) / (xdata(nx) - xdata(1))
    invdx = 1d0 / dx

  end subroutine load_data_1d

  ! *******************
  ! given xdata return variables to be used as arguments for the 1d fit routine
  subroutine fit_prepare(xdata, nx, xmin, xfact, invdx)
    implicit none
    integer,intent(in)::nx
    real*8,intent(in)::xdata(nx)
    real*8::xmin, dx, xfact, invdx

    ! check if grid is uniform
    call check_regular_grid(xdata(:), "unknown")

    xmin = xdata(1)
    dx = xdata(2) - xdata(1)
    xfact = (nx - 1) / (xdata(nx) - xdata(1))
    invdx = 1d0 / dx

  end subroutine fit_prepare

  ! ***************************
  subroutine load_2d_fit_vec(fname, nvec, nx, ny, data, mode, do_log, skip_lines)
    implicit none
    integer,intent(in)::nvec, nx, ny
    character(len=*),intent(in)::fname(nvec)
    character(len=2),intent(in),optional::mode
    character(len=2)::mode_read
    type(fit2d_data_vec(nv=nvec, n1=nx, n2=ny)),intent(inout)::data
    integer::i
    real*8::xdata(nx), fdata(nvec, nx, ny), xmin, invdx, xfact
    real*8::ydata(ny), ymin, invdy, yfact
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    mode_read = "xy"
    if(present(mode)) then
      mode_read = mode
    end if

    do i=1,nvec
      call load_data_2d(fname(i), fdata(i, :, :), &
        xdata, nx, xmin, xfact, invdx, &
        ydata, ny, ymin, yfact, invdy, mode_read, dolog, skipl)
    end do

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

    data%ymin = ymin
    data%ymax = ydata(ny)
    data%dy = ydata(2) - ydata(1)
    data%invdy = invdy
    data%yfact = yfact

  end subroutine load_2d_fit_vec

  ! ***************************
  subroutine load_2d_fit(fname, nx, ny, data, mode, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    character(len=2),intent(in),optional::mode
    character(len=2)::mode_read
    integer,intent(in)::nx, ny
    type(fit2d_data(n1=nx, n2=ny)),intent(inout)::data
    real*8::xdata(nx), fdata(nx, ny), xmin, invdx, xfact
    real*8::ydata(ny), ymin, invdy, yfact
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    mode_read = "xy"
    if(present(mode)) then
      mode_read = mode
    end if

    call load_data_2d(fname, fdata, &
      xdata, nx, xmin, xfact, invdx, &
      ydata, ny, ymin, yfact, invdy, mode_read, dolog, skipl)

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

    data%ymin = ymin
    data%ymax = ydata(ny)
    data%dy = ydata(2) - ydata(1)
    data%invdy = invdy
    data%yfact = yfact

  end subroutine load_2d_fit

  ! ********************
  subroutine load_data_2d(fname, fdata, &
    xdata, nx, xmin, xfact, invdx, &
    ydata, ny, ymin, yfact, invdy, mode, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    character(len=2),intent(in),optional::mode
    character(len=2)::mode_read
    integer,intent(in)::nx, ny
    integer::unit, i, j
    real*8::dx, dy
    real*8,intent(out)::fdata(nx, ny), xmin, xfact, invdx
    real*8,intent(out)::ymin, yfact, invdy
    real*8,intent(out)::xdata(nx), ydata(ny)
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    mode_read = "xy"
    if(present(mode)) then
      mode_read = mode
    end if

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    open(newunit=unit, file=trim(fname), status="old")
    do i=1,skipl
      read(unit, *)  ! read lines to be skipped
    end do
    if(mode_read=="xy") then
      do i=1,nx
        do j=1,ny
          read(unit, *) xdata(i), ydata(j), fdata(i, j)
        end do
      end do
    else if(mode_read=="yx") then
      do j=1,ny
        do i=1,nx
          read(unit, *) xdata(i), ydata(j), fdata(i, j)
        end do
      end do
    else
      print *, "ERROR: read mode in load interp_2d can be only xy or yx, found "//mode_read
      stop
    end if
    close(unit)

    if(dolog) then
      xdata(:) = log10(xdata(:))
      ydata(:) = log10(ydata(:))
      fdata(:, :) = log10(fdata(:, :))
    end if

    call check_regular_grid(xdata(:), fname)
    call check_regular_grid(ydata(:), fname)

    xmin = xdata(1)
    ymin = ydata(1)
    dx = xdata(2) - xdata(1)
    dy = ydata(2) - ydata(1)
    xfact = (nx - 1) / (xdata(nx) - xdata(1))
    yfact = (ny - 1) / (ydata(ny) - ydata(1))
    invdx = 1d0 / dx
    invdy = 1d0 / dy

  end subroutine load_data_2d

  ! ***************************
  subroutine load_3d_fit_vec(fname, nvec, nx, ny, nz, data, mode, do_log, skip_lines)
    implicit none
    integer,intent(in)::nvec, nx, ny, nz
    character(len=3),intent(in),optional::mode
    character(len=3)::mode_read
    character(len=1024),intent(in)::fname(nvec)
    type(fit3d_data_vec(nv=nvec, n1=nx, n2=ny, n3=nz)),intent(out)::data
    real*8::xdata(nx), fdata(nvec, nx, ny, nz), xmin, invdx, xfact
    real*8::ydata(ny), ymin, invdy, yfact
    real*8::zdata(nz), zmin, invdz, zfact
    integer::i
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    mode_read = "xyz"
    if(present(mode)) then
      mode_read = mode
    end if

    do i=1,nvec
      call load_data_3d(fname(i), fdata(i, :, :, :), &
        xdata, nx, xmin, xfact, invdx, &
        ydata, ny, ymin, yfact, invdy, &
        zdata, nz, zmin, zfact, invdz, mode_read, dolog, skipl)
    end do

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

    data%ymin = ymin
    data%ymax = ydata(ny)
    data%dy = ydata(2) - ydata(1)
    data%invdy = invdy
    data%yfact = yfact

    data%zmin = zmin
    data%zmax = zdata(nz)
    data%dz = zdata(2) - zdata(1)
    data%invdz = invdz
    data%zfact = zfact

  end subroutine load_3d_fit_vec

  ! ***************************
  subroutine load_3d_fit(fname, nx, ny, nz, data, mode, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    character(len=3),intent(in),optional::mode
    character(len=3)::mode_read
    integer,intent(in)::nx, ny, nz
    type(fit3d_data(n1=nx, n2=ny, n3=nz)),intent(out)::data
    real*8::xdata(nx), fdata(nx, ny, nz), xmin, invdx, xfact
    real*8::ydata(ny), ymin, invdy, yfact
    real*8::zdata(nz), zmin, invdz, zfact
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    mode_read = "xyz"
    if(present(mode)) then
      mode_read = mode
    end if

    call load_data_3d(fname, fdata, &
      xdata, nx, xmin, xfact, invdx, &
      ydata, ny, ymin, yfact, invdy, &
      zdata, nz, zmin, zfact, invdz, mode_read, dolog, skipl)

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

    data%ymin = ymin
    data%ymax = ydata(ny)
    data%dy = ydata(2) - ydata(1)
    data%invdy = invdy
    data%yfact = yfact

    data%zmin = zmin
    data%zmax = zdata(nz)
    data%dz = zdata(2) - zdata(1)
    data%invdz = invdz
    data%zfact = zfact

  end subroutine load_3d_fit

  ! ********************
  subroutine load_data_3d(fname, fdata, &
    xdata, nx, xmin, xfact, invdx, &
    ydata, ny, ymin, yfact, invdy, &
    zdata, nz, zmin, zfact, invdz, mode, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nx, ny, nz
    integer::unit, i, j, k
    real*8::dx, dy, dz
    real*8,intent(out)::fdata(nx, ny, nz), xdata(nx), xmin, xfact, invdx
    real*8,intent(out)::ydata(ny), ymin, yfact, invdy
    real*8,intent(out)::zdata(nz), zmin, zfact, invdz
    character(len=3),intent(in),optional::mode
    character(len=3)::mode_read
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    mode_read = "xyz"
    if(present(mode)) then
      mode_read = mode
    end if

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    open(newunit=unit, file=trim(fname), status="old")
    do i=1,skipl
      read(unit, *)  ! read lines to be skipped
    end do
    if(mode_read=="xyz") then
      do i=1,nx
        do j=1,ny
          do k=1,nz
            read(unit, *) xdata(i), ydata(j), zdata(k), fdata(i, j, k)
          end do
        end do
      end do
    elseif(mode_read=="zyx") then
      do k=1,nz
        do j=1,ny
          do i=1,nx
            read(unit, *) xdata(i), ydata(j), zdata(k), fdata(i, j, k)
          end do
        end do
      end do
    else
      print *, "ERROR: read mode in load interp_3d can be only xyz or zyx, found "//mode_read
      stop
    end if
    close(unit)

    if(dolog) then
      xdata(:) = log10(xdata(:))
      ydata(:) = log10(ydata(:))
      zdata(:) = log10(zdata(:))
      if(minval(fdata) <= 0d0) then
        print *, "ERROR: negative values in ", trim(fname)
        stop
      end if
      fdata(:, :, :) = log10(fdata(:, :, :))
    end if

    call check_regular_grid(xdata(:), fname, "axis: x")
    call check_regular_grid(ydata(:), fname, "axis: y")
    call check_regular_grid(zdata(:), fname, "axis: z")

    xmin = xdata(1)
    ymin = ydata(1)
    zmin = zdata(1)

    dx = xdata(2) - xdata(1)
    dy = ydata(2) - ydata(1)
    dz = zdata(2) - zdata(1)

    xfact = (nx - 1) / (xdata(nx) - xdata(1))
    yfact = (ny - 1) / (ydata(ny) - ydata(1))
    zfact = (nz - 1) / (zdata(nz) - zdata(1))

    invdx = 1d0 / dx
    invdy = 1d0 / dy
    invdz = 1d0 / dz

  end subroutine load_data_3d

  ! ***************************
  subroutine load_4d_fit_vec(fname, nvec, nx, ny, nz, nu, data, mode, do_log, skip_lines)
    implicit none
    integer,intent(in)::nvec, nx, ny, nz, nu
    character(len=4),intent(in),optional::mode
    character(len=4)::mode_read
    character(len=1024),intent(in)::fname(nvec)
    type(fit4d_data_vec(nv=nvec, n1=nx, n2=ny, n3=nz, n4=nu)),intent(out)::data
    real*8::xdata(nx), fdata(nvec, nx, ny, nz, nu), xmin, invdx, xfact
    real*8::ydata(ny), ymin, invdy, yfact
    real*8::zdata(nz), zmin, invdz, zfact
    real*8::udata(nu), umin, invdu, ufact
    integer::i
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    mode_read = "xyzu"
    if(present(mode)) then
      mode_read = mode
    end if

    do i=1,nvec
      call load_data_4d(fname(i), fdata(i, :, :, :, :), &
        xdata, nx, xmin, xfact, invdx, &
        ydata, ny, ymin, yfact, invdy, &
        zdata, nz, zmin, zfact, invdz, &
        udata, nu, umin, ufact, invdu, &
        mode_read, dolog, skipl)
    end do

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

    data%ymin = ymin
    data%ymax = ydata(ny)
    data%dy = ydata(2) - ydata(1)
    data%invdy = invdy
    data%yfact = yfact

    data%zmin = zmin
    data%zmax = zdata(nz)
    data%dz = zdata(2) - zdata(1)
    data%invdz = invdz
    data%zfact = zfact

    data%umin = umin
    data%umax = udata(nu)
    data%du = udata(2) - udata(1)
    data%invdu = invdu
    data%ufact = ufact

  end subroutine load_4d_fit_vec

  ! ***************************
  subroutine load_4d_fit(fname, nx, ny, nz, nu, data, mode, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    character(len=4),intent(in),optional::mode
    character(len=4)::mode_read
    integer,intent(in)::nx, ny, nz, nu
    type(fit4d_data(n1=nx, n2=ny, n3=nz, n4=nu)),intent(out)::data
    real*8::xdata(nx), fdata(nx, ny, nz, nu), xmin, invdx, xfact
    real*8::ydata(ny), ymin, invdy, yfact
    real*8::zdata(nz), zmin, invdz, zfact
    real*8::udata(nu), umin, invdu, ufact
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    mode_read = "xyzu"
    if(present(mode)) then
      mode_read = mode
    end if

    call load_data_4d(fname, fdata, &
      xdata, nx, xmin, xfact, invdx, &
      ydata, ny, ymin, yfact, invdy, &
      zdata, nz, zmin, zfact, invdz, &
      udata, nu, umin, ufact, invdu, mode_read, dolog, skipl)

    data%fdata = fdata
    data%xmin = xmin
    data%xmax = xdata(nx)
    data%dx = xdata(2) - xdata(1)
    data%invdx = invdx
    data%xfact = xfact

    data%ymin = ymin
    data%ymax = ydata(ny)
    data%dy = ydata(2) - ydata(1)
    data%invdy = invdy
    data%yfact = yfact

    data%zmin = zmin
    data%zmax = zdata(nz)
    data%dz = zdata(2) - zdata(1)
    data%invdz = invdz
    data%zfact = zfact

    data%umin = umin
    data%umax = udata(nu)
    data%du = udata(2) - udata(1)
    data%invdu = invdu
    data%ufact = ufact

  end subroutine load_4d_fit

  ! ********************
  subroutine load_data_4d(fname, fdata, &
    xdata, nx, xmin, xfact, invdx, &
    ydata, ny, ymin, yfact, invdy, &
    zdata, nz, zmin, zfact, invdz, &
    udata, nu, umin, ufact, invdu, &
    mode, do_log, skip_lines)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nx, ny, nz, nu
    integer::unit, i, j, k, ii
    real*8::dx, dy, dz, du
    real*8,intent(out)::fdata(nx, ny, nz, nu), xdata(nx), xmin, xfact, invdx
    real*8,intent(out)::ydata(ny), ymin, yfact, invdy
    real*8,intent(out)::zdata(nz), zmin, zfact, invdz
    real*8,intent(out)::udata(nu), umin, ufact, invdu
    character(len=4),intent(in),optional::mode
    character(len=4)::mode_read
    logical,intent(in),optional::do_log
    logical::dolog
    integer,intent(in),optional::skip_lines
    integer::skipl

    mode_read = "xyzu"
    if(present(mode)) then
      mode_read = mode
    end if

    dolog = .false.
    if(present(do_log)) then
      dolog = do_log
    end if

    skipl = 0
    if(present(skip_lines)) then
      skipl = skip_lines
    end if

    open(newunit=unit, file=trim(fname), status="old")
    do i=1,skipl
      read(unit, *)  ! read lines to be skipped
    end do
    if(mode_read=="xyzu") then
      do i=1,nx
        do j=1,ny
          do k=1,nz
            do ii=1,nu
              read(unit, *) xdata(i), ydata(j), zdata(k), udata(ii), fdata(i, j, k, ii)
            end do
          end do
        end do
      end do
    else
      print *, "ERROR: read mode in load interp_4d can be only xyzu, found "//mode_read
      stop
    end if
    close(unit)

    if(dolog) then
      xdata(:) = log10(xdata(:))
      ydata(:) = log10(ydata(:))
      zdata(:) = log10(zdata(:))
      udata(:) = log10(udata(:))
      if(minval(fdata) <= 0d0) then
        print *, "ERROR: negative values in ", trim(fname)
        stop
      end if
      fdata = log10(fdata)
    end if

    call check_regular_grid(xdata(:), fname, "axis: x")
    call check_regular_grid(ydata(:), fname, "axis: y")
    call check_regular_grid(zdata(:), fname, "axis: z")
    call check_regular_grid(udata(:), fname, "axis: u")

    xmin = xdata(1)
    ymin = ydata(1)
    zmin = zdata(1)
    umin = udata(1)

    dx = xdata(2) - xdata(1)
    dy = ydata(2) - ydata(1)
    dz = zdata(2) - zdata(1)
    du = udata(2) - udata(1)

    xfact = (nx - 1) / (xdata(nx) - xdata(1))
    yfact = (ny - 1) / (ydata(ny) - ydata(1))
    zfact = (nz - 1) / (zdata(nz) - zdata(1))
    ufact = (nu - 1) / (udata(nu) - udata(1))

    invdx = 1d0 / dx
    invdy = 1d0 / dy
    invdz = 1d0 / dz
    invdu = 1d0 / du

  end subroutine load_data_4d

  ! ********************
  subroutine check_regular_grid(xdata, fname, message)
    implicit none
    real*8,intent(in)::xdata(:)
    real*8::dx, err
    integer::n, i
    character(len=*),optional::fname
    character(len=*),optional::message
    character(len=200)::f_name

    f_name = "unknown file"
    if(present(fname)) then
      f_name = fname
    end if

    dx = xdata(2) - xdata(1)
    n = size(xdata)

    do i=3,n
      err = abs(dx - (xdata(i) - xdata(i-1))) / dx
      if(err > 1d-4) then
        print *, "ERROR: grid is not regular in "//trim(f_name)
        print *, "error (should be 0), index:", err, i
        print *, "value at i, value at i-1", xdata(i), xdata(i-1)
        print *, "dx found, expected dx", xdata(i) - xdata(i-1), dx
        if(present(message)) then
          print *, message
        end if
        stop
      end if
    end do

  end subroutine check_regular_grid

  ! ****************
  function interp_1dfit_vec(x, data, nvec, nx) result(f)
    implicit none
    real*8,intent(in)::x
    integer,intent(in)::nvec, nx
    type(fit1d_data_vec(nv=nvec, n1=nx)),intent(in)::data
    real*8::f(nvec), xlow
    integer::i

    i = floor((x - data%xmin) * data%xfact + 1)

    xlow = (i - 1) / data%xfact + data%xmin

    f = (x - xlow) * (data%fdata(:, i+1) - data%fdata(:, i)) * data%invdx + data%fdata(:, i)

  end function interp_1dfit_vec

  ! ****************
  function interp_1dfit(x, data, nx) result(f)
    implicit none
    real*8,intent(in)::x
    integer,intent(in)::nx
    type(fit1d_data(n1=nx)),intent(in)::data
    real*8::f

    f = interp_1d(x, data%fdata, data%xmin, data%xfact, data%invdx)

  end function interp_1dfit

  ! ********************
  function interp_1d(x, fdata, xmin, xfact, invdx) result(f)
    implicit none
    real*8,intent(in)::x, xmin, xfact, invdx, fdata(:)
    real*8::f, xlow
    integer::i

    i = floor((x - xmin) * xfact + 1)

    xlow = (i - 1) / xfact + xmin

    f = (x - xlow) * (fdata(i+1) - fdata(i)) * invdx + fdata(i)

  end function interp_1d

  ! ****************
  function interp_2dfit_vec(x, y, data, nvec, nx, ny) result(f)
    implicit none
    real*8,intent(in)::x, y
    integer,intent(in)::nvec, nx, ny
    type(fit2d_data_vec(nv=nvec, n1=nx, n2=ny)),intent(in)::data
    real*8::f(nvec), f0(nvec), f1(nvec), dxi, xlow, ylow
    integer::i, j

    i = floor((x - data%xmin) * data%xfact + 1)
    j = floor((y - data%ymin) * data%yfact + 1)

    xlow = (i - 1) / data%xfact + data%xmin
    ylow = (j - 1) / data%yfact + data%ymin

    dxi = (x - xlow) * data%invdx

    f0 = dxi * (data%fdata(:, i+1, j) - data%fdata(:, i, j)) + data%fdata(:, i, j)
    f1 = dxi * (data%fdata(:, i+1, j+1) - data%fdata(:, i, j+1)) + data%fdata(:, i, j+1)

    f = (y - ylow) * (f1 - f0) * data%invdy + f0

  end function interp_2dfit_vec

  ! ****************
  function interp_2dfit(x, y, data, nx, ny) result(f)
    implicit none
    real*8,intent(in)::x, y
    integer,intent(in)::nx, ny
    type(fit2d_data(n1=nx, n2=ny)),intent(in)::data
    real*8::f

    f = interp_2d(x, y, data%fdata, data%xmin, data%xfact, data%invdx, &
                  data%ymin, data%yfact, data%invdy)

  end function interp_2dfit

  ! ********************
  function interp_2d(x, y, fdata, xmin, xfact, invdx, &
    ymin, yfact, invdy) result(f)
    implicit none
    real*8,intent(in)::x, xmin, xfact, invdx, fdata(:, :)
    real*8,intent(in)::y, ymin, yfact, invdy
    real*8::f, xlow, ylow, f0, f1, dxi
    integer::i, j

    i = floor((x - xmin) * xfact + 1)
    j = floor((y - ymin) * yfact + 1)

    xlow = (i - 1) / xfact + xmin
    ylow = (j - 1) / yfact + ymin

    dxi = (x - xlow) * invdx

    f0 = dxi * (fdata(i+1, j) - fdata(i, j)) + fdata(i, j)
    f1 = dxi * (fdata(i+1, j+1) - fdata(i, j+1)) + fdata(i, j+1)

    f = (y - ylow) * (f1 - f0) * invdy + f0

  end function interp_2d

  ! ****************
  function interp_3dfit_vec(x, y, z, data, nvec, nx, ny, nz) result(f)
    implicit none
    real*8,intent(in)::x, y, z
    integer,intent(in)::nvec, nx, ny, nz
    type(fit3d_data_vec(nv=nvec, n1=nx, n2=ny, n3=nz)),intent(in)::data
    real*8::f(nvec), xlow, ylow, zlow, dxi, dyi
    real*8::f00(nvec), f10(nvec), f01(nvec), f11(nvec), f0(nvec), f1(nvec)
    integer::i, j, k

    i = floor((x - data%xmin) * data%xfact + 1)
    j = floor((y - data%ymin) * data%yfact + 1)
    k = floor((z - data%zmin) * data%zfact + 1)

    xlow = (i - 1) / data%xfact + data%xmin
    ylow = (j - 1) / data%yfact + data%ymin
    zlow = (k - 1) / data%zfact + data%zmin

    dxi = (x - xlow) * data%invdx
    dyi = (y - ylow) * data%invdy

    f00 = dxi * (data%fdata(:, i+1, j, k) - data%fdata(:, i, j, k)) + data%fdata(:, i, j, k)
    f10 = dxi * (data%fdata(:, i+1, j+1, k) - data%fdata(:, i, j+1, k)) + data%fdata(:, i, j+1, k)
    f01 = dxi * (data%fdata(:, i+1, j, k+1) - data%fdata(:, i, j, k+1)) + data%fdata(:, i, j, k+1)
    f11 = dxi * (data%fdata(:, i+1, j+1, k+1) - data%fdata(:, i, j+1, k+1)) + data%fdata(:, i, j+1, k+1)

    f0 = dyi * (f10 - f00) + f00
    f1 = dyi * (f11 - f01) + f01

    f = (z - zlow) * (f1 - f0) * data%invdz + f0

  end function interp_3dfit_vec

  ! ****************
  function interp_3dfit(x, y, z, data, nx, ny, nz) result(f)
    implicit none
    real*8,intent(in)::x, y, z
    integer,intent(in)::nx, ny, nz
    type(fit3d_data(n1=nx, n2=ny, n3=nz)),intent(in)::data
    real*8::f

    f = interp_3d(x, y, z, data%fdata, &
                  data%xmin, data%xfact, data%invdx, &
                  data%ymin, data%yfact, data%invdy, &
                  data%zmin, data%zfact, data%invdz)

  end function interp_3dfit


  ! ********************
  function interp_3d(x, y, z, fdata, &
    xmin, xfact, invdx, &
    ymin, yfact, invdy, &
    zmin, zfact, invdz) result(f)
    implicit none
    real*8,intent(in)::x, xmin, xfact, invdx, fdata(:, :, :)
    real*8,intent(in)::y, ymin, yfact, invdy
    real*8,intent(in)::z, zmin, zfact, invdz
    real*8::f, xlow, ylow, zlow, f00, f10, f01, f11, f0, f1, dxi, dyi
    integer::i, j, k

    i = floor((x - xmin) * xfact + 1)
    j = floor((y - ymin) * yfact + 1)
    k = floor((z - zmin) * zfact + 1)

    xlow = (i - 1) / xfact + xmin
    ylow = (j - 1) / yfact + ymin
    zlow = (k - 1) / zfact + zmin

    dxi = (x - xlow) * invdx
    dyi = (y - ylow) * invdy

    f00 = dxi * (fdata(i+1, j, k) - fdata(i, j, k)) + fdata(i, j, k)
    f10 = dxi * (fdata(i+1, j+1, k) - fdata(i, j+1, k)) + fdata(i, j+1, k)
    f01 = dxi * (fdata(i+1, j, k+1) - fdata(i, j, k+1)) + fdata(i, j, k+1)
    f11 = dxi * (fdata(i+1, j+1, k+1) - fdata(i, j+1, k+1)) + fdata(i, j+1, k+1)

    f0 = dyi * (f10 - f00) + f00
    f1 = dyi * (f11 - f01) + f01

    f = (z - zlow) * (f1 - f0) * invdz + f0

  end function interp_3d

  ! ****************
  function interp_4dfit_vec(x, y, z, u, data, nvec, nx, ny, nz, nu) result(f)
    implicit none
    real*8,intent(in)::x, y, z, u
    integer,intent(in)::nvec, nx, ny, nz, nu
    type(fit4d_data_vec(nv=nvec, n1=nx, n2=ny, n3=nz, n4=nu)),intent(in)::data
    real*8::f(nvec), xlow, ylow, zlow, ulow, dxi, dyi, dzi
    real*8::f00(nvec), f10(nvec), f01(nvec), f11(nvec), f0(nvec), f1(nvec)
    real*8::f000(nvec), f010(nvec), f001(nvec), f011(nvec)
    real*8::f100(nvec), f110(nvec), f101(nvec), f111(nvec)
    integer::i, j, k, ii

    i = floor((x - data%xmin) * data%xfact + 1)
    j = floor((y - data%ymin) * data%yfact + 1)
    k = floor((z - data%zmin) * data%zfact + 1)
    ii = floor((u - data%umin) * data%ufact + 1)

    xlow = (i - 1) / data%xfact + data%xmin
    ylow = (j - 1) / data%yfact + data%ymin
    zlow = (k - 1) / data%zfact + data%zmin
    ulow = (ii - 1) / data%ufact + data%umin

    dxi = (x - xlow) * data%invdx
    dyi = (y - ylow) * data%invdy
    dzi = (z - zlow) * data%invdz

    f000 = dxi * (data%fdata(:, i+1, j, k, ii) - data%fdata(:, i, j, k, ii)) + data%fdata(:, i, j, k, ii)
    f100 = dxi * (data%fdata(:, i+1, j+1, k, ii) - data%fdata(:, i, j+1, k, ii)) + data%fdata(:, i, j+1, k, ii)
    f010 = dxi * (data%fdata(:, i+1, j, k+1, ii) - data%fdata(:, i, j, k+1, ii)) + data%fdata(:, i, j, k+1, ii)
    f110 = dxi * (data%fdata(:, i+1, j+1, k+1, ii) - data%fdata(:, i, j+1, k+1, ii)) + data%fdata(:, i, j+1, k+1, ii)

    f001 = dxi * (data%fdata(:, i+1, j, k, ii+1) - data%fdata(:, i, j, k, ii+1)) + data%fdata(:, i, j, k, ii+1)
    f101 = dxi * (data%fdata(:, i+1, j+1, k, ii+1) - data%fdata(:, i, j+1, k, ii+1)) + data%fdata(:, i, j+1, k, ii+1)
    f011 = dxi * (data%fdata(:, i+1, j, k+1, ii+1) - data%fdata(:, i, j, k+1, ii+1)) + data%fdata(:, i, j, k+1, ii+1)
    f111 = dxi * (data%fdata(:, i+1, j+1, k+1, ii+1) - data%fdata(:, i, j+1, k+1, ii+1)) + data%fdata(:, i, j+1, k+1, ii+1)

    f00 = dyi * (f100 - f000) + f000
    f10 = dyi * (f110 - f010) + f010

    f01 = dyi * (f100 - f001) + f001
    f11 = dyi * (f111 - f011) + f011

    f0 = dzi * (f10 - f00) + f00
    f1 = dzi * (f11 - f01) + f01

    f = (u - ulow) * (f1 - f0) * data%invdu + f0

  end function interp_4dfit_vec


end module prizmo_fit
