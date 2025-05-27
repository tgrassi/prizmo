module fit_types
    implicit none

    ! Note on the ifdef statements for different compilers:
    ! Intel compiler initializes the arrays in the types by using the integer,len attribute.
    ! This is not standard Fortran, but it works with Intel compilers.
    ! GNU compilers do not support this, so we use allocatable arrays instead.
    ! GNU compilers do not support the integer,len attribute in types, so we use allocatable arrays instead.
    ! This requires also the allocation of the arrays in the code via allocate().

! INTEL COMPILER ALTERNATIVE
#ifdef __INTEL_COMPILER
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

#endif

! GNU COMPILER ALTERNATIVE
#ifdef __GNUC__
    type, public :: fit1d_data
        real*8,allocatable::fdata(:)
        real*8::xmin, invdx, dx, xmax, xfact
    end type

    type, public :: fit2d_data
        real*8,allocatable::fdata(:, :)
        real*8::xmin, invdx, dx, xmax, xfact
        real*8::ymin, invdy, dy, ymax, yfact
    end type

    type, public :: fit3d_data
        integer::n1, n2, n3
        real*8,allocatable::fdata(:, :, :)
        real*8::xmin, invdx, dx, xmax, xfact
        real*8::ymin, invdy, dy, ymax, yfact
        real*8::zmin, invdz, dz, zmax, zfact
    end type

    type, public :: fit4d_data
        real*8,allocatable::fdata(:, :, :, :)
        real*8::xmin, invdx, dx, xmax, xfact
        real*8::ymin, invdy, dy, ymax, yfact
        real*8::zmin, invdz, dz, zmax, zfact
        real*8::umin, invdu, du, umax, ufact
    end type

    type, public :: fit1d_data_vec
        integer::nv, n1
        real*8,allocatable::fdata(:, :)
        real*8::xmin, invdx, dx, xmax, xfact
    end type

    type, public :: fit2d_data_vec
        integer::nv, n1, n2
        real*8,allocatable::fdata(:, :, :)
        real*8::xmin, invdx, dx, xmax, xfact
        real*8::ymin, invdy, dy, ymax, yfact
    end type

    type, public :: fit3d_data_vec
        integer::nv, n1, n2, n3
        real*8,allocatable::fdata(:, :, :, :)
        real*8::xmin, invdx, dx, xmax, xfact
        real*8::ymin, invdy, dy, ymax, yfact
        real*8::zmin, invdz, dz, zmax, zfact
    end type

    type, public :: fit4d_data_vec
        integer::nv, n1, n2, n3, n4
        real*8,allocatable::fdata(:, :, :, :, :)
        real*8::xmin, invdx, dx, xmax, xfact
        real*8::ymin, invdy, dy, ymax, yfact
        real*8::zmin, invdz, dz, zmax, zfact
        real*8::umin, invdu, du, umax, ufact
    end type
#endif


end module fit_types