program test
  use prizmo
  use ifport
  integer,parameter::nt=50, nxpos=100, nstep=20
  real*8::ngas_grid(nt, nxpos), Tgas_grid(nt, nxpos)
  real*8::x(nspecies), Tgas, Tdust, jflux(prizmo_nphoto), xall(nspecies, nt, nxpos), jflux0(prizmo_nphoto)
  real*8::r, dt, spy, dr, au2cm, t, ngas, dz, dt_ref
  real*8::rad_Ncol_H2, rad_Ncol_CO, vert_Ncol_CO(nxpos), zold(nxpos)
  real*8::cools(5), heats(4), energy_ev(prizmo_nphoto), Jscale, krate(prizmo_nreactions)
  real*8::coola(nspecies), xpos, zpos, xpos_max, xpos_min, theta_min, theta_max, theta
  real*8::rold, ngas_old, hscale, rstar, gravity, mstar
  real::time_start, time_stop, tcpu(nt, nxpos), total_time
  integer::i, it, ix, istep, ierr
  character(len=10) :: sstep

  ! coordinates for the 2D disk are
  ! theta: polar angle, 0 is the midplane, pi/2 is the pole
  ! r: radius from the star, polar coordinates
  ! z: vertical cartesian coordinate
  ! x: horizontal cartesian coordinate

  spy = 365 * 3600 * 24.  ! seconds per year
  au2cm = 1.49598e13  ! 1 au in cm
  rstar = 6.957d10  ! 1 Rsun in cm
  mstar = 1.989d33  ! 1 Msun in g
  gravity = 6.67259d-8  ! gravitational constant cgs

  ! limits, cm
  xpos_min = 1d0 * au2cm  ! inner disk rim
  xpos_max = 1d2 * au2cm  ! outer disk rim
  theta_min = 0  ! midplane
  theta_max = pi / 4

  ! initialize prizmo, mandatory, only once
  call prizmo_init()

  ! initialize cpu time to zero
  tcpu = 0e0
  total_time = 0e0


  ! interates to find equilibrium
  do istep=1,nstep

    ! increase timestep
    dt_ref = spy * 1d1**((istep - 1) * (6. + 4) / (nstep - 1) - 4)

    ! first step is just chemistry, then solve also cooling/heating
    if(istep == 1) then
      call prizmo_set_solve_thermo(.false.)
      call prizmo_set_solve_chemistry(.true.)
    else
      call prizmo_set_solve_thermo(.true.)
      call prizmo_set_solve_chemistry(.true.)
    end if

    ! step number in string format
    write(sstep,"(I0.5)") istep

    ! open files for ouput
    open(23, file="output_chem_"//trim(sstep)//".dat", status="replace")
    open(24, file="output_cool_"//trim(sstep)//".dat", status="replace")
    open(25, file="output_heat_"//trim(sstep)//".dat", status="replace")
    open(26, file="output_tcpu_"//trim(sstep)//".dat", status="replace")
    open(55, file="output_ngas_"//trim(sstep)//".dat", status="replace")

    ! loop on angles to prepare the disk structure
    do it=1,nt
      theta = (it - 1) * (theta_max - theta_min) / (nt  - 1) + theta_min
      ! loop on horizontal position
      do ix=1,nxpos
        xpos = 1d1**((ix - 1) * (log10(xpos_max) - log10(xpos_min)) / (nxpos  - 1) + log10(xpos_min))
        r = xpos * sqrt(1d0 + tan(theta)**2)  ! radius
        zpos = xpos * tan(theta)  ! vertical coordinate
        Tgas = 200. * (xpos / au2cm)**(-0.5)  ! midplane temperature
        ! temperature is initialized only in the first step, then computed
        if(istep == 1) then
          Tgas_grid(it, ix) = Tgas
        end if

        ! compute scale height and density
        hscale = sqrt(kboltzmann * Tgas / 2.34 / pmass) / sqrt(gravity * Mstar / xpos**3)
        ngas = 1700.* (xpos / au2cm)**(-1.5) / hscale / sqrt(2*pi) * exp(-zpos**2 / 2. / hscale**2) / pmass / 2.34
        ngas = max(ngas, 1d1)  ! disk envelope density
        ! first step is to initialize temperature
        if (istep==1) then
          ngas_grid(it, ix) = ngas
        !else
          !ngas_grid(it, ix) = ngas_grid(it, ix) + (ngas - ngas_grid(it, ix)) * 1d-4
        end if

        ! initialize chemistry during the first step
        if(istep == 1) then
          x(:) = 0d0
          x(prizmo_idx_H2) = 0d0 !ngas / 2d0
          x(prizmo_idx_H) = ngas
          x(prizmo_idx_C) = ngas * 1d-4
          x(prizmo_idx_O) = ngas * 3d-4
          x(prizmo_idx_He) = ngas * 1d-1
          x(prizmo_idx_E) = x(prizmo_idx_Cj) + x(prizmo_idx_Hj) + x(prizmo_idx_Oj) + x(prizmo_idx_Hej)
          xall(:, it, ix) = x / sum(x) * ngas
        end if

        ! save disk structure to disk :)
        write(55, '(99e17.8e3)') theta, r, ngas_grid(it, ix), Tgas_grid(it, ix), zpos, hscale
      end do
    end do
    close(55)


    ! compute initial z coordinate for each x grid point at theta_max
    do ix=1,nxpos
      xpos = 1d1**((ix - 1) * (log10(xpos_max) - log10(xpos_min)) / (nxpos  - 1) + log10(xpos_min))
      zold(ix) = xpos * tan(theta_max)
    end do

    jflux0 = prizmo_load_radiation_field("runtime_data/radiation_field.dat")

    ! integration is from pole to midplane
    do it=nt,1,-1
      print '(a10,3I5,e17.8e3)', "*********", istep, nt-it+1, nt, dt_ref / spy
      theta = (it - 1) * (theta_max - theta_min) / (nt - 1) + theta_min

      ! set dust/gas mass ratio
      call prizmo_set_d2g(1d-2)

      ! compute min radius at given theta
      r = xpos_min * sqrt(1d0 + tan(theta)**2)
      rold = r

      ngas_old = ngas_grid(it, 1)

      ! init vertical CO column density, 1/cm2
      vert_Ncol_CO(:) = 0d0

      ! init radial column densities
      rad_Ncol_H2 = 0d0  ! cm-2
      rad_Ncol_CO = 0d0  ! cm-2

      ! load radiation from file
      Jscale = 1d0
      jflux(:) = 4d0 * pi**2 * Jscale * jflux0 * rstar**2 / r**2

      ! loop on the x grid points
      do ix=1,nxpos
        ngas = ngas_grid(it, ix)

        ! change absolute tolerance depending on the total density
        call prizmo_set_atol_all(ngas * 1d-25)

        xpos = 1d1**((ix - 1) * (log10(xpos_max) - log10(xpos_min)) / (nxpos  - 1) + log10(xpos_min))
        r = xpos * sqrt(1d0 + tan(theta)**2)
        zpos = sqrt(r**2 - xpos**2)

        ! get abundances and temperature from the grid
        x = xall(:, it, ix)
        Tgas = Tgas_grid(it, ix)

        ! compute cell spacing
        dr = r - rold
        dz = max(zold(ix) - zpos, 1d-40)

        ! scale radtiation with distance
        jflux = jflux * rold**2 / r**2

        ! update column densities
        rad_Ncol_H2 = rad_Ncol_H2 + x(prizmo_idx_H2) * dr
        rad_Ncol_CO = rad_Ncol_CO + x(prizmo_idx_CO) * dr
        vert_Ncol_CO(ix) = vert_Ncol_CO(ix) + x(prizmo_idx_CO) * dz

        ! print '(a5,I5,99e17.8e3)', "eeee", ix, Tgas, ngas, sum(x), abs(sum(cools) - sum(heats)) / abs(sum(cools))

        ! set column densities
        call prizmo_set_radial_Ncol_H2(rad_Ncol_H2)
        call prizmo_set_radial_Ncol_CO(rad_Ncol_CO)
        call prizmo_set_vertical_Ncol_CO(vert_Ncol_CO(ix))
        t = 0d0
        dt = dt_ref

        ! evolve thermochemistry
        call cpu_time(time_start)
        ! if(prizmo_is_thermal_equilibrium(x, Tgas, jflux, dt)) then
        !   print *, "eq"
        ! else
        ! end if
        call prizmo_evolve(x, Tgas, jflux, dt, 0, ierr)
        call cpu_time(time_stop)
        total_time = total_time + time_stop - time_start
        tcpu(it, ix) = tcpu(it, ix) + time_stop - time_start

        ! do RT
        call prizmo_rt(x, Tgas, jflux, dr)

        ! put quantities back to the grid
        xall(:, it, ix) = x
        Tgas_grid(it, ix) = Tgas

        rold = r
        zold(ix) = zpos
        ngas_old = ngas

        !krate = prizmo_get_rates(x, Tgas, Tdust, jflux)

        !this is for output
        Tdust = prizmo_get_tdust(x, Tgas, jflux)
        cools = prizmo_get_cooling_array(x, Tgas, Tdust, jflux)
        heats = prizmo_get_heating_array(x, Tgas, Tdust, jflux)
                ! ! write output to file
        write(23, '(99e17.8e3)') theta, r, Tgas, Tdust, x
        write(24, '(99e17.8e3)') theta, r, Tgas, cools
        write(25, '(99e17.8e3)') theta, r, Tgas, heats
        write(26, '(99e17.8e3)') theta, r, Tgas, ngas, tcpu(it, ix)!, prizmo_get_chi_FUV()

      end do

    end do
    print *, "total time (s):", total_time

    close(23)
    close(24)
    close(25)
    close(26)

  end do

  print *, "done, KTHXBYE!"

end program test
