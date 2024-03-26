program test
  use prizmo
  use ifport
  integer,parameter::nt=70, nr=400, nstep=10
  real*8::ngas_grid(nt, nr), Tgas_grid(nt, nr), thetas(nt, nr), rads(nt, nr)
  real*8::x(nspecies), Tgas, Tdust, jflux(prizmo_nphoto), xall(nspecies, nt, nr)
  real*8::r, dt, spy, dr, au2cm, t, ngas, dz, dt_ref
  real*8::rad_Ncol_H2, rad_Ncol_CO, vert_Ncol_CO(nr)
  real*8::cools(5), heats(4), energy_ev(prizmo_nphoto), Jscale, krate(prizmo_nreactions)
  real*8::coola(nspecies), xpos, zpos, xpos_max, xpos_min, theta_min, theta_max, theta
  real*8::rold, ngas_old, hscale, rstar, gravity, mstar
  real::time_start, time_stop, tcpu(nt, nr)
  integer::i, it, ir, istep
  character(len=10) :: sstep

  ! coordinates for the 2D disk are
  ! theta: polar angle, 0 is the midplane, pi/2 is the pole
  ! r: radius from the star, polar coordinates
  ! z: vertical cartesian coordinate
  ! x: horizontal cartesian coordinate

  spy = 365 * 3600 * 24.  ! seconds per year
  au2cm = 1.49598e13  ! 1 au in cm
  rstar = 6.957d10  ! 1 Rsun in cm

  ! initialize prizmo, mandatory, only once
  call prizmo_init()

  open(77, file="dump.dat", status="old")
  do it=1, nt
    do ir=1, nr
      read(77, *) thetas(it, ir), rads(it, ir), Tgas_grid(it, ir), ngas_grid(it, ir), xall(:, it, ir)
      xall(:, it, ir) = prizmo_frac2n(xall(:, it, ir), ngas_grid(it, ir))
      ngas_grid(it, ir) = ngas_grid(it, ir) / pmass / 2.
    end do
  end do
  close(77)

  print *, minval(ngas_grid)
  t = 0d0

  ! interates to find equilibrium
  do istep=1,nstep

    tcpu = 0d0

    ! increase timestep
    dt_ref = spy * 1d1**((istep - 1) * (6. + 1.) / (nstep - 1) - 1.)

    ! first step is just chemistry, then solve also cooling/heating
    if(istep == 1) then
       call prizmo_set_solve_thermo(.false.)
       call prizmo_set_solve_chemistry(.false.)
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
    open(27, file="output_atomic_"//trim(sstep)//".dat", status="replace")
    open(55, file="output_ngas_"//trim(sstep)//".dat", status="replace")

    do it=1, nt

      print *, istep, it, nt

      !vert_Ncol_CO(:) = 0d0

      ! init radial column densities
      rad_Ncol_H2 = 0d0  ! cm-2
      rad_Ncol_CO = 0d0  ! cm-2

      rold = rstar
      r = rold

      ! load radiation from file
      Jscale = 1d0
      jflux(:) = 4d0 * pi**2 * Jscale * prizmo_load_radiation_field("runtime_data/radiation_field.dat") * rstar**2 / r**2

      !theta = thetas(it, 1)

      !if(theta < .18*pi) cycle

      do ir=1,nr
        ngas = ngas_grid(it, ir)

        ! set dust/gas mass ratio
        call prizmo_set_d2g(1d-2)

        ! change absolute tolerance depending on the total density
        call prizmo_set_atol_all(ngas * 1d-25)

        ! get abundances and temperature from the grid
        x = xall(:, it, ir)
        Tgas = Tgas_grid(it, ir)

        r = rads(it, ir)
        theta = thetas(it, ir)

        !if(r > 2d14) exit

        ! compute cell spacing
        dr = r - rold

        ! scale radtiation with distance
        jflux = jflux * rold**2 / r**2

        ! update column densities
        rad_Ncol_H2 = rad_Ncol_H2 + x(prizmo_idx_H2) * dr
        rad_Ncol_CO = rad_Ncol_CO + x(prizmo_idx_CO) * dr
        !vert_Ncol_CO(ix) = vert_Ncol_CO(ix) + x(prizmo_idx_CO) * dz

        ! print '(a5,I5,99e17.8e3)', "eeee", ix, Tgas, ngas, sum(x), abs(sum(cools) - sum(heats)) / abs(sum(cools))

        ! set column densities
        call prizmo_set_radial_Ncol_H2(rad_Ncol_H2)
        call prizmo_set_radial_Ncol_CO(rad_Ncol_CO)
        call prizmo_set_vertical_Ncol_CO(rad_Ncol_CO)  ! FIXME
        dt = dt_ref
        t = t + dt

        ! evolve thermochemistry
        call cpu_time(time_start)
        call prizmo_evolve(x, Tgas, jflux, dt)
        call cpu_time(time_stop)

        tcpu(it, ir) = tcpu(it, ir) + time_stop - time_start

        ! do RT
        call prizmo_rt(x, Tgas, jflux, dr)

        ! set quantities back to the grid
        xall(:, it, ir) = x
        Tgas_grid(it, ir) = Tgas

        rold = r

        ! this is for output
        Tdust = prizmo_get_tdust(x, Tgas, jflux)
        cools = prizmo_get_cooling_array(x, Tgas, Tdust, jflux)
        heats = prizmo_get_heating_array(x, Tgas, Tdust, jflux)

        ! write output to file
        write(23, '(99e17.8e3)') theta, r, Tgas, Tdust, x
        write(24, '(99e17.8e3)') theta, r, Tgas, cools
        write(25, '(99e17.8e3)') theta, r, Tgas, heats
        write(26, '(99e17.8e3)') theta, r, Tgas, ngas, tcpu(it, ir), prizmo_get_chi_FUV(jflux), t
        write(27, '(99e17.8e3)') theta, r, Tgas, prizmo_get_atomic_cooling_array(x, Tgas)

      end do

    end do

    close(23)
    close(24)
    close(25)
    close(26)
    close(27)

  end do

  print *, "done, KTHXBYE!"

end program test
