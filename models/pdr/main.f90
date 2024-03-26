program test
  use prizmo
  real*8::x(nspecies), Tgas, Tdust, jflux(prizmo_nphoto)
  real*8::rstar, r, dt, spy, dr, au2cm, t, ngas, N2Av, Av
  real*8::rad_Ncol_H2, rad_Ncol_CO
  real*8::cools(5), heats(4), energy_ev(prizmo_nphoto), Jscale, krate(prizmo_nreactions)
  integer::i

  N2Av = 6.289d-22
  spy = 365 * 3600 * 24.
  au2cm = 1.49598e13

  ngas = 1d1**3 !1d3
  Jscale = 1d1

  rstar = 0.00465047d0 * au2cm  ! cm
  r = 1d-6 / N2Av / ngas  !1d- * au2cm  ! cm

  call prizmo_init()

  jflux(:) = 2d0 * pi * Jscale * prizmo_load_radiation_field("runtime_data/radiation_field.dat")! * rstar**2 / r**2

  energy_ev = prizmo_get_energy_ev()

  call prizmo_set_d2g(1d-2)

  x(:) = 0d0
  Tgas = 5d3
  Tdust = 2d1
  x(prizmo_idx_H2) = ngas / 2d0
  x(prizmo_idx_Hj) = 0d0
  x(prizmo_idx_E) = 0d0
  x(prizmo_idx_C) = ngas * 1d-4
  x(prizmo_idx_O) = ngas * 3d-4
  x(prizmo_idx_He) = ngas * 1d-1

  rad_Ncol_H2 = 0d0  ! cm-2
  rad_Ncol_CO = 0d0  ! cm-2

  dt = spy * 1d8
  dr = r / 1d1
  t = 0d0
  do
    dr = dr * 1.2
    !dt = dt * 1.2
    r = r + dr
    t = t + dt
    Av = r * ngas * N2av

    !call prizmo_set_user_Av(Av)
    call prizmo_set_radial_Ncol_H2(rad_Ncol_H2)
    call prizmo_set_radial_Ncol_CO(rad_Ncol_CO)
    call prizmo_set_vertical_Ncol_CO(rad_Ncol_CO)
    call prizmo_evolve(x, Tgas, jflux, dt)
    call prizmo_rt(x, Tgas, jflux, dr)

    Tdust = prizmo_get_tdust(x, Tgas, jflux)

    cools = prizmo_get_cooling_array(x, Tgas, Tdust, jflux)
    heats = prizmo_get_heating_array(x, Tgas, Tdust, jflux)

    krate = prizmo_get_rates(x, Tgas, Tdust, jflux)

    if(Av<1d-4) then
      call prizmo_save_cooling_function(x, jflux, "function_cooling_av1m4.dat")
      call prizmo_save_heating_function(x, jflux, "function_heating_av1m4.dat")
    end if

    if(Av<=1d0) then
      call prizmo_save_cooling_function(x, jflux, "function_cooling_av1m0.dat")
      call prizmo_save_heating_function(x, jflux, "function_heating_av1m0.dat")
    end if

    print '(a20,999e17.8e3)', "Av, NCO, NH2, Tgas", Av, rad_Ncol_CO, rad_Ncol_H2, Tgas
    write(23, '(99e17.8e3)') Av, t / spy, Tgas, Tdust, x
    write(24, '(99e17.8e3)') Av, t / spy, Tgas, cools
    write(25, '(99e17.8e3)') Av, t / spy, Tgas, heats
    write(27, '(99e17.8e3)') Av, t / spy, Tgas, krate(269), krate(270), krate(268)
    write(28, '(99e17.8e3)') Av, t / spy, Tgas, prizmo_get_atomic_cooling_array(x, tgas)
    do i=1,nphoto
      write(26, '(99E17.8e3)') Av, t / spy, energy_ev(i), jflux(i)
    end do
    write(26, *)

    rad_Ncol_H2 = rad_Ncol_H2 + x(prizmo_idx_H2) * dr
    rad_Ncol_CO = rad_Ncol_CO + x(prizmo_idx_CO) * dr

    if(Av > 3d1) exit
    !if(t > spy * 1d8) exit
  end do

end program test
