program test
  use prizmo
  real*8::x(nspecies), Tgas, Tdust, jflux(prizmo_nphoto)
  real*8::r, dt, spy, dr, au2cm, t, ngas, N2Av, Av
  real*8::rad_Ncol_H2, rad_Ncol_CO, rad_Ncol_tot
  real*8::cools(5), heats(4), energy_ev(prizmo_nphoto), Jscale, krate(prizmo_nreactions)
  real*8::coola(nspecies)
  integer::i, ixi

  N2Av = 6.289d-22
  spy = 365 * 3600 * 24.
  au2cm = 1.49598e13

  call prizmo_init()

  do ixi=1,10
    ngas = 1d1**3 !1d3
    Jscale = 1d1**((ixi - 1) * (6.) / (10 - 1) - 8.)

    r = 1d18 / ngas  !1d- * au2cm  ! cm

    jflux(:) = 4d0 * pi * Jscale * prizmo_load_radiation_field("runtime_data/radiation_field.dat")

    energy_ev = prizmo_get_energy_ev()

    call prizmo_set_d2g(1d-10)

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
    rad_Ncol_tot = 0d0  !cm-2

    dt = spy * 1d8
    dr = r / 1d1
    t = 0d0
    do
      dr = dr * 1.2
      !dt = dt * 1.2
      r = r + dr
      t = t + dt
      Av = r * ngas * N2av
      rad_Ncol_H2 = rad_Ncol_H2 + x(prizmo_idx_H2) * dr
      rad_Ncol_CO = rad_Ncol_CO + x(prizmo_idx_CO) * dr
      rad_Ncol_tot = rad_Ncol_tot + ngas * dr

      !call prizmo_set_user_Av(Av)
      call prizmo_set_radial_Ncol_H2(rad_Ncol_H2)
      call prizmo_set_radial_Ncol_CO(rad_Ncol_CO)
      call prizmo_set_vertical_Ncol_CO(rad_Ncol_CO)
      call prizmo_evolve(x, Tgas, Tdust, jflux, dt)
      call prizmo_rt(x, Tgas, Tdust, jflux, dr)

      Tdust = prizmo_get_tdust(x, Tgas, jflux)

      cools = prizmo_get_cooling_array(x, Tgas, Tdust, jflux)
      heats = prizmo_get_heating_array(x, Tgas, Tdust, jflux)

      coola = prizmo_get_atomic_cooling_array(x, Tgas)

      krate = prizmo_get_rates(x, Tgas, Tdust, jflux)

      if(rad_Ncol_tot<=5d21) then
        call prizmo_save_cooling_function(x, jflux, "function_cooling_5d21.dat")
        call prizmo_save_heating_function(x, jflux, "function_heating_5d21.dat")
      end if

      if(rad_Ncol_tot<=1d22) then
        call prizmo_save_cooling_function(x, jflux, "function_cooling_1d22.dat")
        call prizmo_save_heating_function(x, jflux, "function_heating_1d22.dat")
      end if

      if(rad_Ncol_tot<=2d22) then
        call prizmo_save_cooling_function(x, jflux, "function_cooling_2d22.dat")
        call prizmo_save_heating_function(x, jflux, "function_heating_2d22.dat")
      end if

      print '(a30,999e17.8e3)', "r/au, Ncol/cm-2, Tgas/K, kmin, kmax", r/au2cm, rad_Ncol_tot, Tgas, minval(krate), maxval(krate)
      write(23, '(99e17.8e3)') Jscale, Av, rad_Ncol_tot, Tgas, Tdust, x
      write(24, '(99e17.8e3)') Jscale, Av, rad_Ncol_tot, Tgas, cools
      write(25, '(99e17.8e3)') Jscale, Av, rad_Ncol_tot, Tgas, heats
      write(27, '(99e17.8e3)') Jscale, Av, rad_Ncol_tot, Tgas, krate(260), krate(261), krate(263)
      write(28, '(99e17.8e3)') Jscale, Av, rad_Ncol_tot, Tgas, coola(:)
      do i=1,nphoto
        write(26, '(99E17.8e3)') Jscale, Av, rad_Ncol_tot, energy_ev(i), jflux(i)
      end do
      write(26, *)
      if(rad_Ncol_tot > 2d22) exit
    end do
    !if(t > spy * 1d8) exit
  end do

end program test
