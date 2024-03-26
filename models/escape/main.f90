program test
  use prizmo
  real*8::x(nspecies), Tgas, jflux(prizmo_nphoto)
  real*8::dt, spy, t, ngas
  real*8::rad_Ncol_H2, rad_Ncol_CO

  spy = 365 * 3600 * 24.

  call prizmo_init()

  jflux(:) = 2d0 * pi * prizmo_load_radiation_field("runtime_data/radiation_field.dat") * 1d-10! * rstar**2 / r**2

  call prizmo_set_d2g(1d-2)

  x(:) = 0d0
  Tgas = 5d3
  ngas = 1d3
  x(prizmo_idx_H2) = ngas / 2d0

  rad_Ncol_H2 = 0d0  ! cm-2
  rad_Ncol_CO = 0d0  ! cm-2

  call prizmo_set_radial_Ncol_H2(rad_Ncol_H2)
  call prizmo_set_radial_Ncol_CO(rad_Ncol_CO)
  call prizmo_set_vertical_Ncol_CO(rad_Ncol_CO)

  dt = spy * 1d-1
  t = 0d0
  do
    dt = dt * 1.1
    call prizmo_evolve(x, Tgas, jflux, dt)
    !call prizmo_rt(x, Tgas, jflux, dr)
    write(22, '(99e17.8e3)') t / spy, Tgas, x
    t = t + dt
    if(t > 1d4 * spy) exit
  end do

  print *, "done, bye!"

end program test
