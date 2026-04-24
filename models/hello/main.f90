program test
  use prizmo
  implicit none
  integer::ierr
  real*8::x(nspecies), Tgas, jflux(nphoto), dt, seconds_per_year


  ! initialize prizmo, mandatory, only once
  call prizmo_init()

  x(:) = 0d0
  x(prizmo_idx_Cj) = 1d0
  x(prizmo_idx_E) = 1d0

  Tgas = 1d2
  jflux(:) = 0d0
  seconds_per_year = 365. * 24. * 3600
  dt = 1d5 * seconds_per_year

  call prizmo_evolve(x, Tgas, jflux, dt, 0, ierr)

  print *, "C ", x(prizmo_idx_C)
  print *, "C+", x(prizmo_idx_Cj)
  print *, "e-", x(prizmo_idx_E)

  print *, "done, KTHXBYE!"

end program test
