module prizmo_rates_heating
  use prizmo_commons
contains

  ! ************************
  subroutine compute_photorates_heating(x, Tgas, jflux)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto)
    real*8::f(nphoto), kernel(nphoto)

    kall_heat = 0d0

    kernel = jflux / energy / hplanck

    !! PREPROCESS_PHOTOHEATING_RATE
    !! PREPROCESS_END

  end subroutine compute_photorates_heating

end module prizmo_rates_heating
