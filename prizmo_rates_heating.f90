module prizmo_rates_heating
  use prizmo_commons
  use prizmo_shielding
contains

  ! ************************
  subroutine compute_photorates_heating(x, Tgas, jflux)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto)
    real*8::f(nphoto), kernel(nphoto)
    real*8::log_NH2, log_NCO

    kall_heat = 0d0

    log_NH2 = log10(radial_Ncol_H2 + 1d-40)
    log_NCO = log10(radial_Ncol_CO + 1d-40)

    kernel = jflux / energy / hplanck

    !! PREPROCESS_PHOTOHEATING_RATE
    !! PREPROCESS_END

  end subroutine compute_photorates_heating

end module prizmo_rates_heating
