module prizmo_rates_photo
  use prizmo_commons
  use prizmo_shielding
contains

  ! ************************
  subroutine compute_photorates(x, Tgas, jflux)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto)
    real*8::f(nphoto), log_tgas, log_NH2, log_NCO, kernel(nphoto)

    log_tgas = log10(Tgas)
    log_NH2 = log10(radial_Ncol_H2 + 1d-40)
    log_NCO = log10(radial_Ncol_CO + 1d-40)

    kernel = jflux / energy / hplanck

    !! PREPROCESS_PHOTORATES
    !! PREPROCESS_END

  end subroutine compute_photorates

end module prizmo_rates_photo
