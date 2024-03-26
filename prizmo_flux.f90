module prizmo_flux
  use prizmo_commons
  use prizmo_rates
contains

  ! *****************
  function get_flux(x, Tgas, Tdust) result(flux)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust
    real*8::flux(nreactions)

    call compute_rates(x, Tgas, Tdust)

    !! PREPROCESS_FLUXES
    !! PREPROCESS_END

  end function get_flux

end module prizmo_flux
