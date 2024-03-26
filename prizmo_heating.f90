module prizmo_heating
  use prizmo_commons
  use prizmo_heating_photo
  use prizmo_heating_photoelectric
  use prizmo_heating_CR
  !! PREPROCESS_USE_H2DISS
  !! PREPROCESS_END
  integer,parameter::nheating=4
contains

  ! **************************
  function heating(xin, Tgas, Tdust, jflux) result(heat)
    implicit none
    real*8,intent(in)::xin(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::heat, heats(nheating), x(nspecies)

    x = max(xin, 0d0)

    heats = heating_array(x, Tgas, Tdust, jflux)

    heat = sum(heats)

  end function heating

  ! **************************
  function heating_array(x, Tgas, Tdust, jflux) result(heats)
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::heats(nheating), ntot

    ntot = sum(x)

    heats(1) = heating_photo(x, Tgas, Tdust, jflux, ntot)
    heats(2) = heating_photoelectric(x, Tgas, jflux)
    heats(3) = heating_CR(x)
    !! PREPROCESS_H2DISS_HEATING
    !! PREPROCESS_END

  end function heating_array

end module prizmo_heating
