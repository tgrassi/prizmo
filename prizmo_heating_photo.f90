module prizmo_heating_photo
  use prizmo_commons
contains

  ! ***************
  function heating_photo(x, Tgas, Tdust, ntot) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, ntot
    real*8::heat

    heat = 0d0

    !! PREPROCESS_PHOTOHEATING
    !! PREPROCESS_END

  end function heating_photo

end module prizmo_heating_photo
