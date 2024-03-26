module prizmo_heating_CR
  use prizmo_commons
contains

  ! ***************
  function heating_CR(x) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::heat

    !! PREPROCESS_CR_HEATING
    !! PREPROCESS_END

  end function heating_CR

end module prizmo_heating_CR
