module prizmo_cooling_custom
  use prizmo_commons
contains

  ! ***************
  function cooling_custom(x, Tgas, Tdust) result(cool)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust
    real*8::cool

    cool = 0d0

    !! PREPROCESS_CUSTOM_COOLING
    !! PREPROCESS_END

  end function cooling_custom

end module prizmo_cooling_custom
