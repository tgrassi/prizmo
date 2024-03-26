module prizmo_cooling_chemical
  use prizmo_commons
  use prizmo_utils
contains

  ! ***************
  function cooling_chemical(x, fluxes, Tgas) result(cool)
    implicit none
    real*8,intent(in)::x(nspecies), fluxes(nreactions), Tgas
    real*8::coola(nreactions), cool

    coola = get_cooling_chemical_array(x, fluxes, Tgas)

    ! call ranker(abs(coola), 5)  ! DEBUG

    cool = sum(coola)

  end function cooling_chemical

  ! ***********************
  function get_cooling_chemical_array(x, fluxes, Tgas) result(coola)
    implicit none
    real*8,intent(in)::x(nspecies), fluxes(nreactions), Tgas
    real*8::coola(nreactions)

    coola = 0d0

    !! PREPROCESS_RECOMBINATION
    !! PREPROCESS_END

    coola = coola * kboltzmann * Tgas

    !! PREPROCESS_CHEMICAL
    !! PREPROCESS_END

  end function get_cooling_chemical_array

end module prizmo_cooling_chemical
