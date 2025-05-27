module prizmo_cooling
  use prizmo_commons
  use prizmo_cooling_atomic
  use prizmo_cooling_dust
  use prizmo_cooling_chemical
  !! PREPROCESS_USE_MOLECULES
  !! PREPROCESS_END
  integer,parameter::ncooling=5
contains

  ! **********************
  function cooling(xin, Tgas, Tdust, fluxin) result(cool)
    implicit none
    real*8,intent(in)::xin(nspecies), Tgas, Tdust, fluxin(nreactions)
    real*8::cool, cools(ncooling), x(nspecies), fluxes(nreactions)

    fluxes = fluxin
    x = max(xin, 0d0)
    !fluxes = max(fluxin, 0d0)

    cools = cooling_array(x, Tgas, Tdust, fluxes)

    cool = sum(cools)

  end function cooling

  ! **********************
  function cooling_array(x, Tgas, Tdust, fluxes) result(cools)
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, fluxes(nreactions)
    real*8::cools(ncooling), log_Tgas, log_ngas, Eabsorption, log_Hnuclei

    log_ngas = log10(sum(x))
    log_Tgas = log10(Tgas)
    log_Hnuclei = log10(get_Hnuclei(x))

    cools(1) = cooling_atomic(x, log_Tgas, Tgas)
    cools(2) = cooling_chemical(x, fluxes, Tgas)
    cools(3) = cooling_dust(log_Tgas, log_ngas)
    !! PREPROCESS_MOLECULAR_COOLING
    !! PREPROCESS_END

  end function cooling_array

end module prizmo_cooling
