module prizmo_rates
  use prizmo_commons
  use prizmo_shielding
  use prizmo_utils
contains

  ! ************************
  subroutine compute_rates(x, Tgas, Tdust)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust
    real*8::te, lnte, Tgas32, log_tgas, invTgas, sqrTgas32, invsqrTgas32, invsqrTgas
    real*8::invTgas32, Tgas14, sqrTgas, invTe, invsqrTe
    real*8::Tgas32_m06, Tgas32_m05, Tgas32_m03, Tgas32_m04, invTgas14_065
    real*8::sticking, nu_debye

    !! PREPROCESS_PROTOTYPES_DEFINE
    !! PREPROCESS_END

    ! temperature shortcuts
    te = Tgas * 8.617343d-5
    lnte = log(te)
    invTe = 1d0 / te
    invsqrTe = 1d0 / sqrt(te)
    Tgas32 = Tgas / 3d2
    Tgas14 = Tgas / 1d4
    invTgas32 = 1d0 / Tgas32
    log_tgas = log10(Tgas)
    invTgas = 1d0 / Tgas
    sqrTgas32 = sqrt(Tgas32)
    sqrTgas = sqrt(Tgas)
    invsqrTgas32 = 1d0 / sqrTgas32
    invsqrTgas = 1d0 / sqrTgas
    Tgas32_m06 = Tgas32**(-0.6)
    Tgas32_m05 = Tgas32**(-0.5)
    Tgas32_m04 = Tgas32**(-0.4)
    Tgas32_m03 = Tgas32**(-0.3)
    invTgas14_065 = 1d0 / Tgas14**0.65

    ! dust parameters
    sticking = 1d0
    nu_debye = 1d12  ! 1/s

    !! PREPROCESS_PROTOTYPES
    !! PREPROCESS_END

    !! PREPROCESS_RATES
    !! PREPROCESS_END

  end subroutine compute_rates

end module prizmo_rates
