module prizmo_heating_H2diss
  use prizmo_commons
contains

  ! ***************
  function heating_H2diss(x, Tgas, ntot) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, ntot
    real*8::ncrn, ncrd1, ncrd2, yH, yH2, ncr, hf, rdiss
    real*8::heat

    !! PREPROCESS_H2DISS
    !! PREPROCESS_END

    ! H2 photoheating
    ncrn = 1d6 / sqrt(Tgas)
    ncrd1 = 1.6 * exp(-(4d2 / Tgas)**2) + 1d-40
    ncrd2 = 1.4 * exp(-1.2d4 / (Tgas + 1.2d3)) + 1d-40

    yH = x(idx_H) / ntot + 1d-40
    yH2 = x(idx_H2) / ntot + 1d-40

    ncr = ncrn / (ncrd1 * yH + ncrd2 * yH2)
    hf = ntot / (ntot + ncr)

    ! H2 dissociation heating, erg/s/cm3
    heat = (6.4d-13 + 2.7d-11 * hf) * Rdiss * x(idx_H2)

  end function heating_H2diss

end module prizmo_heating_H2diss
