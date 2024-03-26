module prizmo_cooling_H2
  use prizmo_commons
contains

  ! ***************
  function cooling_H2(x, log_tgas_in) result(cool)
    implicit none
    real*8,intent(in)::x(nspecies), log_tgas_in
    real*8::cool, LDL, HDL, coola(cool_H2_vec), log_tgas

    log_tgas = min(log_tgas_in, cooling_H2data%xmax*0.9999999)

    coola(:) = 1d1**interp_1dfit_vec(log_tgas, cooling_H2data, cool_H2_vec, cool_H2_nx)

    cool = 0d0
    cool = cool + coola(1) * x(idx_H)
    cool = cool + coola(2) * x(idx_Hj)
    cool = cool + coola(3) * x(idx_H2)
    cool = cool + coola(4) * x(idx_e)

    ! this to avoid negative, overflow and useless calculations below
    if(cool <= 0d0) then
       cool = 0d0
       return
    end if

    HDL = coola(5) !1d1**interp_1d(log_tgas, cooling_H2_HDL, cool_H2_min, cool_H2_fact, cool_H2_invdx)

    if (HDL == 0d0) then
       cool = 0d0
    else
      LDL = cool !erg/s
      cool = max(x(idx_H2), 0d0) / (1d0 / HDL + 1d0 / LDL) !erg/cm3/s
    endif

  end function cooling_H2

end module prizmo_cooling_H2
