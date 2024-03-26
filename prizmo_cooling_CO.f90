module prizmo_cooling_CO
  use prizmo_commons
contains

  ! ***************
  function cooling_CO(x, log_tgas, log_Hnuclei) result(cool)
    implicit none
    real*8,intent(in)::x(nspecies), log_tgas, log_Hnuclei
    real*8::cool, log_NCO, log_t, log_h

    log_h = min(max(log_Hnuclei, cool_CO_tab_data%ymin), cool_CO_tab_data%ymax*0.9999)
    log_t = min(max(log_tgas, cool_CO_tab_data%xmin), cool_CO_tab_data%xmax*0.9999)
    log_NCO = max(log10(vertical_Ncol_CO + 1d-40), cool_CO_tab_data%zmin)

    cool = 1d1**interp_3dfit(log_t, log_h, log_NCO, &
      cool_CO_tab_data, cool_CO_tab_n1, cool_CO_tab_n2, cool_CO_tab_n3) * x(idx_CO)

  end function cooling_CO

end module prizmo_cooling_CO
