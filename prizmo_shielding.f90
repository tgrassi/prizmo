module prizmo_shielding
  use prizmo_commons
  use prizmo_fit
contains

  ! ********************
  function shielding_H2(log_NH2, log_tgas_in) result(f)
    implicit none
    real*8,intent(in)::log_NH2, log_tgas_in
    real*8::f, log_ncol, log_tgas

    log_ncol = min(max(log_NH2, shielding_H2_data%xmin), shielding_H2_data%xmax*0.9999)
    log_tgas = min(max(log_tgas_in, shielding_H2_data%ymin), shielding_H2_data%ymax*0.9999)

    f = 1d1**interp_2dfit(log_ncol, log_tgas, &
      shielding_H2_data, shielding_H2_n1, shielding_H2_n2)

  end function shielding_H2

  ! ********************
  function shielding_CO(log_NCO, log_NH2) result(f)
    implicit none
    real*8,intent(in)::log_NCO, log_NH2
    real*8::f, log_ncolH2, log_ncolCO

    log_ncolCO = min(max(log_NCO, shielding_CO_data%xmin), shielding_CO_data%xmax*0.9999)
    log_ncolH2 = min(max(log_NH2, shielding_CO_data%ymin), shielding_CO_data%ymax*0.9999)

    f = 1d1**interp_2dfit(log_ncolCO, log_ncolH2, &
      shielding_CO_data, shielding_CO_n1, shielding_CO_n2)

  end function shielding_CO

end module prizmo_shielding
