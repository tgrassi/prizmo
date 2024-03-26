module prizmo_cooling_dust
  use prizmo_commons
  use prizmo_fit
contains

  ! ***************
  function cooling_dust(log_Tgas_in, log_ngas_in) result(cool)
    implicit none
    real*8,intent(in)::log_Tgas_in, log_ngas_in
    real*8::cool, Eabsorption, f(nphoto), log_ngas, log_tgas

    if(d2g < d2g_min) then
      cool = 0d0
      return
    end if

    log_tgas = max(min(log_tgas_in, tdust_table_data%ymax*0.99999), tdust_table_data%ymin)  ! FIXME
    log_ngas = max(min(log_ngas_in, tdust_table_data%zmax*0.99999), tdust_table_data%zmin)  ! FIXME

    ! cooling
    cool = 1d1**interp_3dfit(log_Eabsorption, log_Tgas, log_ngas, dust_cooling_table_data, &
      dust_cooling_table_n1, dust_cooling_table_n2, dust_cooling_table_n3) * rho_dust

    ! heating
    cool = cool - 1d1**interp_3dfit(log_Eabsorption, log_Tgas, log_ngas, dust_heating_table_data, &
      dust_cooling_table_n1, dust_cooling_table_n2, dust_cooling_table_n3) * rho_dust

  end function cooling_dust

end module prizmo_cooling_dust
