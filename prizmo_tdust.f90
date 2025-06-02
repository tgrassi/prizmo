module prizmo_tdust
  use prizmo_commons
  use prizmo_fit
contains

  ! ***************
  function get_tdust(log_Tgas_in, log_ngas_in) result(tdust)
    implicit none
    real*8,intent(in)::log_Tgas_in, log_ngas_in
    real*8::tdust, f(nphoto), log_ngas, log_tgas

    log_Tgas = max(min(log_Tgas_in, tdust_table_data%ymax*0.99999), tdust_table_data%ymin)
    log_ngas = max(min(log_ngas_in, tdust_table_data%zmax*0.99999), tdust_table_data%zmin)

    tdust = 1d1**interp_3dfit(log_Eabsorption, log_Tgas, log_ngas, tdust_table_data, &
      dust_cooling_table_n1, dust_cooling_table_n2, dust_cooling_table_n3)

  end function get_tdust

  ! ********************
  subroutine compute_Eabsorption(jflux)
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::f(nphoto), Eabsorption

    f(:) = Jflux * pre_dust_cooling_table(:) / hplanck
    Eabsorption = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0

    log_Eabsorption = log10(Eabsorption + 1d-40)
    log_Eabsorption = min(max(log_Eabsorption, tdust_table_data%xmin), tdust_table_data%xmax*0.99999)  ! FIXME

  end subroutine compute_Eabsorption

end module prizmo_tdust
