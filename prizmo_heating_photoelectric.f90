module prizmo_heating_photoelectric
  use prizmo_commons
contains

  ! ***************
  ! NOTE: this is not used
  function heating_photoelectric_old(x, tgas, jflux) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies), tgas, jflux(nphoto)
    real*8::heat, jion(zmin:zmax), jele(zmin:zmax), y(zmin:zmax), f(nphoto)
    real*8::jion_cool(zmin:zmax), jele_cool(zmin:zmax), xions, log_tgas
    integer::i

    if(d2g < d2g_min) then
      heat = 0d0
      return
    end if

    log_tgas = log10(tgas)

    xions = x(idx_Hj)  ! FIXME x(idx_Hj) + x(idx_Hj) + x(idx_Cj) + x(idx_Oj) +  x(idx_Cj)

    jion(:) = jie_fit(log_tgas, jtab_fit_nt, jion_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx) * (xions + 1d-40)
    jele(:) = jie_fit(log_tgas, jtab_fit_nt, jele_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx) * (x(idx_E) + 1d-40)

    y(0) = 1d0
    do i=1,zmax
      y(i) = min(y(i-1) * (phterm_jpe(i-1) + jion(i-1) + 1d-40) / (jele(i) + 1d-40), 1d80)
    end do

    do i=-1,zmin,-1
      y(i) = y(i+1) / (phterm_jpe(i) + jion(i) + 1d-40) * (jele(i+1) + 1d-40)
    end do

    y = y / sum(y)

    !print '(99e17.8e3)', y
    jion_cool(:) = jie_fit(log_tgas, jtab_fit_nt, jion_cool_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx) * (xions + 1d-40)
    jele_cool(:) = jie_fit(log_tgas, jtab_fit_nt, jele_cool_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx) * (x(idx_E) + 1d-40)

    ! print *, sum(y * jpe_heating), sum(y * jion_cool), sum(y * jele_cool)

    heat = rho_dust * (sum(y * phterm_jpe_heating) - sum(y * jion_cool) - sum(y * jele_cool))

  end function heating_photoelectric_old


  ! ************************
  function jie_fit(log_tgas, nt, jfit_data, jfit_xmin, jfit_dx, jfit_invdx) result(f)
    implicit none
    real*8,intent(in)::log_tgas, jfit_data(zmin:zmax, nt), jfit_xmin, jfit_dx, jfit_invdx
    integer,intent(in)::nt
    real*8::f(zmin:zmax), x0
    integer::idx

    idx = floor((log_tgas - jfit_xmin) * jfit_invdx) + 1
    x0 = (idx - 1) * jfit_dx + jfit_xmin
    f(:) = (log_tgas - x0) * jfit_invdx * (jfit_data(:, idx+1) - jfit_data(:, idx)) + jfit_data(:, idx)

    f = 1d1**f

  end function jie_fit

  ! *******************
  subroutine compute_photoelectric_terms(jflux)
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::f(nphoto)
    integer::i, idx1, idx2

    idx1 = fuv_idx1
    idx2 = fuv_idx2

    ! FIXME
    ! do i=zmin, zmax
    !   f(:) = Jflux * jpe_table(:, i)
    !   phterm_jpe(i) = sum((f(idx1+1:idx2) + f(idx1:idx2-1)) * delta_energy(idx1:idx2-1)) / 2d0
    ! end do
    !
    ! do i=zmin, zmax
    !   f(:) = Jflux * jpe_heating_table(:, i)
    !   phterm_jpe_heating(i) = sum((f(idx1+1:idx2) + f(idx1:idx2-1)) * delta_energy(idx1:idx2-1)) / 2d0
    ! end do

    ! ratio between the integrated flux in the FUV range, and the Habing flux in the same range
    f = jflux / energy
    chi_FUV = sum((f(idx1+1:idx2) + f(idx1:idx2-1)) * delta_energy(idx1:idx2-1)) / 2. / habing_flux

  end subroutine compute_photoelectric_terms

  ! ****************************
  function heating_photoelectric(x, tgas_in) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies), tgas_in
    real*8::heat
    real*8::xx, y, eps, sigma, xe, tgas, xfuv

    xe = min(max(x(idx_E), 1d-5), 1d5)
    tgas = min(max(tgas_in, 1d1), 1d4)
    xfuv = min(max(chi_FUV, 1d-5), 1d5)

    xx = sqrt(tgas) * xfuv / (xe + 1d-40)

    if(xx <= 1d-4) then
      y = 0.7
    elseif(xx > 1d-4 .and. xx <= 1d0) then
      y = 0.0171875*xx**3 + 0.103125*xx**2 + 0.15 !0.36
    else
      y = 0.15
    endif

    eps = 0.06 / (1d0 + 1.8d-3 * xx**.91) + y * (1d-4 * tgas)**1.2 / (1d0 + 1d-2 * xx)

    sigma = kabs_integral * rho_dust

    heat = 2.5d-4 * sigma * eps * xfuv

  end function heating_photoelectric

end module prizmo_heating_photoelectric
