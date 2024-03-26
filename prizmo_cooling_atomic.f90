module prizmo_cooling_atomic
  use prizmo_commons
  use prizmo_fit
contains

  ! ***************
  function cooling_atomic(x, log_Tgas_in, Tgas) result(cool)
    implicit none
    real*8,intent(in)::x(nspecies), log_Tgas_in, Tgas
    real*8::cool !, coola2d(atomic_cooling_nvec2d), coola3d(atomic_cooling_nvec3d), coola4d(atomic_cooling_nvec4d)
    real*8::xlimit, tmin, tmax, log_tgas
    !real*8::log_xnr, log_xpr, log_xer, log_tgas
    !real*8::nmin, nmax, tmin, tmax, xlimit

    xlimit = -1d99

    tmin = 0d0
    tmax = 6d0 * 0.99999

    log_tgas = min(max(log_Tgas_in, tmin), tmax)

    cool = 0d0

    ! Ly-alpha H cooling
    cool = cool + atomic_cooling_H(x, Tgas)

    !! PREPROCESS_ATOMIC_COOLING
    !! PREPROCESS_END

  end function cooling_atomic

  ! **********************
  function cooling_atomic_array(x, log_Tgas) result(cools)
    use prizmo_commons
    implicit none
    real*8,intent(in)::x(nspecies), log_Tgas
    real*8::cools(natomic_cools)

    !! PREPROCESS_ATOMIC_COOLING_ARRAY
    !! PREPROCESS_END

  end function cooling_atomic_array

  ! ***************************
  function atomic_cooling_H(x, tgas) result(cool)
    use prizmo_commons
    implicit none
    real*8,intent(in)::x(nspecies), tgas
    real*8::cool

    cool = 7.39d-19 * x(idx_H) * x(idx_E) * exp(-118400.d0 / tgas)

  end function atomic_cooling_H

  !! PREPROCESS_ATOMIC_COOLING_FUNCTIONS
  !! PREPROCESS_END

  ! ! **********************
  ! function cooling_atomic_array2d(log_xer, log_tgas) result(coola)
  !   implicit none
  !   real*8,intent(in)::log_xer, log_tgas
  !   real*8::coola(atomic_cooling_nvec2d)
  !
  !   coola(:) = 1d1**interp_2dfit_vec(log_xer, log_Tgas, atomic_cooling_table_2d, &
  !         atomic_cooling_nvec2d, atomic_cooling2d_n1, atomic_cooling2d_n2)
  !
  ! end function cooling_atomic_array2d
  !
  ! ! **********************
  ! function cooling_atomic_array3d(log_xnr, log_xer, log_tgas) result(coola)
  !   implicit none
  !   real*8,intent(in)::log_xnr, log_xer, log_tgas
  !   real*8::coola(atomic_cooling_nvec3d)
  !
  !   coola(:) = 1d1**interp_3dfit_vec(log_xnr, log_xer, log_Tgas, atomic_cooling_table_3d, &
  !         atomic_cooling_nvec3d, atomic_cooling3d_n1, atomic_cooling3d_n2, atomic_cooling3d_n3)
  !
  ! end function cooling_atomic_array3d
  !
  ! ! **********************
  ! function cooling_atomic_array4d(log_xpr, log_xnr, log_xer, log_tgas) result(coola)
  !   implicit none
  !   real*8,intent(in)::log_xpr, log_xnr, log_xer, log_tgas
  !   real*8::coola(atomic_cooling_nvec4d)
  !
  !   coola(:) = 1d1**interp_4dfit_vec(log_xnr, log_xpr, log_xer, log_Tgas, atomic_cooling_table_4d, &
  !         atomic_cooling_nvec4d, atomic_cooling4d_n1, atomic_cooling4d_n2, atomic_cooling4d_n3, atomic_cooling4d_n4)
  !
  ! end function cooling_atomic_array4d

end module prizmo_cooling_atomic
