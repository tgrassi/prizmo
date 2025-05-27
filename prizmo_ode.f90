module prizmo_ode
  use prizmo_commons
  use prizmo_flux
  use prizmo_heating
  use prizmo_cooling
  use prizmo_tdust
contains

  ! ***************************
  subroutine fex(neq, tt, yin, dy)
    integer::neq, i
    real*8::yin(neq), tt, dy(neq), y(neq)
    real*8::flux(nreactions), Tdust, ntot, log_tgas, log_ngas
    real*8:: heat, cool, Tgas

    !y = max(yin, 0d0)
    y = yin
    Tgas = max(y(idx_tgas), 3d0)  ! FIXME

    ntot = sum(max(y(1:nspecies), 1d-40))

    log_tgas = log10(Tgas)
    log_ngas = log10(ntot)
    Tdust = get_tdust(log_Tgas, log_ngas)

    flux(:) = get_flux(y(1:nspecies), Tgas, Tdust)

    if(minval(kall) < 0d0) then
      print *, "ERROR: negative kall!!!"
      print *, "temperature, K", Tgas
      do i=1,nreactions
        print *, i, kall(i)
      end do
      do i=1,nspecies
        print *, i, y(i)
      end do
      stop
    end if

    if(solve_thermo) then
      heat = heating(y(1:nspecies), Tgas, Tdust)
      cool = cooling(y(1:nspecies), Tgas, Tdust, flux)

      dy(idx_Tgas) = (gamma_ad - 1d0) * (heat - cool) / kboltzmann / ntot
    else
      dy(idx_Tgas) = 0d0
    end if

    if(solve_chemistry) then
      !! PREPROCESS_ODE
      !! PREPROCESS_END
    end if

  end subroutine fex

  ! ***************************
  ! Jacobian, pd(i,j)=df(i)/dx(j), see DLSODES documentation
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use prizmo_commons
    implicit none
    integer::neq, j, ian, jan
    real*8::tt, n(neq), pdj(neq)

  end subroutine jes

  ! *******************
  function get_fex(x, Tgas) result(dy)
    use prizmo_commons
    implicit none
    integer,parameter::neq=nspecies+1
    real*8,intent(in)::x(nspecies), Tgas
    real*8::dy(neq), tt, y(neq)

    tt = 0d0

    y(1:nspecies) = x
    y(idx_Tgas) = Tgas

    call fex(neq, tt, y, dy)

  end function get_fex

end module prizmo_ode
