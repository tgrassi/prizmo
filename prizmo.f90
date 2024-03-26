module prizmo
  use prizmo_commons
  use prizmo_core
  use prizmo_loaders
  use prizmo_attenuate
  integer,parameter::prizmo_nspecies=nspecies
  integer,parameter::prizmo_nphoto=nphoto
  integer,parameter::prizmo_nreactions=nreactions

  !! PREPROCESS_INDEXES
  !! PREPROCESS_END

contains

  ! ************************
  subroutine prizmo_init()
    implicit none

    runtime_data_folder = "runtime_data/"

    call load_energy()
    call load_all_photo_xsecs()
    call load_energy_thresholds()
    call load_all_atomic_cooling_tables()
    call load_dust_cooling_table()
    !call load_photoelectric_tables() ! NOTE: not used
    call load_dust_kappa_opacity()
    !! PREPROCESS_H2
    !! PREPROCESS_END
    !! PREPROCESS_CO
    !! PREPROCESS_END
    call load_verbatim_reactions()

    print *, "everything loaded!"

    gamma_ad = 7./5.
    d2g = 1d-2
    ortho_to_para = 3. / 1.
    user_Av = 0d0
    radial_Ncol_H2 = 0d0  ! 1/cm2
    radial_Ncol_CO = 0d0  ! 1/cm2
    vertical_Ncol_CO = 0d0  ! 1/cm2
    user_cr = 5d-17  ! 1/s
    solve_thermo = .true.
    solve_chemistry = .true.

    ode_atol = 1d-20  ! default absolute tolerance
    ode_rtol = 1d-8  ! default relative tolerance

  end subroutine prizmo_init

  ! ************************
  subroutine prizmo_evolve(x, Tgas, jflux, dt, verboseChem, errState)
    implicit none
    real*8,intent(inout)::x(nspecies), Tgas, jflux(nphoto)
    real*8,intent(in)::dt
    integer,intent(in)::verboseChem
    integer,intent(out)::errState

    call init(x, Tgas, jflux)
    call evolve(x, Tgas, jflux, dt, verboseChem, errState)

  end subroutine prizmo_evolve

  ! ************************
  subroutine prizmo_evolve_rho(xx, rho, Tgas, jflux, dt, verboseChem, errState)
    use prizmo_commons
    implicit none
    real*8,intent(inout)::xx(nspecies), Tgas, jflux(nphoto)
    real*8,intent(in)::dt, rho
    real*8::x(nspecies)
    integer,intent(in)::verboseChem
    integer,intent(out)::errState

    x = rho * xx / masses

    call prizmo_evolve(x, Tgas, jflux, dt, verboseChem, errState)

    xx = x * masses / rho

  end subroutine prizmo_evolve_rho

  ! **********************
  subroutine prizmo_rt(x, Tgas, jflux, ds)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, ds
    real*8,intent(inout)::jflux(nphoto)

    call attenuate(x, Tgas, jflux, ds)

  end subroutine prizmo_rt

  ! **********************
  subroutine prizmo_rt_rho(xx, rho, Tgas, jflux, ds)
    implicit none
    real*8,intent(in)::xx(nspecies), rho, Tgas, ds
    real*8,intent(inout)::jflux(nphoto)
    real*8::x(nspecies)

    x = rho * xx / masses

    call attenuate(x, Tgas, jflux, ds)

  end subroutine prizmo_rt_rho

  ! **************
  function prizmo_is_thermal_equilibrium(x, Tgas, jflux, dt) result(is_eq)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto), dt
    real*8::delta_tgas, cools(ncooling), heats(nheating), Tdust, rel_error
    logical::is_eq

    rel_error = 1d-3

    Tdust = prizmo_get_tdust(x, Tgas, jflux)
    cools = prizmo_get_cooling_array(x, Tgas, Tdust, jflux)
    heats = prizmo_get_heating_array(x, Tgas, Tdust, jflux)

    delta_tgas = (gamma_ad - 1d0) * abs(sum(heats) - sum(cools)) / kboltzmann / sum(x) * dt
    is_eq = delta_tgas / Tgas < rel_error

  end function prizmo_is_thermal_equilibrium

  ! **************
  function prizmo_get_rho(x) result(rho)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::rho

    rho = sum(masses * x)

  end function prizmo_get_rho

  ! **********************
  function prizmo_n2frac(x) result(xx)
    use prizmo_commons
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::xx(nspecies)

    xx = x * masses / sum(x * masses)

  end function prizmo_n2frac

  ! **********************
  function prizmo_frac2n(xx, rho) result(x)
    use prizmo_commons
    implicit none
    real*8,intent(in)::xx(nspecies), rho
    real*8::x(nspecies)

    x = xx * rho / masses

  end function prizmo_frac2n

  ! **************
  subroutine prizmo_set_atol_all(val)
    implicit none
    real*8,intent(in)::val

    ode_atol = val

  end subroutine prizmo_set_atol_all

  ! **************
  subroutine prizmo_set_atol_idx(val, idx)
    implicit none
    real*8,intent(in)::val
    integer,intent(in)::idx

    ode_atol(idx) = val

  end subroutine prizmo_set_atol_idx

  ! **************
  function prizmo_load_radiation_field(filename) result(field)
    implicit none
    character(len=*),intent(in)::filename
    real*8::field(nphoto)
    integer::i, unit

    open(newunit=unit, file=trim(filename), status="old")
    do i=1,nphoto
      read(unit, *) field(i)
    end do
    close(unit)

  end function prizmo_load_radiation_field

  ! ********************
  function prizmo_get_atomic_cooling_array(x, Tgas) result(coola)
    use prizmo_cooling_atomic
    implicit none
    real*8,intent(in)::x(nspecies), Tgas
    real*8::coola(natomic_cools)

    coola = cooling_atomic_array(x, log10(Tgas))

  end function prizmo_get_atomic_cooling_array

  ! ******************
  function prizmo_get_cooling_array(x, Tgas, Tdust, jflux) result(cools)
    use prizmo_cooling
    use prizmo_flux
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::fluxes(nreactions), cools(ncooling)

    call init(x, Tgas, jflux)

    fluxes(:) = get_flux(x, Tgas, Tdust)
    cools(:) = cooling_array(x, Tgas, Tdust, jflux, fluxes)

  end function prizmo_get_cooling_array

  ! ******************
  function prizmo_get_heating_array(x, Tgas, Tdust, jflux) result(heats)
    use prizmo_heating
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::heats(nheating)

    call init(x, Tgas, jflux)
    heats(:) = heating_array(x, Tgas, Tdust, jflux)

  end function prizmo_get_heating_array

  ! ***************************
  subroutine prizmo_save_cooling_function(x, jflux, fname)
    use prizmo_cooling
    use prizmo_flux
    implicit none
    integer,parameter::ntemp=100
    real*8,intent(in)::x(nspecies), jflux(nphoto)
    character(len=*),intent(in)::fname
    real*8::fluxes(nreactions), cools(ncooling), tmin, tmax, Tgas, Tdust
    integer::unit, i

    tmin = log10(1d0)
    tmax = log10(1d5)

    open(newunit=unit, file=trim(fname), status="replace")
      do i=1, ntemp
        tgas = 1d1**((i - 1) * (tmax - tmin) / (ntemp - 1) + tmin)
        tdust = tgas
        fluxes(:) = get_flux(x, Tgas, Tdust)
        cools(:) = cooling_array(x, Tgas, Tdust, jflux, fluxes)
        write(unit, '(99e17.8e3)') tgas, cools
      end do
    close(unit)

  end subroutine prizmo_save_cooling_function

  ! ***************************
  subroutine prizmo_save_heating_function(x, jflux, fname)
    use prizmo_heating
    implicit none
    integer,parameter::ntemp=100
    real*8,intent(in)::x(nspecies), jflux(nphoto)
    character(len=*),intent(in)::fname
    real*8::heats(nheating), tmin, tmax, Tgas, Tdust
    integer::unit, i

    tmin = log10(1d0)
    tmax = log10(1d5)

    open(newunit=unit, file=trim(fname), status="replace")
      do i=1, ntemp
        tgas = 1d1**((i - 1) * (tmax - tmin) / (ntemp - 1) + tmin)
        tdust = tgas
        heats(:) = heating_array(x, Tgas, Tdust, jflux)
        write(unit, '(99e17.8e3)') tgas, heats
      end do
    close(unit)

  end subroutine prizmo_save_heating_function

  ! ****************************
  function prizmo_get_energy()
    implicit none
    real*8::prizmo_get_energy(nphoto)

    prizmo_get_energy = energy

  end function prizmo_get_energy

  ! ****************************
  function prizmo_get_energy_ev()
    implicit none
    real*8::prizmo_get_energy_ev(nphoto)

    prizmo_get_energy_ev = energy / ev2erg

  end function prizmo_get_energy_ev

  ! ****************************
  subroutine prizmo_set_user_Av(Av)
    implicit none
    real*8,intent(in)::Av

    user_Av = Av

  end subroutine prizmo_set_user_Av

  ! ****************************
  subroutine prizmo_set_d2g(val)
    implicit none
    real*8,intent(in)::val

    d2g = val

  end subroutine prizmo_set_d2g

  ! ****************************
  subroutine prizmo_set_crate(val)
    implicit none
    real*8,intent(in)::val

    user_cr = val

  end subroutine prizmo_set_crate

  ! ****************************
  subroutine prizmo_set_radial_Ncol_H2(val)
    implicit none
    real*8,intent(in)::val

    radial_Ncol_H2 = val

  end subroutine prizmo_set_radial_Ncol_H2

  ! ****************************
  subroutine prizmo_set_radial_Ncol_CO(val)
    implicit none
    real*8,intent(in)::val

    radial_Ncol_CO = val

  end subroutine prizmo_set_radial_Ncol_CO

  ! ****************************
  subroutine prizmo_set_vertical_Ncol_CO(val)
    implicit none
    real*8,intent(in)::val

    vertical_Ncol_CO = val

  end subroutine prizmo_set_vertical_Ncol_CO

  ! ********************
  function prizmo_get_rates(x, Tgas, Tdust, jflux) result(k)
    use prizmo_rates
    use prizmo_rates_photo
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::k(nreactions)

    call compute_rates(x, Tgas, Tdust)
    call compute_photorates(x, Tgas, jflux)

    k = kall

  end function prizmo_get_rates

  ! ***************
  function prizmo_get_tdust(x, Tgas, jflux) result(tdust)
    implicit none
    real*8,intent(in)::jflux(nphoto), Tgas, x(nspecies)
    real*8::tdust, log_ngas, log_tgas

    log_ngas = log10(sum(x))
    log_tgas = log10(Tgas)

    call compute_Eabsorption(jflux)

    tdust = get_tdust(log_Tgas, log_ngas)

  end function prizmo_get_tdust

  ! ******************
  subroutine prizmo_set_solve_thermo(val)
    implicit none
    logical,intent(in)::val

    solve_thermo = val

  end subroutine prizmo_set_solve_thermo

  ! ******************
  subroutine prizmo_set_solve_chemistry(val)
    implicit none
    logical,intent(in)::val

    solve_chemistry = val

  end subroutine prizmo_set_solve_chemistry

  ! ***********************
  subroutine print_ranked_fluxes(x, Tgas, jflux)
    use prizmo_rates
    use prizmo_rates_photo
    use prizmo_flux
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto)
    real*8::fluxes(nreactions), Tdust

    Tdust = prizmo_get_tdust(x, Tgas, jflux)

    call compute_rates(x, Tgas, Tdust)
    call compute_photorates(x, Tgas, jflux)

    fluxes = get_flux(x, Tgas, Tdust)

    call ranker(fluxes, 60, reactions_verbatim)

  end subroutine print_ranked_fluxes

  ! ***************
  function prizmo_get_Xnuclei(x, atom) result(nX)
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies)
    character(len=*),intent(in)::atom
    real*8::nX

    if(trim(atom) == "H") then
      nX = get_Hnuclei(x)
    else if(trim(atom) == "C") then
      nX = get_Cnuclei(x)
    else if(trim(atom) == "O") then
      nX = get_Onuclei(x)
    else if(trim(atom) == "He") then
      nX = get_HEnuclei(x)
    else
      print *, "ERROR: unknown atom "//trim(atom)//" in prizmo_get_Xnuclei"
      stop
    end if

  end function prizmo_get_Xnuclei

  ! ************************
  ! ratio between the integrated flux in the FUV range, and the Habing flux in the same range
  function prizmo_get_chi_FUV(jflux) result(xfuv)
    use prizmo_heating_photoelectric
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::xfuv

    call compute_photoelectric_terms(jflux)

    xfuv = chi_FUV

  end function prizmo_get_chi_FUV

end module prizmo
