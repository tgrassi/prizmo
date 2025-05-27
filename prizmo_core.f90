module prizmo_core
  use prizmo_commons
  use prizmo_rates
  use prizmo_ode
  use prizmo_rates_heating
  use prizmo_rates_photo
  use prizmo_tdust
  use prizmo_heating_photoelectric
  use prizmo_utils
contains

  ! *************************
  subroutine init(x, Tgas, jflux)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto)

    rho_gas = get_rho(x)
    rho_dust = rho_gas * d2g

    call compute_photorates(x, Tgas,  jflux)
    call compute_photorates_heating(x, Tgas, jflux)
    call compute_Eabsorption(jflux)
    call compute_photoelectric_terms(jflux)

  end subroutine init

  ! *************************
  subroutine evolve(x, Tgas, jflux, dt, verboseChem, errState)
    implicit none
    real*8,intent(inout)::x(nspecies), Tgas
    real*8,intent(in)::jflux(nphoto), dt
    real*8::y(nspecies+1), Tdust
    integer,parameter::neq=nspecies+1
    external::dlsodes, xsetf
    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer,parameter::lrw=20+3*neq**2+20*neq
    integer,parameter::liw=30!+nmols+nmols**2
    real*8::rtol(neq)
    real*8::atol(neq)
    integer::neqa(1),itol,itask,istate,iopt,mf
    integer::iwork(liw), i, j, unit
    real*8::rwork(lrw), tloc, ytot
    integer,parameter::maxsteps=200, checksteps=10
    real*8::x_steps(nspecies, maxsteps), tgas_steps(maxsteps), dy_steps(nspecies+1, maxsteps), delWeight_steps(nspecies+1, maxsteps)
    real*8::cool_steps(ncooling, maxsteps), heat_steps(nheating, maxsteps), flux(nreactions)
    integer,intent(in)::verboseChem
    logical::MAXSTEPstop=.True.
    integer,intent(out):: errState

    !if(verboseChem==1) print*, "Begin evolve"

    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are arrays

    !number of equations is a single sized array
    neqa(:) = neq

    !  FROM DLSODES manual
    !  Name    Location   Meaning and default value
    !  ------  ---------  -----------------------------------------------
    !  H0      RWORK(5)   Step size to be attempted on the first step.
    !                     The default value is determined by the solver.
    !  HMAX    RWORK(6)   Maximum absolute step size allowed.  The
    !                     default value is infinite.
    !  HMIN    RWORK(7)   Minimum absolute step size allowed.  The
    !                     default value is 0.  (This lower bound is not
    !                     enforced on the final step before reaching
    !                     TCRIT when ITASK = 4 or 5.)
    !  MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
    !                     is 12 if METH = 1, and 5 if METH = 2. (See the
    !                     MF description above for METH.)  If MAXORD
    !                     exceeds the default value, it will be reduced
    !                     to the default value.  If MAXORD is changed
    !                     during the problem, it may cause the current
    !                     order to be reduced.
    !  MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
    !                     allowed during one call to the solver.  The
    !                     default value is 500.
    !  MXHNIL  IWORK(7)   Maximum number of messages printed (per
    !                     problem) warning that T + H = T on a step
    !                     (H = step size).  This must be positive to
    !                     result in a nondefault value.  The default
    !                     value is 10.


    itask = 1
    iopt = 1
    iwork(6) = int(1e4)  ! maximum number of iteration before warning
    iwork(5) = 0  ! maximum integration order, zero means default
    call XSETF(0)  ! 0 = print no output
    ! 222: JAC and IA/JA both generated
    ! 121: IA/JA generated, JAC provided
    !  21: JAC and IA/JA both provided
    MF = 222
    !rwork(1) = dt

    ! tolerances
    atol = ode_atol
    rtol = ode_rtol

    !if(verboseChem==1) print*, "Inititated evolve"

    y(1:nspecies) = x(:)
    y(idx_Tgas) = Tgas
    y(idx_E) = get_electrons(y)
    tloc = 0d0
    neqa = nspecies+1

    jflux_common = jflux

    Tdust = get_tdust(log10(tgas), log10(sum(x)))

    if(verboseChem==1) print*, "Got inputs to evolve:"
    if(verboseChem==1) print*, "Tgas:", Tgas
    !if(verboseChem==1) print*, "sum(jflux):", sum(jflux, dim=1)
    if(verboseChem==1) print*, "tracer densities:"
    if(verboseChem==1) print*, x
    if(verboseChem==1) print*, "Advance for dt =", dt
    errState = 0

    ytot = sum(max(y(1:nspecies), 0e0))
    if(minval(y(1:nspecies))  < -ytot*1d-10) then
      print *, "ERROR BEFORE: negative species!"
      print *, "Minloc/nspecies", minloc(y), "/", nspecies
      do i=1,nspecies
        print *, i, y(i)
      end do
      stop
    end if
    if(y(idx_Tgas)  <= 0.d0) then
      print *, "WARNING BEFORE: non-positive temperature:", y(idx_Tgas)
    end if

    istate = 1
    !call solver, see DLSODES documentation
    do i=1,maxsteps
      CALL DLSODES(fex, neqa(:), y(:), tloc, dt, &
            ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
            LIW, JES, MF)

      ! istate 2 means solution found, so no need for more steps
      if(istate == 2) exit

      ! store data at each large-step to check solution if maxstep is reached
      tgas = y(idx_tgas)
      flux(:) = get_flux(y(1:nspecies), Tgas, Tdust)
      x_steps(:, i) = y(1:nspecies)
      tgas_steps(i) = y(idx_tgas)
      heat_steps(:, i) = heating_array(y(1:nspecies), Tgas, Tdust)
      cool_steps(:, i) = cooling_array(y(1:nspecies), Tgas, Tdust, jflux, flux)
      dy_steps(:, i) = get_fex(y(1:nspecies), Tgas)

      ! Relative contributions to the initial delta
      delWeight_steps(:, i) = dy_steps(:, i) / ( y(:) * rtol + atol )

      ! reached max number of large-steps
      if(mod(i, maxsteps) == checksteps) then

        ! check if solution is oscillating around cool-heating = 0
        if(minval(sum(heat_steps(:, 1:i), dim=1) - sum(cool_steps(:, 1:i), dim=1)) &
           * maxval(sum(heat_steps(:, 1:i), dim=1) - sum(cool_steps(:, 1:i), dim=1)) < 0d0) then
          print *, "WARNING: MAXSTEPS with oscillating solution!"
          exit
        else
          if(i == maxsteps) then
              if(MAXSTEPstop) then
                  print *, "ERROR: MAXSTEPS but no oscillating solution!", istate, dt
                  stop
              else
                  print *, "WARNING: MAXSTEPS but no oscillating solution!", istate, dt
              ! save iterations to file
              open(newunit=unit, file="maxsteps.out", status="replace")
              do j=1,maxsteps
                write(unit, "(99E17.8e3)") tgas_steps(j), cool_steps(:, j), heat_steps(:, j), x_steps(:, j), dy_steps(:, j)
              end do
              close(unit)
              open(newunit=unit, file="delWeight.out", status="replace")
              do j=1,maxsteps
                write(unit, "(99E17.8e3)") delWeight_steps(:, j)
              end do
              close(unit)
              endif
          else
            print *, "WARNING: CHECKSTEPS but no oscillating solution! continue..."
          end if
        end if

      end if

      istate = 3
    end do

    ytot = sum(max(y(1:nspecies), 0e0))
    if(minval(y(1:nspecies))  < -ytot*1d-10) then
      print *, "ERROR AFTER: negative species!"
      print *, "Minloc/nspecies", minloc(y), "/", nspecies
      do i=1,nspecies
        print *, i, y(i)
      end do
      stop
    end if
    if(y(idx_Tgas)  <= 0.d0) then
      print *, "WARNING AFTER: non-positive temperature:", y(idx_Tgas)
      errState = 1
    end if

    y = max(y, 0d0)

    x = y(1:nspecies)
    x(idx_e) = get_electrons(x)
    Tgas = y(idx_Tgas)

  end subroutine evolve

end module prizmo_core
