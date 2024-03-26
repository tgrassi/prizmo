module prizmo_commons
use prizmo_fit
implicit none

!! PREPROCESS_COMMON_VARS
!! PREPROCESS_END

integer,parameter::idx_Tgas=nspecies + 1
integer,parameter::zmin=-2, zmax=2

integer,parameter::jtab_fit_nt=30
real*8::jion_fit_data(zmin:zmax, jtab_fit_nt), jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx
real*8::jele_fit_data(zmin:zmax, jtab_fit_nt)
real*8::jion_cool_fit_data(zmin:zmax, jtab_fit_nt)
real*8::jele_cool_fit_data(zmin:zmax, jtab_fit_nt)

!! PREPROCESS_ATOMIC_COOLING_NVEC
!! PREPROCESS_END

integer,parameter::atomic_cooling_n1=10000  ! tgas, K
integer,parameter::atomic_cooling_5lev_nvec=20  ! 5 levels
integer,parameter::atomic_cooling_3lev_nvec=6  ! 3 levels
integer,parameter::atomic_cooling_2lev_nvec=2  ! 2 levels

!! PREPROCESS_ATOMIC_COOLING_COMMONS
!! PREPROCESS_END

! ! atomic cooling table, 2D
! integer,parameter::atomic_cooling2d_n1=30  ! e- abundance, cm-3
! integer,parameter::atomic_cooling2d_n2=100  ! tgas, K
! type(fit2d_data_vec(nv=atomic_cooling_nvec2d, n1=atomic_cooling2d_n1, n2=atomic_cooling2d_n2))::atomic_cooling_table_2d
!
! ! atomic cooling table, 3D
! integer,parameter::atomic_cooling3d_n1=30  ! H abundance, cm-3
! integer,parameter::atomic_cooling3d_n2=30  ! e- abundance, cm-3
! integer,parameter::atomic_cooling3d_n3=100  ! tgas, K
! type(fit3d_data_vec(nv=atomic_cooling_nvec3d, n1=atomic_cooling3d_n1, n2=atomic_cooling3d_n2, &
! n3=atomic_cooling3d_n3))::atomic_cooling_table_3d
!
! ! atomic cooling table, 4D
! integer,parameter::atomic_cooling4d_n1=30  ! H abundance, cm-3
! integer,parameter::atomic_cooling4d_n2=30  ! H+ abundance, cm-3
! integer,parameter::atomic_cooling4d_n3=30  ! e- abundance, cm-3
! integer,parameter::atomic_cooling4d_n4=100  ! tgas, K
! type(fit4d_data_vec(nv=atomic_cooling_nvec4d, n1=atomic_cooling4d_n1, n2=atomic_cooling4d_n2, &
!                     n3=atomic_cooling4d_n3, n4=atomic_cooling4d_n4))::atomic_cooling_table_4d

! dust cooling table, 3D
integer,parameter::dust_cooling_table_n1=100
integer,parameter::dust_cooling_table_n2=100
integer,parameter::dust_cooling_table_n3=50
type(fit3d_data(n1=dust_cooling_table_n1, n2=dust_cooling_table_n2, n3=dust_cooling_table_n3))::dust_cooling_table_data
type(fit3d_data(n1=dust_cooling_table_n1, n2=dust_cooling_table_n2, n3=dust_cooling_table_n3))::dust_heating_table_data
type(fit3d_data(n1=dust_cooling_table_n1, n2=dust_cooling_table_n2, n3=dust_cooling_table_n3))::tdust_table_data

! dust cooling pre-jflux
real*8::pre_dust_cooling_table(nphoto)

! photoelectic heating Jpe table
real*8::jpe_table(nphoto, zmin:zmax)
real*8::jpe_heating_table(nphoto, zmin:zmax)
real*8::phterm_jpe(zmin:zmax), phterm_jpe_heating(zmin:zmax)

! H2 cooling tables
integer,parameter::cool_H2_vec=5
integer,parameter::cool_H2_nx=1000
type(fit1d_data_vec(nv=cool_H2_vec, n1=cool_H2_nx))::cooling_H2data

! dust kabs table, cm2/g
real*8::dust_kappa_opacity(nphoto)

! shielding H2 tables
integer,parameter::shielding_H2_n1=30, shielding_H2_n2=100
type(fit2d_data(n1=shielding_H2_n1, n2=shielding_H2_n2))::shielding_H2_data

! shielding CO tables
integer,parameter::shielding_CO_n1=50, shielding_CO_n2=50
type(fit2d_data(n1=shielding_CO_n1, n2=shielding_CO_n2))::shielding_CO_data

! cooling CO tables
integer,parameter::cool_CO_tab_n1=40, cool_CO_tab_n2=40, cool_CO_tab_n3=40
type(fit3d_data(n1=cool_CO_tab_n1, n2=cool_CO_tab_n2, n3=cool_CO_tab_n3))::cool_CO_tab_data


real*8::kall(nreactions)
real*8::kall_heat(nreactions)
real*8::photo_xsecs(nphoto, nreactions)
real*8::energy_threshold(nreactions)

real*8::radial_Ncol_H2
real*8::radial_Ncol_CO
real*8::vertical_Ncol_CO

real*8::energy(nphoto)
real*8::delta_energy(nphoto-1)

real*8::ode_atol(nspecies+1)
real*8::ode_rtol(nspecies+1)

real*8::gamma_ad, d2g, user_Av, user_cr, ortho_to_para
real*8::chi_FUV  ! habing flux in range 912-1100 AA
real*8::rho_gas, rho_dust  ! gas and dust mass densities, g/cm3, do not change during integration


!! PREPROCESS_RADIATION_CONSTANTS
!! PREPROCESS_END

!! PREPROCESS_KABS_INTEGRAL
!! PREPROCESS_END

real*8,parameter::d2g_min=1d-8  ! below this limit no dust cooling/heating are calculated

real*8,parameter::hplanck=6.6260755d-27  ! erg * s
real*8,parameter::kboltzmann=1.380658d-16  ! erg / K
real*8,parameter::erg2ev=6.24150647996d11  ! 1 erg in eV
real*8,parameter::ev2erg=1d0 / erg2ev  ! 1 eV in erg
real*8,parameter::pmass=1.6726219d-24  ! proton mass in g
real*8,parameter::pi=acos(-1d0)

logical::solve_thermo
logical::solve_chemistry
logical::kall_done

character(len=1024)::runtime_data_folder
character(len=128)::reactions_verbatim(nreactions)

! FIXME: could cause problems with OPENMP
real*8::jflux_common(nphoto), log_Eabsorption

!! PREPROCESS_MASSES
!! PREPROCESS_END

end module prizmo_commons
