#executable name
exe = test

-include Makefile_pragmas

# test if ifort is present
wres = $(shell which ifort > /dev/null; echo $$?)
ifeq "$(wres)" "0"
	fc = ifort
	switchOPT = -O3 -g -traceback -xHost -fp-model=precise -fpp
	#switchOPT += -no-prec-sqrt
	switchDBG = -O0 -check all -warn all -fpp -save-temps
	switchDBG += -fpe0 -u -traceback -warn nounused -g
	switchDBG += -init=snan,zero,arrays -ftrapuv -check noarg_temp_created $(pragmas)
	switchOMP = -qopenmp
	nowarn = -nowarn
else
	fc = gfortran
	switchOPT = -ffree-line-length-none -O3 -g
	switchDBG = -fbacktrace -g
	switchDBG += -ffpe-trap=zero,overflow,invalid
	switchDBG += -fbounds-check -ffree-line-length-none -O0
	nowarn = -w
endif

# test if icc is present
wres_c = $(shell which icc > /dev/null; echo $$?)
ifeq "$(wres_c)" "0"
	cc = icc
	switchCLIB = -lifcore -limf -lsvml -lintlc -ldl
	switchCDBG = -O0
	cswitch =
else
	cc = gcc
	switchCLIB = -lgfortran -lm
	switchCDBG = -O0
	cswitch =
endif

#default switch
fswitch = $(switchOPT)

#objects
objs = opkda1.o
objs += opkda2.o
objs += opkdmain.o
objs += prizmo_fit.o
objs += prizmo_linear_solver.o
#objs += prizmo_spline.o
objs += prizmo_commons.o
objs += prizmo_utils.o
objs += prizmo_loaders.o
objs += prizmo_shielding.o
objs += prizmo_rates.o
objs += prizmo_rates_photo.o
objs += prizmo_rates_heating.o
objs += prizmo_flux.o
objs += prizmo_tdust.o
objs += prizmo_heating_photo.o
objs += prizmo_heating_photoelectric.o
objs += prizmo_heating_H2diss.o
objs += prizmo_heating_CR.o
objs += prizmo_heating.o
objs += prizmo_cooling_H2.o
objs += prizmo_cooling_CO.o
objs += prizmo_cooling_atomic.o
objs += prizmo_cooling_dust.o
objs += prizmo_cooling_chemical.o
objs += prizmo_cooling.o
objs += prizmo_ode.o
objs += prizmo_attenuate.o
objs += prizmo_core.o
objs += prizmo.o

obj_main = main.o
obj_main_c = main_c.o

#default target
all:	$(objs) $(obj_main)
	$(fc) $(fswitch) $(objs) $(obj_main) -o $(exe) $(lib)

#full debug target
debug: fswitch = $(switchDBG)
debug: all

profile: fswitch = $(switchOPT) -pg
profile: all

cbind:	$(objs) prizmo_c.o $(obj_main_c)
	$(cc) $(objs) prizmo_c.o $(obj_main_c) -o $(exe) $(switchCLIB) $(lib)

cbind_debug: fswitch = $(switchDBG)
cbind_debug: cswitch = $(switchCDBG)
cbind_debug: cbind

cbind_profile: fswitch = -O0 -pg
cbind_profile: cswitch = -O0 -pg
cbind_profile: cbind

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ *.i90 *.i $(exe)

.PHONY: clean

#rule for f90
%.o:%.f90
	$(fc) $(fswitch) -c $^ -o $@ $(lib)

#rule for f
%.o:%.f
	$(fc) $(fswitch) -c $^ -o $@ $(lib) $(nowarn)

#rule for c
%.o:%.c
	$(cc) $(cswitch) -c $^ -o $@
