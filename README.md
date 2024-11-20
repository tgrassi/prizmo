# PRIZMO
- PRIZMO is a library-like code to advance time-dependetly chemistry and temperature of a protoplanetary disk (M)HD simulation.
- It preprocesses the input information to write optimized FORTRAN code.
- It has a C interface that allows it to be coupled with codes like PLUTO.
- The earlier code is described in [https://arxiv.org/abs/2004.04748] (Grassi et al. 2020)
- The newer version is discussed in [https://arxiv.org/abs/2408.00848] (Sellek et al. 2024)

![plot](./assets/disk.png)

## Basic usage
- Run the preprocessor
```
cd src_py
python prizmo.py
```
- Compile the code
```
cd ..
make
```
- Run the code
```
./test
```


## Preprocessor
Preprocessing is needed:

* After cloning this repository
* When changing chemistry (i.e. chemical network)
* After changing data in the `data` folder
* When changing dust or radiation parameters in `src_py/prizmo_commons.py`

How to preprocess:
```
cd src_py
python prizmo.py
```
wait...

### Changing preprocessor inputs (e.g. radiation or dust properties)
PRIZMO's preprocessor has default values for many choices. However, these can be configured by the user either by passing an input file with the flag -i (see test.ini for an example), or by setting any of the following command line arguments directly:
* chemNet - the chemical network specified as a list of reactions
* atomData - the file containing the details of level energies and fits for the de-exciation rates for the atomic cooling
* radiation_type - details of the spectrum to use
* nphoto - the number of energy bins to use
* energy_minmax - the minimum and maximum energies to use (eV)
* dust_minmax - the minimum and maximum dust grain sizes to use (cm)
* refInd_file - the file containing the refractive indices for the dust
The command line arguments that were used are logged in a readme file in the runtime_data folder

Note: When you change radiation or dust properties, it is recommended that the contents of the runtime_data folder be deleted!
This also applies if you experience weird behaviors during the runtime stage.

## Compile and run
### Fortran
The first example is `main.f90`, it is written in FORTRAN and simulates a static disk    

```
make
./test
```
The makefile automatically searches for Intel Fortran, otherwise uses `gfortran`.    

### C
The example test `main_c.c` evolves a single cell, and it is intended to show how to call PRIZMO from C.
```
make cbind
./test
```

## Known bugs/errors/warnings
#### Missing XUVTOP warning
Ignore it unless you want to produce new cooling tables using CHIANTI.    
Download the latest CHIANTI release [https://www.chiantidatabase.org/chianti_download.html].     
Export the path of where you unzipped the data (the folder containing the atom folders), e.g.
```
export XUVTOP=~/chianti
```
#### Early segmentation fault    
Segmentation fault at the beginning when running `./test`.       
PRIZMO uses large tables (especially atomic cooling), hence you need to increase the stack size to 
```
ulimit -s unlimited
```
#### DLSODES warnings
Warnings similar to this, but the code continues to run    
```
 DLSODES- At current T (=R1), MXSTEP (=I1) steps             
       taken on this call before reaching TOUT               
      In above message,  I1 =     10000
      In above message,  R1 =  0.5040993677130D+03
```
The solver is taking too many iterations to advance, but the solution is found anyway.    
The cell is probably close to thermochemical equilibrium.     
Ignore, if you don't have any clear strategy on how to improve the convergence (e.g., producing finer tables, changing tolerances). 


#### MAXSTEPS warning
Message `WARNING: MAXSTEPS with oscillating solution`.   
Ignore, the solver is oscillating around the thermochemical equilibrium, hence it stops the integration earlier to avoid useless calculations.   




