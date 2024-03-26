# PRIZMO
## Preprocess
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

### Changing radiation or dust properties
When you change radiation or dust properties it is recommended to delete the contents of the runtime_data folder!     
This applies also if you experience weired behaviours during the runtime stage.  


## Compile and run
### Fortran
The first example is `main.f90`, it is written in FORTRAN and simulates a static disk    

```
make
./test
```
The compiler automatically searches for Intel Fortran, otherwise uses `gfortran`.    

### C
The example test `main_c.c` evolves a single cell, and it is intentended to show how to call PRIZMO from C.
```
make cbind
./test
```

### Knwon bugs/errors/warnings
#### Missing XUVTOP warning
Ignore, unless you want to produce new cooling tables using CHIANTI.    
Download the latest CHIANTI release [https://www.chiantidatabase.org/chianti_download.html].     
Export the path of where you unzipped the data (the folder containing the atom folders), e.g.
```
export XUVTOP=~/chianti
```
#### Early segmentaion fault    
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
Ignore, if you don't have any clear strategy on how to improve the convergence (e.g. producing finer tables, changing tolerances). 


#### MAXSTEPS warning
Message `WARNING: MAXSTEPS with oscillating solution`.   
Ignore, the solver is oscillating around the thermochemical equilibrium, hence it stops the integration earlier to avoid useless calculations.   




