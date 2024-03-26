Verner, D. A., & Ferland, G. J. 1996, ApJS, in press. 
Atomic data for astrophysics. I. Radiative recombination rates for 
H-like, He-like, Li-like and Na-like ions over a broad range of temperature. 

                    DESCRIPTION OF THE TABLES

Table 1. Fit parameters for radiative recombination rates (table1.dat,
6048 byte).
Byte-per-byte description of file: table1.dat (108 lines)
-------------------------------------------------------------------------------
   Bytes Format  Units   Label        Explanations
-------------------------------------------------------------------------------
   1-  8  A8     ---     Ion          Ion label    
  10- 11  I2     ---     Z            Atomic number
  13- 14  I2     ---     N            Number of electrons
  17- 25  E9.3   cm^3/s  a            Fit parameter, Eq. (4)
  28- 33  F6.4   ---     b            Fit parameter, Eq. (4)
  36- 44  E9.3   K       T_0          Fit parameter, Eq. (4)
  47- 55  E9.3   K       T_1          Fit parameter, Eq. (4)
-------------------------------------------------------------------------------
Note that there are two entries for HeI:
  a - fit is valid at 3 K < T < 10^6 K, rms error 2.5%
  b - fit is valid at 3 K < T < 10^10 K, rms error 4.7%

Table 2. The partial radiative recombination rate coefficients for He I
(he1.dat, 8065 byte).
The partial radiative recombination rate coefficients,
in cm^3/s, for the ground state and 40 excited states of He I at 
20 temperatures, from log T = 0.5 K to log T = 10.0 K. 
Format: a11,20(1x,e8.2). The rates are calculated by use of the 
Opacity Project photoionization cross sections (Fernley, Taylor, & Seaton,
1987, J. Phys. B., 20, 6457). Note that the higher excited levels are
also contribute to the total rates, especially at low temperatures. 

Table 3. The partial radiative recombination rate coefficients for He-like
ions (he-like.dat, 71065 byte).
The partial radiative recombination rate 
coefficients, in cm^3/s, for the He-like ions at 20 temperatures, 
from log T = 0.5 K to log T = 10.0 K. 
Format: i2,1x,i2,1x,a2,20(1x,e8.2). The rates are calculated by use of the 
fits to the photoionization cross sections (ground states: Verner et al.,
1995, in preparation; excited states: Clark, Cowan, & Bobrowicz, 1986,
Atomic Data Nucl. Data Tables, 34, 415). Note that the higher excited levels 
are also contribute to the total rates. 
Byte-per-byte description of file: he-like.dat
-------------------------------------------------------------------------------
Bytes   Format      Units  Label   Explanations
-------------------------------------------------------------------------------
1-  2     I2         ---     Z     Atomic number
4-  5     I2         ---     N     Number of electrons
7-  8     A2         ---     Sh    Subshell label
9-188 20(1x,E8.2)  cm^3/s    --    Partial radiative recombination rates
-------------------------------------------------------------------------------

Table 4. The partial radiative recombination rate coefficients for Li-like
ions (li-like.dat, 63694 byte).
The same as Table 3, fo Li-like ions.

Table 5. The partial radiative recombination rate coefficients for Na-like
ions (na-like.dat, 36478 byte).
The same as Table 3, fo Na-like ions.

Note: Tables 2-5 are available in electronic form only.


