#include <stdio.h>

// interfaces
extern void prizmo_init_c();
extern void prizmo_evolve_c(double *, double *, double *, double *, int *, int *);

#define NTRACER 33
#define NPHOTO 1000
#define NFLX 0
#define NIONS 0

#define IDX_CHEM_C (NFLX + NIONS + 0)
#define IDX_CHEM_CH (NFLX + NIONS + 1)
#define IDX_CHEM_CH2 (NFLX + NIONS + 2)
#define IDX_CHEM_CH2j (NFLX + NIONS + 3)
#define IDX_CHEM_CH3 (NFLX + NIONS + 4)
#define IDX_CHEM_CH3j (NFLX + NIONS + 5)
#define IDX_CHEM_CH4 (NFLX + NIONS + 6)
#define IDX_CHEM_CH4j (NFLX + NIONS + 7)
#define IDX_CHEM_CH5j (NFLX + NIONS + 8)
#define IDX_CHEM_CHj (NFLX + NIONS + 9)
#define IDX_CHEM_CO (NFLX + NIONS + 10)
#define IDX_CHEM_CO_DUST (NFLX + NIONS + 11)
#define IDX_CHEM_COj (NFLX + NIONS + 12)
#define IDX_CHEM_Cj (NFLX + NIONS + 13)
#define IDX_CHEM_E (NFLX + NIONS + 14)
#define IDX_CHEM_H (NFLX + NIONS + 15)
#define IDX_CHEM_H2 (NFLX + NIONS + 16)
#define IDX_CHEM_H2O (NFLX + NIONS + 17)
#define IDX_CHEM_H2O_DUST (NFLX + NIONS + 18)
#define IDX_CHEM_H2Oj (NFLX + NIONS + 19)
#define IDX_CHEM_H2j (NFLX + NIONS + 20)
#define IDX_CHEM_H3Oj (NFLX + NIONS + 21)
#define IDX_CHEM_H3j (NFLX + NIONS + 22)
#define IDX_CHEM_HCOj (NFLX + NIONS + 23)
#define IDX_CHEM_He (NFLX + NIONS + 24)
#define IDX_CHEM_Hej (NFLX + NIONS + 25)
#define IDX_CHEM_Hj (NFLX + NIONS + 26)
#define IDX_CHEM_O (NFLX + NIONS + 27)
#define IDX_CHEM_O2 (NFLX + NIONS + 28)
#define IDX_CHEM_O2j (NFLX + NIONS + 29)
#define IDX_CHEM_OH (NFLX + NIONS + 30)
#define IDX_CHEM_OHj (NFLX + NIONS + 31)
#define IDX_CHEM_Oj (NFLX + NIONS + 32)

int main (void){
  double x[NTRACER] = {0e0};
  double jflux[NPHOTO] = {0e0};
  double ngas, tgas, dt, t, tend, spy;
  int ierr, verb;

  prizmo_init_c();
  spy = 3.1536e7;

  tgas = 1e3;
  ngas = 1e4;
  tend = spy * 1e8;

  x[IDX_CHEM_H2] = 1e0;
  x[IDX_CHEM_Cj] = 1e-4;
  x[IDX_CHEM_E] = 1e-4;
  x[IDX_CHEM_O] = 2e-4;
  x[IDX_CHEM_He] = 0.08;

  for(int i = 0; i < NTRACER; i++){
    x[i] *= ngas;
  }

  dt = spy;
  verb = 0;

  t = 0e0;

  FILE *f = fopen("output.dat", "w");

  while (1){
    dt *= 1.2;
    prizmo_evolve_c(x, &tgas, jflux, &dt, &verb, &ierr);
    t += dt;

    fprintf(f, "%e ", t);
    fprintf(f, "%e ", tgas);
    for(int i = 0; i < NTRACER; i++){
      fprintf(f, "%e ", x[i]);
    }
    fprintf(f, "\n");

    if(t > tend){
      break;
    }
  }

  fclose(f);

}