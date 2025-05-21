#include <stdio.h>
#include <math.h>


// interfaces
extern void prizmo_init_c();
extern void prizmo_evolve_c(double *, double *, double *, double *, int *, int *);
extern void prizmo_rt_c(double *, double *, double *, double *);

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

#define nr 20
#define ntheta 21

int main (void){
    double x[nr][ntheta][NTRACER] = {0e0};
    double ngas[nr][ntheta];
    double tgas[nr][ntheta];
    double jflux[NPHOTO] = {0e0};
    double dt, t, tend, r[nr], theta[ntheta], ds;
    int ierr, verb;

    double spy = 3.1536e7;
    double au2cm = 1.496e13;
    double pmass = 1.67e-24; // g
    double msun = 1.989e33; // g
    double gravity = 6.67e-8; // cm^3/g/s^2
    double kboltzmann = 1.3807e-16; // erg/K

    tend = spy * 1e8;

    double rmin = 1e0 * au2cm;
    double rmax = 1e2 * au2cm;
    double sigma0 = 1700.; // g/cm^2
    double tgas0 = 1e2; // K
    double mass = msun; // g
    double mu = 2.34;
    double ngas_min = 1e-2; // 1/cm^3

    prizmo_init_c();


    // Initialize the gas density and temperature
    FILE *f = fopen("disk.dat", "w");
    for(int i = 0; i < nr; i++){
        r[i] = pow(1e1, log10(rmin) + (log10(rmax) - log10(rmin)) * i / (nr - 1));
        double omega = sqrt(gravity * mass / pow(r[i], 3));

        for(int j = 0; j < ntheta; j++){
            theta[j] = M_PI / 4.0 + M_PI / 4.0 * j / (ntheta - 1);
            tgas[i][j] = tgas0 * pow(r[i] / au2cm, -0.5);
            double cs = sqrt(kboltzmann * tgas[i][j] / mu / pmass);
            double ng = sigma0 * pow(r[i] / au2cm, -1.5) / pmass / mu;
            double h = cs / omega;
            double z = r[i] * tan(M_PI / 2.0 - theta[j]);
            ngas[i][j] = ng / h / sqrt(2e0 * M_PI) * exp(-0.5 * pow(z / h, 2));
            ngas[i][j] = fmax(ngas[i][j], ngas_min);
            fprintf(f, "%e %e %e %e %e\n", r[i], theta[j], z, ngas[i][j], tgas[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);


    // Initialize the chemical abundances
    for(int i = 0; i < nr; i++){
        for(int j = 0; j < ntheta; j++){
            x[i][j][IDX_CHEM_H2] = 1e0;
            x[i][j][IDX_CHEM_Cj] = 1e-4;
            x[i][j][IDX_CHEM_E] = 1e-4;
            x[i][j][IDX_CHEM_O] = 2e-4;
            x[i][j][IDX_CHEM_He] = 0.08;
            for(int k = 0; k < NTRACER; k++){
                x[i][j][k] *= ngas[i][j];
            }
        }
    }

    dt = spy * 1e3;
    verb = 0;

    f = fopen("output.dat", "w");

    for(int j = 0; j < ntheta; j++){
        printf("%d %d\n", j, ntheta);
        double rold = 0e0;
        for(int i = 0; i < nr; i++){
            prizmo_evolve_c(x[i][j], &tgas[i][j], jflux, &dt, &verb, &ierr);
            double z = r[i] * tan(M_PI / 2.0 - theta[j]);
            fprintf(f, "%e %e %e %e %e ", r[i], theta[j], z, ngas[i][j], tgas[i][j]);
            for(int k = 0; k < NTRACER; k++){
                fprintf(f, "%e ", x[i][j][k]);
            }
            ds = r[i] - rold;
            prizmo_rt_c(x[i][j], &tgas[i][j], jflux, &ds);
            rold = r[i];
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }

    fclose(f);

    /*FILE *f = fopen("output.dat", "w");

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

    fclose(f);*/
    return 0;
}