#include <stdio.h>
#include <math.h>
#include <time.h>

// interfaces
extern void prizmo_init_c();
extern void prizmo_evolve_c(double *, double *, double *, double *, int *, int *);
extern void prizmo_rt_c(double *, double *, double *, double *);
extern void prizmo_set_radial_ncol_h2_c(double *);
extern void prizmo_set_radial_ncol_co_c(double *);
extern void prizmo_set_vertical_ncol_co_c(double *);
extern void prizmo_get_tdust_c(double *, double *, double *, double *);
extern void prizmo_load_radiation_field_c(double *);

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
#define ntheta 24

int main(void)
{
    // define variables
    double x[nr][ntheta][NTRACER] = {0e0};
    double ngas[nr][ntheta];
    double tgas[nr][ntheta];
    double tdust[nr][ntheta];
    double cpu[nr][ntheta];
    double jflux[NPHOTO];
    double jflux0[NPHOTO];
    double dt, t, r[nr], theta[ntheta], ds, dtheta;
    int ierr, verb;

    // constants
    double spy = 3.1536e7;
    double au2cm = 1.496e13;
    double pmass = 1.67e-24;        // g
    double msun = 1.989e33;         // g
    double rsun = 6.96e10;         // cm
    double gravity = 6.67e-8;       // cm^3/g/s^2
    double kboltzmann = 1.3807e-16; // erg/K

    // disk parameters
    double rmin = 1e0 * au2cm; // min radius, cm
    double rmax = 1e2 * au2cm; // max radius, cm
    double tend = 1e5 * spy; // max integration time, s
    double sigma0 = 1700.; // surface density at 1AU, g/cm^2
    double tgas0 = 1e2;    // temperature at 1AU, K
    double mstar = msun;    // star mass, g
    double rstar = rsun; // star radius, cm
    double mu = 2.34;  // mean molecular weight, no units
    double ngas_min = 1e-2; // min density, 1/cm^3

    // -----------------------------
    // initialize prizmo
    prizmo_init_c();

    // -----------------------------
    // Initialize the gas density and temperature
    // loop over radius
    for (int i = 0; i < nr; i++)
    {
        // compute radius
        r[i] = rmin + (rmax - rmin) * i / (nr - 1);
        // compute orbital frequency
        double omega = sqrt(gravity * mstar / pow(r[i], 3));

        // loop over theta (note theta=0 is the pole)
        for (int j = 0; j < ntheta; j++)
        {
            theta[j] = M_PI / 4.0 + M_PI / 4.0 * j / (ntheta - 1);
            tgas[i][j] = tgas0 * pow(r[i] / au2cm, -0.5);
            double cs = sqrt(kboltzmann * tgas[i][j] / mu / pmass);
            double ng = sigma0 * pow(r[i] / au2cm, -1.5) / pmass / mu;
            double h = cs / omega;
            double z = r[i] * tan(M_PI / 2.0 - theta[j]);
            ngas[i][j] = ng / h / sqrt(2e0 * M_PI) * exp(-0.5 * pow(z / h, 2));
            ngas[i][j] = fmax(ngas[i][j], ngas_min);
        }
    }

    // -----------------------------
    // Initialize the chemical abundances all grid cells
    // loop over radius
    for (int i = 0; i < nr; i++)
    {
        // loop over theta (note theta=0 is the pole)
        for (int j = 0; j < ntheta; j++)
        {
            // fractional abundances wrt total gas density
            x[i][j][IDX_CHEM_H2] = 1e0;
            x[i][j][IDX_CHEM_Cj] = 1e-4;
            x[i][j][IDX_CHEM_E] = 1e-4;
            x[i][j][IDX_CHEM_O] = 2e-4;
            x[i][j][IDX_CHEM_He] = 0.08;
            // multiply by gas density to get number density, 1/cm^3
            for (int k = 0; k < NTRACER; k++)
            {
                x[i][j][k] *= ngas[i][j];
            }
        }
    }

    // -----------------------------
    // Initialize the radiation field read from file
    prizmo_load_radiation_field_c(jflux0);

    // -----------------------------
    // time integration
    dt = spy * 1e-3;  // initial time step, s
    t = 0e0; // initial time, s
    verb = 0;  // verbosity level, 0=none
    ds = r[1] - r[0];  // cell size, cm
    dtheta = theta[1] - theta[0];  // cell size, rad

    // loop over time
    do
    {

        // compute the time step
        dt = fmin(dt * 1.5, tend - t);
        t += dt;
        double ncol_co_vertical[nr] = {0e0}; // column density of CO, from pole to equator, cm^-2

        // loop over theta (note theta=0 is the pole), from pole to equator
        for (int j = ntheta - 1; j >= 0; j--)
        {
            double rold = rmin;
            double ncol_h2 = 0e0; // column density of H2, cm^-2
            double ncol_co = 0e0; // column density of CO, cm^-2

            // initialize the flux to the star radiation field
            for(int k=0; k<NPHOTO; k++)
            {
                jflux[k] = jflux0[k] * 4. * M_PI * M_PI * rstar * rstar / rmin / rmin;
            }

            // loop over radius
            for (int i = 0; i < nr; i++)
            {
                // print out the progress
                printf("\r t=%.1f%%; theta=%.1f%%; r=%.1f%%;    ",
                    1e2 * log10(t+1e-40) / log10(tend),
                    1e2 * (1e0 - (float) j / (ntheta - 1)),
                    1e2 * i / (nr - 1));

                // compute geomeric dilution of the radiation field
                for(int k=0; k<NPHOTO; k++)
                {
                    jflux[k] *= rold * rold / r[i] / r[i];
                }

                // set the column densities
                prizmo_set_radial_ncol_h2_c(&ncol_h2);
                prizmo_set_radial_ncol_co_c(&ncol_co);

                // set the vertical column density of CO
                prizmo_set_vertical_ncol_co_c(&ncol_co_vertical[i]);

                // set start clock time
                cpu[i][j] = clock();

                // compute the chemical evolution
                // x: abundances, 1/cm^3, array of size NTRACER, input/output
                // tgas: gas temperature, K, input/output
                // jflux: radiation field, erg/cm^2/s, array of size NPHOTO, input
                // dt: time step, s, input
                // verb: verbosity level, 0=none, input
                // ierr: error code, 0=ok, output
                prizmo_evolve_c(x[i][j], &tgas[i][j], jflux, &dt, &verb, &ierr);

                // compute elapsed time
                cpu[i][j] = ((double) (clock() - cpu[i][j])) / CLOCKS_PER_SEC;

                // get the dust temperature for output reasons only
                prizmo_get_tdust_c(x[i][j], &tgas[i][j], &tdust[i][j], jflux);

                // radiative transfer
                // x: abundances, 1/cm^3, array of size NTRACER, input
                // tgas: gas temperature, K, input
                // jflux: radiation field, erg/cm^2/s, array of size NPHOTO, input/output
                // ds: cell size, cm, input
                prizmo_rt_c(x[i][j], &tgas[i][j], jflux, &ds);

                // update the radial column densities of H2 and CO
                ncol_h2 += x[i][j][IDX_CHEM_H2] * ds;
                ncol_co += x[i][j][IDX_CHEM_CO] * ds;

                // compute approximately the height of the cell
                double dz = r[i] * (tan(M_PI / 2.0 - theta[j]) - tan(M_PI / 2.0 - theta[j] - dtheta));

                // update the vertical column density of CO
                ncol_co_vertical[i] += x[i][j][IDX_CHEM_CO] * dz;

                // store the radius for geometric dilution
                rold = r[i];
            }
        }

        // Write the results to a file
        FILE *f = fopen("output.dat", "w");
        for (int j = 0; j < ntheta; j++)
        {
            for (int i = 0; i < nr; i++)
            {
                double z = r[i] * tan(M_PI / 2.0 - theta[j]);
                fprintf(f, "%e %e %e %e %e %e %e ", r[i], theta[j], z, ngas[i][j], tgas[i][j], tdust[i][j], cpu[i][j]);
                for (int k = 0; k < NTRACER; k++)
                {
                    fprintf(f, "%e ", x[i][j][k]);
                }
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
        fclose(f);
    } while (t < tend);
    printf("\n");

    return 0;
}