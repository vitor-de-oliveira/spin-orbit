#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>

#include "aux_vmo.h"

// Dynamical System
typedef struct DynSys{
    char *name;
    int (*field)    (double t, 
                    const double y[], 
                    double dydt[], 
                    void *params);
    int (*jac)  (double t,
                const double y[],
                double *dfdy,
                double dfdt[],
                void *params);
    size_t dim;
    void *params;
} dynsys;

// System analysis
typedef struct AnlSis{
    size_t nc, nv;
	int number_of_cycles;
	double cycle_period;
	double coordinate_min;
	double coordinate_max;
	double velocity_min;
	double velocity_max;
} anlsis;

int evolve_cycle(double *y, double cycle_period, 
                 double *t, dynsys system);

int evolve_orbit(double *ic, double cycle_period, 
                 int number_of_cycles, double ***orbit, 
                 int *orbit_size, dynsys system);

int trace_orbit_map(double *ic, 
                    dynsys system,
                    anlsis analysis);

int draw_phase_space(dynsys system,
                     anlsis analysis);

/**
 * functions to handle GSL
**/

int set_driver(gsl_odeiv2_driver **d, 
               gsl_odeiv2_system *sys, 
			   const gsl_odeiv2_step_type *T, 
			   double h, double h_max, double h_min,
			   double error_abs, double error_rel,
			   char *control);

int set_integrator(const gsl_odeiv2_step_type **T, 
                   char *integrator);

int set_system(gsl_odeiv2_system *sys, 
               dynsys system);

#endif