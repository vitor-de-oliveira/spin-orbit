#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    int nc, nv;
    int grid_resolution;
	int number_of_cycles;
	double cycle_period;
	double coordinate_min;
	double coordinate_max;
	double velocity_min;
	double velocity_max;
    double grid_coordinate_min;
	double grid_coordinate_max;
	double grid_velocity_min;
	double grid_velocity_max;
    double evolve_box_size;     // maximum size of a variable
    double evolve_basin_eps;    // distance from reference for which we say an orbit converged
    double time_series_delta;
} anlsis;

anlsis copy_anlsis(anlsis analysis);

int evolve_cycle(double *y, double *t,
                dynsys system, anlsis analysis);

int evolve_orbit(double *ic, double ***orbit, 
                 int *orbit_size, dynsys system,
                 anlsis analysis);

// returns the position of a double on a grid
int double_to_grid  (int grid[2], 
                    double x[2],
                    anlsis analysis);

// returns the first double on a grid square
int grid_to_double  (int grid[2], 
                    double x[2],
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