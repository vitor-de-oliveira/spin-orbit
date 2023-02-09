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

dynsys copy_dynsys(dynsys system);

// System analysis
typedef struct AnlSis{

    // orbit evolution
    int nc, nv;
	int number_of_cycles;
	double cycle_period;
    double evolve_box_size;         // maximum size of a variable
	double coordinate_min;
	double coordinate_max;
	double velocity_min;
	double velocity_max;

    // grid variables
    int grid_resolution;
    double grid_coordinate_min;
	double grid_coordinate_max;
	double grid_velocity_min;
	double grid_velocity_max;

    // basins of attraction    
    int spin_period_min;
    int spin_period_max;
    int orbit_period_min;
    int orbit_period_max;
    int evolve_basin_time_tol;      // time close to the reference for which we say an orbit converged
    double evolve_basin_eps;        // distance from reference for which we say an orbit converged

    // basin entropy
    int sqrt_orbits_on_box;         // square root of the number of orbits on a box
    
    // time series
    int number_of_time_series;      // number of ICs for multiple time series
    double time_series_delta;       // distance between ICs for multiple time series

    // periodic orbits
    int po_max_step;                // maximum numbers of steps for the iterative method
    double po_tol;                  // maximum error allowed for the periodic orbit calculation

    // monte carlo simulation
    int number_of_rand_orbits;
    int convergence_window;
    double convergence_precision;

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
// fixes the end positions
int grid_to_double  (int grid[2], 
                     double x[2],
                     anlsis analysis);

// returns the position of a double on a grid
int double_to_grid_v2   (int grid[2], 
                         double x[2],
                         anlsis analysis);

// returns the middle double on a grid square
int grid_to_double_v2   (int grid[2], 
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