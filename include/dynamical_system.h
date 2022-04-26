#ifndef DYN_SYS_H
#define DYN_SYS_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

#include "aux_vmo.h"

// field for rotational frame
int field_two_body(double t, const double y[], double f[], 
                   void *params);

// jacobian for rotational frame
int jacobian_two_body(double t, const double y[], 
            double *dfdy, double dfdt[], void *params);

// field for rotational frame
int field_rigid(double t, const double y[], double f[], 
            void *params);

// jacobian for rotational frame
int jacobian_rigid(double t, const double y[], double *dfdy, 
            double dfdt[], void *params);

// field for rotational frame
int field_rigid_kepler(double t, const double y[], 
            double f[], void *params);

// jacobian for rotational frame
int jacobian_rigid_kepler(double t, const double y[], 
            double *dfdy, double dfdt[], void *params);

// field for rotational frame
double angular_momentum_two_body(double y[4]);

double vis_viva_two_body(double y[4], double T, double a);

int init_orbital(double y[4], double e);

// Dynamical System
typedef struct DynSys{
    char *name;
    int (*field) (double t, const double y[], double dydt[], void *params);
    int (*jac) (double t, const double y[], double *dfdy, double dfdt[],
                   void *params);
    size_t dim;
    void *params;
} dynsys;

dynsys init_rigid(void *params);

dynsys init_rigid_kepler(void *params);

dynsys init_two_body(void *params);

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

/**
 * functions to handle GSL
**/

// sets driver for chosen error and control
int set_driver(gsl_odeiv2_driver **d, 
               gsl_odeiv2_system *sys, 
			   const gsl_odeiv2_step_type *T, 
			   double h, double h_max, double h_min,
			   double error_abs, double error_rel,
			   char *control);

// sets the chosen integrator
int set_integrator(const gsl_odeiv2_step_type **T, 
                   char *integrator);

// sets the chosen system
int set_system(gsl_odeiv2_system *sys, 
               dynsys system);

// structure with parameters for root calculation
struct root_params
{
	double e, t;
};

// function to be solved for root calculation
// with seed close to the third Lagrangian point
double root_function_kepler(double u, void *params);

// derivative of function to be solved
// for root calculation
// with seed close to the third Lagrangian point
double root_derivative_kepler(double u, void *params);

// both function and derivative to be solved
// for root calculation
// with seed close to the third Lagrangian point
void root_fdf_kepler(double u, void *params, double *y, double *dy);

// maybe move it or change lib name
double kepler_equation(double e, double t);

#endif