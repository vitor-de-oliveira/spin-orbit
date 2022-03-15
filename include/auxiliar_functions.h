#ifndef AUX_FUNC_H
#define AUX_FUNC_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "dynamical_system.h"
#include "auxiliar_functions_vmo.h"

/****************** high-level functions *************************/

// calculates an orbit with given initial condition in the rotational system
// using gsl with rk8pd as the integrator of choice
// defines passed arrays with the complete orbit evenly spaced in time
int evolve_cycle(double *y, void *params, double cycle_period, double *t);

int evolve_orbit(void *params, double *ic, double cycle_period, int number_of_cycles, double ***orbit);

/******************* zero-velocity curve *************************/

// field for tracing zero velocity curves
// on the rotational frame
int field_zero_velociy_curves_rotational(double t, const double y[], double f[], void *params);

// structure with parameters for root calculation
struct root_params
{
	double mu, J;
};

// function to be solved for root calculation
// with seed close to the third Lagrangian point
double root_function_L3(double x, void *params);

// derivative of function to be solved
// for root calculation
// with seed close to the third Lagrangian point
double root_derivative_L3(double x, void *params);

// both function and derivative to be solved
// for root calculation
// with seed close to the third Lagrangian point
void root_fdf_L3(double x, void *params, double *y, double *dy);

// function to be solved for root calculation
// with seed close to the fourth and fifth Lagrangian points
double root_function_triangular_points(double y, void *params);

// derivative of function to be solved
// for root calculation
// with seed close to the fourth and fifth Lagrangian points
double root_derivative_triangular_points(double y, void *params);

// both function and derivative to be solved
// for root calculation
// with seed close to the fourth and fifth Lagrangian points
void root_fdf_triangular_points(double y, void *params, double *z, double *dz);

// determines the root for zero velocity curve tracing
int root_zvc(double J, void *params, double *y, double x);

#endif