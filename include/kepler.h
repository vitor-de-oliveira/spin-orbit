#ifndef KEPLER_H
#define KEPLER_H

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

#include "auxiliar_functions_vmo.h"

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