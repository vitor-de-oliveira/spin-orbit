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
#include "auxiliar_functions_gsl.h"

/****************** high-level functions ********************/

// calculates an orbit with given initial condition in the rotational system
// using gsl with rk8pd as the integrator of choice
// defines passed arrays with the complete orbit evenly spaced in time
int evolve_cycle(double *y, void *params, double cycle_period, 
                 double *t, dynsys system);

int evolve_orbit(void *params, double *ic, double cycle_period, 
                 int number_of_cycles, double ***orbit, 
                 int *orbit_size, dynsys system);

#endif