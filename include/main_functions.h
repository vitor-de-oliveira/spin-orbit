#ifndef MAIN_FUNC_H
#define MAIN_FUNC_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "dynamical_system.h"
#include "aux_vmo.h"

// calculates and writes an orbit with given initial 
// condition in the rotational system
// also calculates and writes orbit on poincare map
int trace_orbit_map(double *ic, 
                    dynsys system,
                    anlsis analysis);

// calculates and writes the phase space  
// for a given jacobi constant value
int draw_phase_space(dynsys system,
                     anlsis analysis);

#endif