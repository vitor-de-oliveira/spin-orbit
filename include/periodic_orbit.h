/**
 * This lib calculates periodic orbits in autonomous
 * 3d flows and 2d maps using the Levenberg–Marquardt
 * algorithm and was adapted from a software called
 * Automan written by David Ciro
**/

#ifndef PER_ORB_H
#define PER_ORB_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "aux_vmo.h"
#include "dynamical_system.h"

typedef struct PerOrb{
    int     period;
    double  seed[2];
    double  initial_condition[2];
    double  **orbit;

    // absolute value of both eigenvalues
    double  eigenvalues_absolute_value[2];

    // winding number
    int     winding_number_numerator;
    int     winding_number_denominator;

    // evolves an initial condition over n cycles
    int     (*evolve_n_cycles)  (double y0[2],
                                 int n,
                                 dynsys system,
                                 anlsis analysis,
                                 rngkta rk);

    // distance between points on 2d phase space
    double  (*dist_on_phase_space)  (double x[2], 
                                     double y[2]);

} perorb;

// numerically approximate jacobian on a periodic orbit
void jacobian_periodic_orbit(double **J,
                             perorb po,
                             dynsys system, 
                             anlsis analysis,
                             rngkta rk);

// minimization step for Levenberg-Marquardt method
void minimization_step  (double *x,
                         double *err,
                         double lamb,
                         perorb po,
                         dynsys system, 
                         anlsis analysis,
                         rngkta rk);

// calculates the periodic orbit's IC
int calculate_periodic_orbit_ic(perorb *po,
                                dynsys system, 
                                anlsis analysis,
                                rngkta rk);

// magnitude of eigenvalues for the jacobian
void jacobian_eigenvalues_magnitude  (perorb *po,
                                      dynsys system, 
                                      anlsis analysis,
                                      rngkta rk);                                

/**
 * mathematical functions
**/

// solves Ax = b for x using Gaussian elimination
void gauss_solve (double *x, double **A, double *b, int dim);

#endif
