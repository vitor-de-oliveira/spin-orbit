/**
 * This lib calculates periodic orbits in autonomous
 * 3d flows and 2d maps using the Levenberg–Marquardt
 * algorithm and was adapted from a software called
 * Automan written in c++ by David Ciro
**/

#ifndef PER_ORB_H
#define PER_ORB_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "aux_vmo.h"
#include "dynamical_system.h"

// I have set a limit of 100 for the orbit period
typedef struct PerOrb{
    int     period;
    double  initial_condition[2];
    double  *orbit[2];
} perorb;

// numerically approximate jacobian on a periodic orbit
void jacobian_periodic_orbit(double **J,
                             perorb po,
                             dynsys system, 
                             anlsis analysis);

// minimization step for Levenberg-Marquardt method
void minimization_step  (double *x,
                         double lamb,
                         perorb po,
                         dynsys system, 
                         anlsis analysis);

// calculates the periodic orbit
int periodic_orbit  (double seed[2],
                     perorb po,
                     dynsys system, 
                     anlsis analysis);

// fills in orbit
int fills_periodic_orbit(perorb po,
                         dynsys system, 
                         anlsis analysis);

/**
 * additional dynamical systems functions
**/

// evolves an initial condition for n cycles
// considers t = 0
int evolve_n_cycles (int    n,
                     double *y,
                     dynsys system,
                     anlsis analysis);

/**
 * additional mathematical functions
**/

// solves Ax = b for x using Gaussian elimination
void gauss_solve (double **A, double *x, double *b, int dim);

#endif
