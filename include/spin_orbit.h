#ifndef SPIN_ORB_H
#define SPIN_ORB_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_roots.h>

#include "dynamical_system.h"
#include "aux_vmo.h"

/**
 * vector fields and jacobians related 
 * to the spin-orbit dynamics
**/

// two-body problem
int field_two_body  (double t,
                    const double y[],
                    double f[], 
                    void *params);

int jacobian_two_body   (double t,
                        const double y[],
                        double *dfdy,
                        double dfdt[],
                        void *params);

// rigid secondary body
// orbital motion is integrated directly
int field_rigid (double t,
                const double y[],
                double f[],
                void *params);

int jacobian_rigid  (double t,
                    const double y[],
                    double *dfdy,
                    double dfdt[],
                    void *params);

// rigid secondary body
// orbital motion is determined using kepler equation
int field_rigid_kepler  (double t,
                        const double y[],
                        double f[],
                        void *params);

int jacobian_rigid_kepler   (double t,
                            const double y[],
                            double *dfdy,
                            double dfdt[],
                            void *params);

// field considering a linear tidal torque
int field_linear(double t,
                const double y[],
                double f[],
                void *params);

// field considering an averaged linear tidal torque
int field_linear_average(double t,
                        const double y[],
                        double f[],
                        void *params);

struct root_params_kepler
{
	double e, t;
};

double root_function_kepler (double u,
                            void *params);

double root_derivative_kepler   (double u,
                                void *params);

void root_fdf_kepler    (double u,
                        void *params,
                        double *y,
                        double *dy);

double kepler_equation  (double e,
                        double t);

/**
 * conserved functions
**/

double angular_momentum_two_body(double y[4]);

double vis_viva_two_body    (double y[4]);

/**
 * initializers
**/

dynsys init_two_body(void *params);

dynsys init_rigid(void *params);

dynsys init_rigid_kepler(void *params);

dynsys init_linear(void *params);

dynsys init_linear_average(void *params);

int init_orbital(double y[4],
                double e);

/**
 * implementation
**/

int trace_orbit_map (double *ic,
                    dynsys system,
                    anlsis analysis);

int phase_space    (dynsys system,
                    anlsis analysis);

int draw_phase_space    (dynsys system);

int mean_resonance  (dynsys system,
                    anlsis analysis);


double dist_from_ref(double x[2],
                    double ref[2]);

int evolve_basin(double *ic, double *ref, bool *converged,
                 double ***orbit, int *orbit_size,
                 dynsys system, anlsis analysis);

int evolve_basin_union  (double *ic, double ref[][2], 
                        int num_of_basins, bool *converged,
                        double ***orbit, int *orbit_size,
                        dynsys system, anlsis analysis);

int basin_of_attraction (double *ref,
                        dynsys system,
                        anlsis analysis);

int union_basin_of_attraction   (double ref[][2],
                                int num_of_basins,
                                dynsys system,
                                anlsis analysis);

int multiple_basins_of_attraction   (dynsys system,
                                    anlsis analysis);

int all_basins_of_attraction(dynsys system,
                            anlsis analysis);

#endif