#ifndef SPIN_ORB_H
#define SPIN_ORB_H

#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_roots.h>

#include "aux_vmo.h"
#include "dynamical_system.h"
#include "periodic_orbit.h"

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

int orbit_map (double *ic,
               dynsys system,
               anlsis analysis);

int phase_space (dynsys system,
                anlsis analysis);

int time_series(dynsys system,
                anlsis analysis);

int multiple_time_series(dynsys system,
                        anlsis analysis);

int multiple_time_series_delta_theta_dot(dynsys system,
										anlsis analysis);

int multiple_time_series_delta_theta(dynsys system,
									anlsis analysis);

double dist_from_ref(double x[2],
                     double ref[2]);

int evolve_basin(double *ic, double ref[][2], 
                int ref_period, bool *converged,
                double ***orbit, int *orbit_size,
                dynsys system, anlsis analysis);

int basin_of_attraction (double ref[][2], int ref_period,
                        dynsys system,
                        anlsis analysis);

/**
 * periodic orbit
**/

// search the phase space for a resonance
// of orbital period given by analysis.number_of_cycles
int look_for_resonance	(int orbital_period,
                         int number_of_spins,
                         dynsys system, 
						 anlsis analysis);

// evolves an initial condition for n cycles
int evolve_n_cycles_po  (double y0[2],
                         int n,
                         dynsys system,
                         anlsis analysis);

// calculates periodic orbit 
// and prints it on an exit file
int periodic_orbit	(perorb *po,
                     dynsys system,
                     anlsis analysis);

/**
 * Gnuplot pipe
**/

int draw_orbit_map  (dynsys system);

int draw_phase_space(dynsys system);

int draw_phase_space_clean(dynsys system);

int draw_phase_space_latex(dynsys system);

int draw_orbit_on_phase_space(dynsys system);

int draw_orbit_on_phase_space_latex(dynsys system);

int draw_time_series(dynsys system);

int draw_time_series_latex(dynsys system);

int draw_time_series_union_e(dynsys system);

int draw_time_series_union_e_latex(dynsys system);

int draw_time_series_union_e_eps(dynsys system);

int draw_time_series_union_K(dynsys system);

int draw_time_series_union_K_latex(dynsys system);

int draw_multiple_time_series(dynsys system);

int draw_multiple_time_series_delta_theta_dot(dynsys system,
                                              anlsis analysis);

int draw_multiple_time_series_delta_theta_dot_latex(dynsys system,
                                                    anlsis analysis);

int draw_multiple_time_series_delta_theta   (dynsys system,
                                             anlsis analysis);

int draw_basin_of_attraction(double ref[][2], int ref_period,
                            dynsys system, anlsis analysis);

int draw_basin_of_attraction_clean	(double ref[][2], int ref_period,
                            		 dynsys system, anlsis analysis);

int draw_periodic_orbit_on_phase_space  (perorb po,
                                         dynsys system);

int draw_periodic_orbit_on_phase_space_clean(perorb po,
                                         	 dynsys system);

#endif