#ifndef SPIN_ORB_H
#define SPIN_ORB_H

#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

double vis_viva_two_body(double y[4],
                         dynsys system);

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

int complete_orbital_part   (double y[],
                             dynsys system);

/**
 * poincare map
**/

int orbit_map (double *ic,
               dynsys system,
               anlsis analysis);

int phase_space (dynsys system,
                anlsis analysis);

/**
 * time series
**/

int time_series(dynsys system,
                anlsis analysis);

int multiple_time_series(dynsys system,
                        anlsis analysis);

int multiple_time_series_delta_theta_dot(dynsys system,
										anlsis analysis);

int multiple_time_series_delta_theta(dynsys system,
									anlsis analysis);

/**
 * periodic orbit
**/

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

// search the phase space for number_of_candidates 
// resonances of type spin_period / orbit_period
int look_for_resonance	(int number_of_candidates,
						 double candidates[][2],
						 int spin_period,
                         int orbit_period,
                         dynsys system, 
						 anlsis analysis);

// uses look_for_resonance to find all the stable spin-orbit 
// resonances on the phase space that have spin and orbit 
// periods inside the range given in analysis
// (can be easily modified to find all upos as well)
int find_all_periodic_attractors(int *number_of_pos,
							 	 perorb **multiple_pos,
							 	 dynsys system,
                         	 	 anlsis analysis);

// fills in the array of attractors multiple_pos[] looking for 
// written files or using find_all_periodic_attractors()
int fill_attractor_array(int *number_of_pos,
						 perorb **multiple_pos,
                         dynsys system,
                         anlsis analysis);

// retrieves the information of the basin of attraction
// data file and produces a matrix of cantor values
// representing the basins of SORs
int fill_basin_matrix	(double	***basin_matrix,
						 dynsys system,
                         anlsis analysis);

// retrieves the information of the basin of attraction
// data file and produces a matrix of po ids
// representing the basins of SORs
int fill_control_matrix	(int	***control_matrix,
						 dynsys system,
                         anlsis analysis);

// retrieves the information of the basin of attraction
// data file for Monte Carlo and produces an array of po
// ids representing the basins of SORs
int fill_control_monte_carlo(int	**control,
						 	 dynsys system,
                         	 anlsis analysis);

int fill_control_monte_carlo_with_break(int *control_size,
										int	**control,
						 	 			dynsys system,
                         	 			anlsis analysis);

/**
 * basin of attraction
**/

double dist_from_ref(double x[2],
                     double ref[2]);

int evolve_basin(double *ic,
				 bool *converged,
                 int *convergence_time,
                 perorb po,
                 dynsys system,
				 anlsis analysis);

int basin_of_attraction (perorb po,
                         dynsys system,
                         anlsis analysis);

int evolve_multiple_basin_determined(double *ic,
									 int number_of_po,
									 int *converged,
									 int *convergence_time,
									 perorb po[],
									 dynsys system,
									 anlsis analysis);

int multiple_basin_of_attraction_determined (int number_of_po,
											 perorb po[],
                         					 dynsys system,
                         					 anlsis analysis);

int multiple_basin_of_attraction_determined_monte_carlo (int number_of_po,
											 			 perorb po[],
                         					 			 dynsys system,
                         					 			 anlsis analysis);

int multiple_basin_of_attraction_determined_monte_carlo_with_break	(int number_of_po,
											 			 			 perorb po[],
                         					 			 			 dynsys system,
                         					 			 			 anlsis analysis);

int comparison_entropy_grid_vs_monte_carlo	(int number_of_po,
											 perorb po[],
                         					 dynsys system,
                         					 anlsis analysis);

int basin_entropy_vs_box_size	(int number_of_po,
								 perorb po[],
                         		 dynsys system,
                         		 anlsis analysis);

int basin_size_from_data(int number_of_po,
						 perorb po[],
                         dynsys system,
                         anlsis analysis);

int basin_size_from_data_monte_carlo(int number_of_po,
						 			 perorb po[],
                         			 dynsys system,
                         			 anlsis analysis);

int basin_size_from_data_monte_carlo_with_break(int number_of_po,
						 			 			perorb po[],
                         			 			dynsys system,
                         			 			anlsis analysis);

int basin_entropy_from_data (dynsys system,
                         	 anlsis analysis);

int basin_entropy_from_data_monte_carlo (dynsys system,
                         	 			 anlsis analysis);

int basin_entropy_from_data_monte_carlo_with_break	(dynsys system,
                         	 			 			 anlsis analysis);

int basin_entropy_progress_from_data(int number_of_po,
						 			 perorb po[],
                         			 dynsys system,
                         			 anlsis analysis);

int basin_entropy_progress_from_data_monte_carlo(int number_of_po,
						 			 			 perorb po[],
                         			 			 dynsys system,
                         			 			 anlsis analysis);

int basin_entropy_progress_from_data_monte_carlo_with_break(int number_of_po,
						 			 			 			perorb po[],
                         			 			 			dynsys system,
                         			 			 			anlsis analysis);

int evolve_multiple_basin_undetermined  (double *ic,
									     bool *converged,
                                         int *attractor_period,
									     int *convergence_time,
									     dynsys system,
									     anlsis analysis);

int multiple_basin_of_attraction_undetermined   (dynsys system,
                         					     anlsis analysis);

double basin_entropy(int number_of_orbits,
					 int number_of_po,
                     int *basin_size,
                     perorb po[],
					 anlsis analysis);

/**
 * benchmark tests
**/

// comparison between field_linear and field_linar_average
int linear_average_benchmark();

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

int draw_periodic_orbit_on_phase_space  (perorb po,
                                         dynsys system);

int draw_periodic_orbit_on_phase_space_clean(perorb po,
                                         	 dynsys system);

int draw_basin_of_attraction(perorb po,
                             dynsys system,
                             anlsis analysis);

int draw_basin_of_attraction_clean	(int ref_period, double ref[][2],
                            		 dynsys system, anlsis analysis);

int draw_multiple_basin_of_attraction_determined(dynsys system,
                                        		 anlsis analysis);

int draw_multiple_basin_of_attraction_determined_clean  (dynsys system,
                                        		         anlsis analysis);

int plot_size_multiple_basin_of_attraction_determined_range_e	(int number_of_e,
																 double e_initial,
																 double e_final,
																 dynsys system,
																 anlsis analysis);

int plot_size_multiple_basin_of_attraction_determined_range_e_latex	(int number_of_e,
																 	 double e_initial,
																 	 double e_final,
																 	 dynsys system,
																 	 anlsis analysis);

int draw_multiple_basin_of_attraction_undetermined  (dynsys system,
                                        		     anlsis analysis);

int plot_basin_entropy_vs_box_size	(dynsys system,
                         			 anlsis analysis);

int plot_slope_basin_entropy_range_e(int number_of_e,
									 double e_initial,
									 double e_final,
									 dynsys system,
									 anlsis analysis);

int plot_basin_entropy_range_e	(int number_of_e,
								 double e_initial,
								 double e_final,
								 dynsys system,
								 anlsis analysis);

int plot_size_multiple_basin_of_attraction_determined_plus_basin_entropy_range_e(int number_of_e,
																 	 			 double e_initial,
																 	 			 double e_final,
																 	 			 dynsys system,
																 	 			 anlsis analysis);

int plot_size_multiple_basin_of_attraction_determined_plus_basin_entropy_monte_carlo_with_break_range_e(int number_of_e,
																 	 			 						double e_initial,
																 	 			 						double e_final,
																 	 			 						dynsys system,
																 	 			 						anlsis analysis);

int plot_entropy_comparison_monte_carlo_range_e	(int number_of_e,
												 double e_initial,
												 double e_final,
												 dynsys system,
												 anlsis analysis);

int plot_comparison_entropy_grid_vs_monte_carlo	(dynsys system,
												 anlsis analysis);

/**
 * Python pipe
**/

int plot_histogram_python	(dynsys system,
							 anlsis analysis);
                                                                                   
#endif