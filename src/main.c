#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dynamical_system.h"
#include "periodic_orbit.h"
#include "spin_orbit.h"

int main(int argc, char **argv)
{
	/******************** Start clock **********************/

	clock_t begin_time = clock(), end_time;
	#ifdef _OMP_H
		double begin_time_omp = omp_get_wtime(), end_time_omp;
	#endif

	/******************* Create output **********************/

	struct stat st = {0};
	if (stat("output", &st) == -1) {
		mkdir("output", 0700);
	}

	/*************** Some predefined values ****************/

	/** Moon
	 * e = 0.0549
	**/

	double e_moon = 0.0549;

	/** Hyperion (Wisdom 1987)
	 * e = 0.1 
	 * gamma = (.89 * .89) / 3.
	**/

	double e_hyperion = 0.1;
	double gamma_hyperion = ((.89 * .89) / 3.); //~0.264

	/****************** System parameters *******************/

    double gamma;			// equatorial flattening
    double e;				// eccentricity
	double m_primary;		// mass of primary
	double m_secondary;		// mass of secondary
	double G;				// gravitational constant
	double a;				// semimajor axis
	double K;				// dissipation parameter
	double T;				// system period

	gamma = (2.0/3.0) * 1e-4;
	e = 0.0549;
	m_secondary = 0.0;
	m_primary = 1.0 - m_secondary;
	G = 1.0;
	a = 1.0;
 	K = 1e-6;		//1e-4
	T = kepler_period(m_primary, m_secondary, G, a);

	// gamma = gamma_hyperion;
	// e = e_hyperion;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-2;
	// T = kepler_period(m_primary, m_secondary, G, a);

	/*************** System initialization ****************/

	dynsys	system;
	int		number_of_params = 8;
	double	params[8] = {gamma, e, m_primary, m_secondary,
						 G, a, K, T};

	dynsys	system_two_body = init_two_body(params);
	dynsys	system_rigid = init_rigid(params);
	dynsys	system_rigid_kepler = init_rigid_kepler(params);
	dynsys	system_linear = init_linear(params);
	dynsys	system_linear_average = init_linear_average(params);

	// system = system_rigid;
	system = system_linear;
	// system = system_linear_average;

	/**************** Simulation parameters ******************/

	anlsis	analysis;

	analysis.number_of_cycles = (int) (1.0e7 / T);	// (int) (1.0e7 / T)
	// analysis.number_of_cycles = 6e3;	// 1e3 5e3 1e4 6e3 (hyp - strong)
	analysis.cycle_period = T;
	analysis.evolve_box_size = 1e8;

	analysis.nc = 3;					// 3
	analysis.nv = 50;					// 50;
	analysis.coordinate_min = 0.0; 		// 0.0
	analysis.coordinate_max = M_PI; 	// M_PI
	analysis.velocity_min = 0.8;		// 0.0
	analysis.velocity_max = 1.2;		// 3.0
	// analysis.velocity_max = 3.0;		// 3.0

	analysis.grid_resolution = 300;			// 600
	analysis.grid_coordinate_min = -M_PI;	// -M_PI
	analysis.grid_coordinate_max = M_PI;	// M_PI
	analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 5.0;		// 3.0
	analysis.grid_velocity_max = 3.0;		// 3.0

	analysis.sqrt_orbits_on_box = 10;
	
	analysis.spin_period_min = 1;
	analysis.orbit_period_min = 1;
	analysis.spin_period_max = 5;			// 5
	analysis.orbit_period_max = 4;			// 4
	analysis.evolve_basin_time_tol = 500;
	analysis.evolve_basin_eps = 1e-1;

	analysis.po_max_step = 1000;			// 1000
	analysis.po_tol = 1e-8;					// 1e-13

	analysis.number_of_rand_orbits_mc = 10000;
	analysis.convergence_window_mc = 1e3;		// 1e3 (moon) 5e2 (hyp - weak) 1e3 (hyp - strong)
	analysis.convergence_precision_mc = 1e-2;

	analysis.convergence_transient_wn = 1e4; 	// 1e4 (moon) 1e3 (hyp - weak) 4e3 (hyp - strong)
	analysis.convergence_window_wn = 1e3;	 	// 1e3 (moon) 5e2 (hyp - weak) 1e3 (hyp - strong)
	analysis.convergence_precision_wn = 1e-2; 	// 1e-2

	/*********** Numerical integrator parameters *************/

	rngkta	rk;

	rk.method = "rk8pd";
	rk.control = "fixed";
	rk.n = 100;
	rk.h = T * (1.0 / (double) rk.n) * sign(analysis.cycle_period);
	rk.error_abs = 1e-8;
	rk.error_rel = 1e-8;

	// rk.method = "rk8pd";
	// rk.control = "adaptive";
	// rk.h = 1e-3 * sign(analysis.cycle_period);
	// rk.error_abs = 1e-14;
	// rk.error_rel = 0.0;
	// rk.h_max = 1e-1;
	// rk.h_min = 1e-11;

	/***************** Declared variables *******************/

	perorb	po;
	perorb 	*multiple_pos;

	int 	number_of_pos;
	int 	number_of_e;

	double	ic[system.dim];
	double 	e_initial;
	double 	e_final;
	double 	e_step;

	/////////////////////////////////////////////////////////
	/*				   		   Orbit		   	           */
	/////////////////////////////////////////////////////////

	// time_t t;
	// srand((unsigned) time(&t));
	// ic[0] = rand_number_in_interval(analysis.grid_coordinate_min, analysis.grid_coordinate_max);
	// ic[1] = rand_number_in_interval(analysis.grid_velocity_min, analysis.grid_velocity_max);
	// complete_orbital_part(ic, system);
	// orbit_map(ic, system, analysis, rk);

	// init_orbital(ic, system);
	// orbit_two_body(ic, system, analysis);

	// draw_orbit_map(system);
	// draw_orbit_on_phase_space(system);
	// draw_orbit_on_phase_space_latex(system);

	/////////////////////////////////////////////////////////
	/*				   	Periodic Orbit		   	           */
	///////////////////////////////	//////////////////////////

	// po.period = 4;

	// alloc_2d_double(&po.orbit, po.period, system.dim);
	
	// po.seed[0] = -1.630068286090758e+00;
	// po.seed[1] =  1.249822122223460e+00;
	// periodic_orbit(&po, system, analysis);

	// ic[0] = po.initial_condition[0];
	// ic[1] = po.initial_condition[1];
	// complete_orbital_part(ic, system);
	// orbit_map(ic, system, analysis);

	// dealloc_2d_double(&po.orbit, po.period);

	// draw_periodic_orbit_on_phase_space (po, system);
	// draw_periodic_orbit_on_phase_space_clean (po, system);

	/////////////////////////////////////////////////////////
	/*				   		Time series		   	           */
	/////////////////////////////////////////////////////////

	// analysis.number_of_cycles = 4e3; //1e3 6e3
	// analysis.cycle_period = 2.0 * M_PI; // 1e-3
	// analysis.evolve_box_size = 1e8;
	// analysis.evolve_basin_eps = 1e-1;
	// analysis.number_of_time_series = 20;
	// analysis.time_series_delta = 10;

	// multiple_time_series(system, analysis);

	// draw_multiple_time_series(system);

	// multiple_time_series_delta_theta_dot(system, analysis);
	// draw_multiple_time_series_delta_theta_dot(system, analysis);
	// draw_multiple_time_series_delta_theta_dot_latex(system, analysis);

	// // multiple_time_series_delta_theta(system, analysis);
	// draw_multiple_time_series_delta_theta(system);

	// time_series(system, analysis);

	// // draw_time_series(system);
	// draw_time_series_latex(system);

	// analysis.number_of_cycles = 2e2;
	// for (e = 0.00; e < 0.205; e += 0.01)
	// {
	// 	time_series(system, analysis);
	// }

	// draw_time_series_union_e(system);
	// draw_time_series_union_e_latex(system);
	// draw_time_series_union_e_eps(system);

	// analysis.number_of_cycles = 1e3;
	// for (K = 0.01; K > 0.00095; K -= 0.001)
	// {
	// 	time_series(system, analysis);
	// }

	// draw_time_series_union_K(system);
	// draw_time_series_union_K_latex(system);

	/////////////////////////////////////////////////////////
	/*				   		Phase space		   	           */
	/////////////////////////////////////////////////////////

	// phase_space(system, analysis, rk);
	// draw_phase_space(system, analysis);
	// draw_phase_space_latex(system);
	// draw_phase_space_clean(system);

	// e_initial = 0.0;
	// e_final = 0.2;
	// e_step = 0.01;

	// for (double e_loop = e_initial; e_loop < e_final; e_loop += e_step)
	// {
	// 	dynsys system_loop = system;
	// 	double params_loop[number_of_params];
	// 	copy(params_loop, params, number_of_params);
	// 	params_loop[1] = e_loop;
	// 	system_loop.params = params_loop;
	// 	phase_space(system_loop, analysis);
	// 	draw_phase_space(system_loop, analysis);
	// }

	/////////////////////////////////////////////////////////
	/*				Basin of attraction		   	           */
	/////////////////////////////////////////////////////////

	/* One periodic orbit */

	// // po.period = 1;
	// // po.seed[0] = 0.0; po.seed[1] = 0.551537; // e = 0.1 SFP 1/1 resonance
	// po.period = 2;
	// po.seed[0] = -1.50359; po.seed[1] = 0.860586; // e = 0.1 period 2 SPO 1/2 system linear
	// alloc_2d_double(&po.orbit, po.period, system.dim);
	// periodic_orbit(&po, system, analysis);
  	// basin_of_attraction (po, system, analysis);
	// draw_basin_of_attraction (po, system, analysis);
	// // draw_basin_of_attraction_clean (po.period, po.initial_condition, system, analysis);
	// dealloc_2d_double(&po.orbit, po.period);

	/* Multiple periodic orbits */

	anlsis analysis_res_300 = analysis;
	analysis_res_300.grid_resolution = 300;

	fill_attractor_array(&number_of_pos, &multiple_pos, system, analysis_res_300, rk);

	multiple_basin_of_attraction_determined (number_of_pos, multiple_pos, system, analysis, rk);
	basin_size_from_data (number_of_pos, multiple_pos, system, analysis);
	draw_multiple_basin_of_attraction_determined (system, analysis);

	// if (number_of_pos > 0)
	// {
	// 	for (int j = 0; j < number_of_pos; j++)
	// 	{
	// 		dealloc_2d_double(&multiple_pos[j].orbit, multiple_pos[j].period);
	// 	}
	// 	free(multiple_pos);
	// }

	// multiple_basin_of_attraction_undetermined_monte_carlo_with_break(system, analysis);

	/* Multiple periodic orbits - loop over e */

	// number_of_e = 10;	// 50 (moon) 20 (hyp)
	// e_initial = 0.0;	// 0.0 (moon) 0.0 (hyp)
	// e_final = 0.1;		// 0.5 (moon) 0.2 (hyp)

	// e_step = (e_final - e_initial) / (double)(number_of_e);

	// for(int i = 0; i <= number_of_e; i++)	// (int i = 0; i <= number_of_e; i++)
	// {
	// 	double e_loop = e_initial + (double)i * e_step;
	// 	printf("e = %1.3f\n", e_loop);

	// 	dynsys system_loop = system;
	// 	double params_loop[number_of_params];
	// 	copy(params_loop, params, number_of_params);
	// 	params_loop[1] = e_loop;
	// 	system_loop.params = params_loop;

	// 	anlsis analysis_res_300 = analysis;
	// 	analysis_res_300.grid_resolution = 300;

	// 	fill_attractor_array(&number_of_pos, &multiple_pos, system_loop, analysis_res_300, rk);

	// 	if (number_of_pos > 0)
	// 	{
			/* calculation on grid */
	
			// multiple_basin_of_attraction_determined (number_of_pos, multiple_pos, system_loop, analysis, rk);
			// draw_multiple_basin_of_attraction_determined (system_loop, analysis);
			// draw_multiple_basin_of_attraction_determined_clean (system_loop, analysis);
			// basin_size_from_data (number_of_pos, multiple_pos, system_loop, analysis);
			// basin_entropy_from_data (system_loop, analysis);
			// basin_entropy_progress_from_data (number_of_pos, multiple_pos, system_loop, analysis);
			// basin_entropy_vs_box_size (number_of_pos, multiple_pos, system_loop, analysis);
			// plot_basin_entropy_vs_box_size (system_loop, analysis);

			/* calculation monte carlo */

			// multiple_basin_of_attraction_determined_monte_carlo (number_of_pos, multiple_pos, system_loop, analysis);
			// basin_size_from_data_monte_carlo (number_of_pos, multiple_pos, system_loop, analysis);
			// basin_entropy_from_data_monte_carlo (system_loop, analysis);
			// basin_entropy_progress_from_data_monte_carlo (number_of_pos, multiple_pos, system_loop, analysis);
			// multiple_basin_of_attraction_determined_monte_carlo_with_break (number_of_pos, multiple_pos, system_loop, analysis);
			// basin_size_from_data_monte_carlo_with_break (number_of_pos, multiple_pos, system_loop, analysis);
			// basin_entropy_from_data_monte_carlo_with_break (system_loop, analysis);
			// basin_entropy_progress_from_data_monte_carlo_with_break (number_of_pos, multiple_pos, system_loop, analysis);
			// multiple_basin_of_attraction_undetermined_monte_carlo_with_break(system_loop, analysis);
			
			/* comparison between grid and monte carlo */

			// comparison_entropy_grid_vs_monte_carlo (number_of_pos, multiple_pos, system_loop, analysis);
			// plot_comparison_entropy_grid_vs_monte_carlo (system_loop, analysis);

		// 	for (int j = 0; j < number_of_pos; j++)
		// 	{
		// 		dealloc_2d_double(&multiple_pos[j].orbit, multiple_pos[j].period);
		// 	}
		// 	free(multiple_pos);
		// }
		// else
		// {
		// 	printf("Warning: null number of attractors.\n");
		// }

		// plot_histogram_python (system_loop, analysis);
		// plot_histogram_python_monte_carlo_with_break (system_loop, analysis);

	// }

	// plot_size_multiple_basin_of_attraction_determined_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// analysis.grid_resolution = 600;
	// plot_size_multiple_basin_of_attraction_determined_range_e_latex(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// plot_size_multiple_basin_of_attraction_determined_plus_basin_entropy_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// plot_size_multiple_basin_of_attraction_determined_plus_basin_entropy_monte_carlo_with_break_range_e(0, 0.01,
	// 	0.0, 0.25, system, analysis);

	// plot_size_multiple_basin_of_attraction_undetermined_plus_basin_entropy_monte_carlo_with_break_range_e(0, 0.01,
	// 	e_initial, e_final, system, analysis);

	// analysis.grid_resolution = 600;
	// plot_slope_basin_entropy_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// analysis.grid_resolution = 600;
	// plot_basin_entropy_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// plot_entropy_comparison_monte_carlo_range_e(number_of_e, e_initial, e_final, system, analysis);

	// plot_entropy_comparison_monte_carlo_v2_range_e(number_of_e, e_initial, e_final, system, analysis);

	/////////////////////////////////////////////////////////
	/*						Benchmark		   	           */
	/////////////////////////////////////////////////////////

	// linear_average_benchmark();

	// trace_ellipse();

	/////////////////////////////////////////////////////////
	/*				   		 Testing		   	           */
	/////////////////////////////////////////////////////////

	/******************** Stop clock ***********************/

	end_time = clock();
	double time_spent 
		= (double)(end_time - begin_time) / CLOCKS_PER_SEC;
	if (time_spent < 60.0)
	{
		printf("time spent = %1.2e seconds\n", 
				time_spent);
	}
	else if (time_spent < 3600.0)
	{
		printf("time spent = %1.2e minutes\n", 
				time_spent/60.0);
	}
	else
	{
		printf("time spent = %1.2e hours\n", 
				time_spent/3600.0);
	}

	#ifdef _OMP_H
		end_time_omp = omp_get_wtime();
		double time_spent_omp
			= (double)(end_time_omp - begin_time_omp);
		if (time_spent_omp < 60.0)
		{
			printf("time_spent_omp = %1.2e seconds\n", 
					time_spent_omp);
		}
		else if (time_spent_omp < 3600.0)
		{
			printf("time_spent_omp = %1.2e minutes\n", 
					time_spent_omp/60.0);
		}
		else
		{
			printf("time_spent_omp = %1.2e hours\n", 
					time_spent_omp/3600.0);
		}
	#endif

	/*******************************************************/

	return 0;
}