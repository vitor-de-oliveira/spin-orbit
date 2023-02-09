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

	/***************** Declared variables ******************/

    double gamma;			// equatorial flattening
    double e;				// eccentricity
	double m_primary;		// mass of primary
	double m_secondary;		// mass of secondary
	double G;				// gravitational constant
	double a;				// semimajor axis
	double K;				// dissipation parameter
	double T;				// system period

	double *params[8] = {&gamma,
						 &e,
						 &m_primary,
						 &m_secondary,
						 &G,
						 &a,
						 &K,
						 &T};

	dynsys system;
	dynsys system_two_body = init_two_body(*params);
	dynsys system_rigid = init_rigid(*params);
	dynsys system_rigid_kepler = init_rigid_kepler(*params);
	dynsys system_linear = init_linear(*params);
	dynsys system_linear_average = init_linear_average(*params);

	anlsis analysis;

	perorb po;

	double orbital[4];

	/********************* Some values **********************/

	/** Moon
	 * e = 0.0549
	 * m_secondary = 1.215e-2
	**/

	double e_moon = 0.0549;

	/** Hyperion (Wisdom 1987)
	 * e = 0.1 
	 * gamma = (.89 * .89) / 3.
	 * m_secondary = 0. (assuming for now)
	**/

	double e_hyperion = 0.1;
	double gamma_hyperion = ((.89 * .89) / 3.); //~0.264
	  
	/////////////////////////////////////////////////////////
	/*				   		   Orbit		   	           */
	/////////////////////////////////////////////////////////

	// system = system_rigid;
	// system = system_linear_average;
	// system = system_linear;
	// double ic[system.dim];

	// gamma = gamma_hyperion;
	// e = 0.105;
	// m_secondary = 0.;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
	// K = 1e-2;
	// T = 2.0 * M_PI;

	// analysis.number_of_cycles = 5e3;		//1e3 6e3
	// analysis.cycle_period = T; 				// 1e-3
	// analysis.evolve_box_size = 1e8;

	// ic[0] = -6.696425553095295e-01 + M_PI;
	// ic[1] = 2.102460925123672e+00;
	// init_orbital(orbital, e);
	// for (int i = 0; i < 4; i++) ic[i+2] = orbital[i];
	// orbit_map(ic, system, analysis);

	/////////////////////////////////////////////////////////
	/*				   	Periodic Orbit		   	           */
	///////////////////////////////	//////////////////////////

	// // system = system_rigid;
	// // system = system_linear_average;
	// system = system_linear;

	// gamma = gamma_hyperion;
	// e = 0.02;
	// m_secondary = 0.;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
	// K = 1e-2;
	// T = 2.0 * M_PI;

	// analysis.cycle_period = T;
	// analysis.evolve_box_size = 1e8;

	// analysis.po_max_step = 1000;			// 1000
	// analysis.po_tol = 1e-8;					// 1e-13

	// po.period = 2;

	// alloc_2d_double(&po.orbit, po.period, system.dim);
	
	// // po.seed[0] = 0.0; po.seed[1] = 0.551537; // e = 0.1 SFP 1/1 resonance
	// // po.seed[0] = 0.0; po.seed[1] = 2.32185; // e = 0.1 SFP 2/1 resonance  
	// // po.seed[0] = -1.56892; po.seed[1] = 0.868688; // e = 0.1 period 2 SPO 1/2 resonance 
	// // po.seed[0] = 0.0; po.seed[1] = 1.87878; // e = 0.1 period 2 UPO 2/2 resonance
	// // po.seed[0] = -1.57310; po.seed[1] =  1.71059; // e = 0.1 UFP 2/1 resonance 
	// // po.seed[0] = -1.57246; po.seed[1] =  2.14877; // e = 0.1 UPO 5/2 resonance
	// // po.seed[0] = -1.94124; po.seed[1] =  1.46147; // e = 0.1 period 2 UPO resonance 3/2
	// // po.seed[0] = 1.35558; po.seed[1] =  1.08285; // e = 0.1 UPO 2/2
	// // po.seed[0] = 1.94124; po.seed[1] =  1.46147; // e = 0.1 period 2 UPO resonance 3/2
	// // po.seed[0] = 0.0; po.seed[1] =  2.72177; // e = 0.1 SPO 5/2
	// // po.seed[0] = -1.57079; po.seed[1] =  1.95929; // e = 0.1 SPO 9/4
	// // po.seed[0] = -0.234070; po.seed[1] = 2.29795; // e = 0.1 SFP 2/1 system linear
	// po.seed[0] = -1.50359; po.seed[1] = 0.860586; // e = 0.1 period 2 SPO 1/2 system linear
	// // po.seed[0] = -0.0257629; po.seed[1] = 0.484803; // e = 0.140 period 3 UPO around 1/1 resonance
	// // po.seed[0] = 0.7; po.seed[1] = 1.1; // e = 0.2 ?
	// // po.seed[0] = 0.0; po.seed[1] = 0.380929; // e = 0.2 UFP 1/1 resonance
	// // po.seed[0] = 0.0; po.seed[1] = 0.722967; // e = 0.2 SFP 2/2 resonance
	// po.seed[0] = -2.335247192600534e-01;
	// po.seed[1] = 1.906315169310492;
	// periodic_orbit(&po, system, analysis);

	// // // draw_periodic_orbit_on_phase_space (po, system);
	// // // // draw_periodic_orbit_on_phase_space_clean (po, system);

	// dealloc_2d_double(&po.orbit, po.period);

	// analysis.number_of_cycles = 1e2;
	// // analysis.cycle_period = 1e-3;

	// double ic_po[system.dim];
	// ic_po[0] = po.initial_condition[0];
	// ic_po[1] = po.initial_condition[1];
	// init_orbital(orbital, e);
	// for (int i = 0; i < 4; i++) ic_po[i+2] = orbital[i];
	// orbit_map(ic_po, system, analysis);

	// draw_orbit_map(system);

	// draw_orbit_on_phase_space(system);
	// draw_orbit_on_phase_space_latex(system);

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

	// system = system_rigid;
	// // system = system_linear_average;

	// gamma = gamma_hyperion;
	// e = 0.105;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-2;
	// T = 2.0 * M_PI;

	// analysis.nc = 3;						// 3
	// analysis.nv = 50;						// 50;
	// analysis.number_of_cycles = 2e3; 		// 1e3
	// analysis.cycle_period = T;
	// analysis.evolve_box_size = 1e8;
	// analysis.coordinate_min = 0.0; 	// 0.0
	// analysis.coordinate_max = M_PI; 	// M_PI
	// analysis.velocity_min = 0.0;		// 0.0
	// analysis.velocity_max = 3.0;		// 3.0

	// phase_space(system, analysis);
	// draw_phase_space(system);
	// draw_phase_space_latex(system);

	// gamma = gamma_hyperion;

	// for (e = 0.0; e < 0.405; e += 0.01)
	// {
	// 	// phase_space(system, analysis);
	// 	// draw_phase_space(system);
	// 	draw_phase_space_clean(system);
	// }

	// e = e_hyperion;

	// for (gamma = 0.0; 
	// 	 gamma < gamma_hyperion * 2. + gamma_hyperion / 40.;
	// 	 gamma += gamma_hyperion / 20.)
	// {
	// 	// phase_space(system, analysis);
	// 	// draw_phase_space(system);
	// 	draw_phase_space_clean(system);
	// }

	/////////////////////////////////////////////////////////
	/*				Basin of attraction		   	           */
	/////////////////////////////////////////////////////////

	// system = system_linear_average;
	system = system_linear;

	gamma = gamma_hyperion;
	e = 0.14;
	m_secondary = 0.0;
	m_primary = 1.0 - m_secondary;
	G = 1.0;
	a = 1.0;
 	K = 1e-2;
	T = 2.0 * M_PI;

	// analysis.number_of_cycles = 1e3; //1e3
	// analysis.cycle_period = T;
	// analysis.evolve_box_size = 1e8;

	// analysis.grid_resolution = 100;
	// analysis.grid_coordinate_min = -M_PI;
	// analysis.grid_coordinate_max = M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;
	
	// analysis.evolve_basin_time_tol = 100;
	// analysis.evolve_basin_eps = 1e-1;

	// // // po.period = 1;
	// // // po.seed[0] = 0.0; po.seed[1] = 0.551537; // e = 0.1 SFP 1/1 resonance
	// // po.period = 2;
	// // po.seed[0] = -1.50359; po.seed[1] = 0.860586; // e = 0.1 period 2 SPO 1/2 system linear
	// // alloc_2d_double(&po.orbit, po.period, system.dim);
	// // periodic_orbit(&po, system, analysis);
  	// // basin_of_attraction (po, system, analysis);
	// // draw_basin_of_attraction (po, system, analysis);
	// // // draw_basin_of_attraction_clean (po.period, po.initial_condition, system, analysis);
	// // dealloc_2d_double(&po.orbit, po.period);

	// int number_of_po = 2;
	// perorb multiple_po[number_of_po];

	// multiple_po[0].period = 1;
	// multiple_po[0].seed[0] = 0.0; multiple_po[0].seed[1] = 0.551537; // e = 0.1 SFP 1/1 resonance
	
	// multiple_po[1].period = 2;
	// multiple_po[1].seed[0] = -1.50359; multiple_po[1].seed[1] = 0.860586; // e = 0.1 period 2 SPO 1/2 system linear
	
	// // multiple_po[2].period = 1;
	// // multiple_po[2].seed[0] = M_PI; multiple_po[2].seed[1] = 0.551537;

	// for (int i = 0; i < number_of_po; i++)
	// {
	// 	alloc_2d_double(&multiple_po[i].orbit, multiple_po[i].period, system.dim);
	// 	periodic_orbit(&multiple_po[i], system, analysis);
	// }

	// multiple_basin_of_attraction_determined (number_of_po, multiple_po, system, analysis);
	// draw_multiple_basin_of_attraction_determined (number_of_po, system, analysis);
	
	// for (int i = 0; i < number_of_po; i++)
	// {
	// 	dealloc_2d_double(&multiple_po[i].orbit, multiple_po[i].period);
	// }

	// analysis.number_of_cycles = 1e3; //1e3
	// analysis.cycle_period = T;
	// analysis.evolve_box_size = 1e8;

	// analysis.grid_resolution = 600;
	// analysis.grid_coordinate_min = -M_PI;
	// analysis.grid_coordinate_max = M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;
	
	// analysis.evolve_basin_time_tol = 100;

	// analysis.po_max_step = 1000;
	// analysis.po_tol = 1e-10;	

	// multiple_basin_of_attraction_undetermined (system, analysis);
	// draw_multiple_basin_of_attraction_undetermined (system, analysis);

	int 	number_of_pos;
	perorb 	*multiple_pos;
	int 	number_of_e;
	double 	e_initial;
	double 	e_final;
	double 	e_step;

	analysis.number_of_cycles = 3e3;		// 1e3
	analysis.cycle_period = T;
	analysis.evolve_box_size = 1e8;

	analysis.grid_resolution = 600;
	analysis.grid_coordinate_min = -M_PI;	// -M_PI
	analysis.grid_coordinate_max = M_PI;	// M_PI
	analysis.grid_velocity_min = 0.0;
	analysis.grid_velocity_max = 3.0;

	analysis.sqrt_orbits_on_box = 10;
	
	analysis.spin_period_min = 1;
	analysis.orbit_period_min = 1;
	analysis.spin_period_max = 5;
	analysis.orbit_period_max = 4;
	analysis.evolve_basin_time_tol = 100;
	analysis.evolve_basin_eps = 1e-1;

	analysis.po_max_step = 1000;			// 1000
	analysis.po_tol = 1e-8;					// 1e-13

	analysis.number_of_rand_orbits = 360000;
	analysis.convergence_window = 1e4;
	analysis.convergence_precision = 1e-2;

	number_of_e = 20;
	e_initial = 0.0;
	e_final = 0.2;

	e_step = (e_final - e_initial) / (double)(number_of_e);

	// for (int i = 0; i <= number_of_e; i++)
	// {
	// 	e = e_initial + (double)i * e_step;
	// 	printf("e = %1.3f\n", e);

		// fill_attractor_array(&number_of_pos, &multiple_pos, system, analysis);

		// if (number_of_pos > 0)
		// {
			/* calculation on grid */
	
			// multiple_basin_of_attraction_determined (number_of_pos, multiple_pos, system, analysis);
			// draw_multiple_basin_of_attraction_determined (system, analysis);
			// basin_size_from_data (number_of_pos, multiple_pos, system, analysis);
			// basin_entropy_from_data (system, analysis);
			// basin_entropy_progress_from_data (number_of_pos, multiple_pos, system, analysis);
			// basin_entropy_vs_box_size (number_of_pos, multiple_pos, system, analysis);
			// plot_basin_entropy_vs_box_size (system, analysis);

			/* calculation monte carlo */

			// multiple_basin_of_attraction_determined_monte_carlo (number_of_pos, multiple_pos, system, analysis);
			// basin_size_from_data_monte_carlo (number_of_pos, multiple_pos, system, analysis);
			// basin_entropy_from_data_monte_carlo (system, analysis);
			// basin_entropy_progress_from_data_monte_carlo (number_of_pos, multiple_pos, system, analysis);

			/* comparison between grid and monte carlo */

			// comparison_entropy_grid_vs_monte_carlo (number_of_pos, multiple_pos, system, analysis);
			// plot_comparison_entropy_grid_vs_monte_carlo (system, analysis);

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

		// plot_histogram_python (system, analysis);

	// }

	// plot_size_multiple_basin_of_attraction_determined_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// analysis.grid_resolution = 600;
	// plot_size_multiple_basin_of_attraction_determined_range_e_latex(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// plot_size_multiple_basin_of_attraction_determined_plus_basin_entropy_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// analysis.grid_resolution = 600;
	// plot_slope_basin_entropy_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// analysis.grid_resolution = 600;
	// plot_basin_entropy_range_e(number_of_e,
	// 	e_initial, e_final, system, analysis);

	// for (e = 0.0; e < 0.21; e += 0.01)
	// {
	// 	draw_multiple_basin_of_attraction_determined_clean (system, analysis);
	// }

	// plot_entropy_comparison_monte_carlo_range_e(number_of_e, e_initial, e_final, system, analysis);

	// plot_entropy_comparison_monte_carlo_v2_range_e(number_of_e, e_initial, e_final, system, analysis);

	/////////////////////////////////////////////////////////
	/*						Benchmark		   	           */
	/////////////////////////////////////////////////////////

	// linear_average_benchmark();

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