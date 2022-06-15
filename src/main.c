#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dynamical_system.h"
#include "spin_orbit.h"

int main(int argc, char **argv)
{
	/******************** Start clock **********************/

	clock_t begin_time = clock(), end_time;
	#ifdef _OMP_H
		double begin_time_omp = omp_get_wtime(), end_time_omp;
	#endif

	/***************** Declared variables ******************/

    double gamma;			// equatorial flattening
    double e;				// eccentricity
	double m_primary;		// mass of primary
	double m_secondary;		// mass of secondary
	double G;				// gravitational constant
	double a;				// semimajor axis
	double K;				// dissipation parameter

	double *params[7] = {&gamma,
						 &e,
						 &m_primary,
						 &m_secondary,
						 &G,
						 &a,
						 &K};

	dynsys system;
	dynsys system_two_body = init_two_body(*params);
	dynsys system_rigid = init_rigid(*params);
	dynsys system_rigid_kepler = init_rigid_kepler(*params);
	dynsys system_linear = init_linear(*params);
	dynsys system_linear_average = init_linear_average(*params);

	anlsis analysis;

	double orbital[4];

	/********************* Some values **********************/

	/** Moon
	 * e = 0.0549
	 * m_secondary = 1.215e-2
	**/

	/** Hyperion (Wisdom 1987)
	 * e = 0.1 
	 * gamma = (.89 * .89) / 3.
	 * m_secondary = 0. (assuming for now)
	**/
	  
	/////////////////////////////////////////////////////////
	/*				   		   Orbit		   	           */
	/////////////////////////////////////////////////////////

	// system = system_rigid;
	system = system_linear_average;
	double ic[system.dim];

	gamma = (.89 * .89) / 3.;
	e = 0.180; // e = 0.1;
	m_secondary = 0.;
	m_primary = 1.0 - m_secondary;
	G = 1.0;
	a = 1.0;
	K = 1e-2;

	analysis.number_of_cycles = 2e2; //1e3 6e3
	analysis.cycle_period = 2.0 * M_PI; // 1e-3

	// ic[0] = 0.0, ic[1] = 100.;
	// // //near the 1:1 stable fp in the rigid case
	// // ic[0] = M_PI; ic[1] = 0.551537;
	// init_orbital(orbital, e);
	// for (int i = 0; i < 4; i++) ic[i+2] = orbital[i];

	ic[0] = 0.0, ic[1] = 1000.;
	init_orbital(orbital, e);
	for (int i = 0; i < 4; i++) ic[i+2] = orbital[i];
	orbit_map(ic, system, analysis);

	// draw_orbit_map(system);

	// draw_orbit_on_phase_space(system);

	// multiple_time_series(system, analysis);

	// draw_multiple_time_series(system);

	// time_series(system, analysis);

	// draw_time_series(system);

	// for (e = 0.00; e < 0.205; e += 0.01)
	// {
	// 	time_series(system, analysis);
	// }

	// draw_time_series_union_e(system);

	// analysis.number_of_cycles = 1e3;
	// for (K = 0.01; K > 0.00095; K -= 0.001)
	// {
	// 	time_series(system, analysis);
	// }

	// draw_time_series_union_K(system);

	/////////////////////////////////////////////////////////
	/*				   		Phase space		   	           */
	/////////////////////////////////////////////////////////

	// system = system_rigid;
	// // system = system_linear_average;

	// gamma = (.89 * .89) / 3.;
	// e = 0.180; // e = 0.1;
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-2;

	// analysis.nc = 3, analysis.nv = 50; //nc = 3, nv = 50;
	// analysis.number_of_cycles = 2e3; //1e3
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.coordinate_min = 0.0; // M_PI
	// analysis.coordinate_max = M_PI; // M_PI 2.0* M_PI
	// analysis.velocity_min = 0.0;
	// analysis.velocity_max = 3.0;

	// phase_space(system, analysis);
	// draw_phase_space(system);

	// // for (e = 0.01; e < 0.205; e += 0.01)
	// // {
	// // 	// phase_space(system, analysis);
	// // 	draw_phase_space(system);
	// // }

	// // double gamma_hyperion = ((.89 * .89) / 3.); //~0.264
	// // double delta_gamma =  gamma_hyperion / 10.;

	// // for (gamma = 0.0; 
	// // 	 gamma < gamma_hyperion * 2.;
	// // 	 gamma += delta_gamma)
	// // {
	// // 	phase_space(system, analysis);
	// // 	draw_phase_space(system);
	// // }

	/////////////////////////////////////////////////////////
	/*				   		Mean resonance	   	           */
	/////////////////////////////////////////////////////////

// 	system = system_rigid;

// 	gamma = (.89 * .89) / 3.;
// 	e = 0.1;
// 	m_secondary = 0.0;
// 	m_primary = 1.0 - m_secondary;
// 	G = 1.0;
// 	a = 1.0;

// 	analysis.number_of_cycles = 1e2;
// 	analysis.cycle_period = 2.0 * M_PI;
// 	analysis.grid_resolution = 300;
// 	analysis.grid_coordinate_min = -M_PI; // 0.0
// 	analysis.grid_coordinate_max = M_PI; // 2.0*M_PI
// 	analysis.grid_velocity_min = 0.0;
// 	analysis.grid_velocity_max = 3.0;

// 	mean_resonance(system, analysis);

	/////////////////////////////////////////////////////////
	/*				Basin of attraction		   	           */
	/////////////////////////////////////////////////////////

	// system = system_linear_average;
	// double ref[2], ref_union[2][2];
	// int num_of_basins;

	// gamma = (.89 * .89) / 3.;
	// e = 0.1; //0.1
	// m_secondary = 0.0;
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;
	// a = 1.0;
 	// K = 1e-2;

	// analysis.number_of_cycles = 1e3; //1e3
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.grid_resolution = 600;
	// analysis.grid_coordinate_min = -M_PI;
	// analysis.grid_coordinate_max = M_PI;
	// analysis.grid_velocity_min = 0.0;
	// analysis.grid_velocity_max = 3.0;

  	// // ref[0] = 0.0; ref[1] = 0.551540; // e = 0.1
 	// // ref[0] = M_PI; ref[1] = 0.551540; // e = 0.1
  	// // ref[0] = 0.0; ref[1] = 1.0; // e = 0.0
	// // ref[0] = M_PI; ref[1] = 1.0; // e = 0.0
  	// // basin_of_attraction (ref, system, analysis);

	// num_of_basins = 2;
	// ref_union[0][0] = 0.0; ref_union[0][1] = 0.551540; // e = 0.1
	// ref_union[1][0] = M_PI; ref_union[1][1] = 0.551540; // e = 0.1
	// // ref_union[0][0] = 0.0; ref_union[0][1] = 1.0; // e = 0.0
	// // ref_union[1][0] = M_PI; ref_union[1][1] = 1.0; // e = 0.0
	// union_basin_of_attraction (ref_union, num_of_basins, system, analysis);

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