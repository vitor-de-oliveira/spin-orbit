#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main_functions.h"
#include "dynamical_system.h"

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

	double *params[5] 
		= {&gamma, &e, &m_primary, &m_secondary, &G};

	dynsys system;
	dynsys system_rigid = init_rigid(*params);
	dynsys system_rigid_kepler = init_rigid_kepler(*params);
	dynsys system_two_body = init_two_body(*params);

	anlsis analysis;

	double orbital[4];

	/*******************************************************/
	
	/////////////////////////////////////////////////////////
	/*				   		   Orbit		   	           */
	/////////////////////////////////////////////////////////

	system = system_rigid;
	double ic[system.dim];

	gamma = 0.01;
	e = 0.0549; 					//Moon
	m_secondary = 1.215e-2; 		//Moon
	m_primary = 1.0 - m_secondary;
	G = 1.0;

	analysis.number_of_cycles = 1e3; //1e3 6e3
	analysis.cycle_period = 2.0 * M_PI;

	ic[0] = 0.1, ic[1] = 0.1;
	init_orbital(orbital, e);
	for (int i = 0; i < 4; i++) ic[i+2] = orbital[i];

	trace_orbit_map(ic, system, analysis);

	/////////////////////////////////////////////////////////
	/*				   		Phase space		   	           */
	/////////////////////////////////////////////////////////

	// system = system_rigid;

	// gamma = 0.01;
	// e = 0.0549; 					//Moon
	// m_secondary = 1.215e-2; 		//Moon
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;

	// analysis.nc = 5, analysis.nv = 15; //nc = 11, nv = 30;
	// analysis.number_of_cycles = 2e3; //1e3
	// analysis.cycle_period = 2.0 * M_PI;
	// analysis.coordinate_min = 0.0;
	// analysis.coordinate_max = 2.0 * M_PI; // M_PI 2.0* M_PI
	// analysis.velocity_min = 0.6;
	// analysis.velocity_max = 1.6;

	// draw_phase_space(system, analysis);

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