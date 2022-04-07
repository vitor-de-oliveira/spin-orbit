#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main_functions.h"

int main(int argc, char **argv)
{
	/******************** Start clock **********************/

	clock_t begin = clock(), end;

	/***************** Declared variables ******************/

    double gamma;			// equatorial flattening
    double e;				// eccentricity
	double m_primary;		// mass of primary
	double m_secondary;		// mass of secondary
	double G;				// gravitational constant

	double *params[5] 
		= {&gamma, &e, &m_primary, &m_secondary, &G};

	char *system;
	int nc, nv;
	int number_of_cycles;
	double ic[2];
	double cycle_period;
	double coordinate_min;
	double coordinate_max;
	double velocity_min;
	double velocity_max;

	/*******************************************************/
	
	/////////////////////////////////////////////////////////
	/*				   		   Orbit		   	           */
	/////////////////////////////////////////////////////////

	// gamma = 0.01;
	// e = 0.0549; 					//Moon
	// m_secondary = 1.215e-2; 		//Moon
	// m_primary = 1.0 - m_secondary;
	// G = 1.0;

	// number_of_cycles = 6e3; //1e3 6e3
	// cycle_period = 2.0 * M_PI;
	// system = "rigid_kepler";
    // ic[0] = 0.1, ic[1] = 0.1;

	// trace_orbit_map(ic, *params, cycle_period, number_of_cycles, 
	// 				system);

	/////////////////////////////////////////////////////////
	/*				   		Phase space		   	           */
	/////////////////////////////////////////////////////////

	gamma = 0.01;
	e = 0.0549; 					//Moon
	m_secondary = 1.215e-2; 		//Moon
	m_primary = 1.0 - m_secondary;
	G = 1.0;

	number_of_cycles = 2e3; //1e3
	cycle_period = 2.0 * M_PI;
	system = "rigid";

	nc = 11, nv = 30; //nc = 5, nv = 30;
	coordinate_min = 0.0;
	coordinate_max = 2.0 * M_PI; // M_PI 2.0* M_PI
	velocity_min = 0.6;
	velocity_max = 1.6;

	draw_phase_space(*params, cycle_period, 
		number_of_cycles, coordinate_min, 
		coordinate_max, velocity_min, 
		velocity_max, nc, nv, system);

	/******************** Stop clock ***********************/

	end = clock();
	double time_spent 
		= (double)(end - begin) / CLOCKS_PER_SEC;
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
	

	/*******************************************************/

	return 0;
}