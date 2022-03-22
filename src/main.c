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

    double e;           // eccentricity
    double t;           // time

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
	/*				   		Phase space		   	           */
	/////////////////////////////////////////////////////////

	e = 0.1;

	number_of_cycles = 1e3;
	cycle_period = 2.0 * M_PI;
    ic[0] = 0.1, ic[1] = 0.1;
	
	nc = 10, nv = 10;
	coordinate_min = 0.0;
	coordinate_max = 2.0 * M_PI;
	velocity_min = -0.1;
	velocity_max = 0.1;

	trace_orbit_map(ic, &e, cycle_period, number_of_cycles);

	// draw_phase_space(&e, cycle_period, 
	// 	number_of_cycles, coordinate_min, 
	// 	coordinate_max, velocity_min, 
	// 	velocity_max, nc, nv);

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