#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main_functions.h"
#include "auxiliar_functions.h"
#include "kepler.h"

int main(int argc, char **argv)
{
	/******************** Start clock **********************/

	clock_t begin = clock(), end;

	/***************** Declared variables ******************/

    double e;           // eccentricity
    double t;           // time

	/*******************************************************/

	/////////////////////////////////////////////////////////
	/*				   	Auxiliar functions    			  */
	/////////////////////////////////////////////////////////

	int number_of_cycles;
	int system_dimension;
	int orbit_size;
	double cycle_period;
	double **orbit;
	char *system;

    e = 0.0549; //Moon
	cycle_period = 2.0 * M_PI;
	// cycle_period = M_PI / 1000.0;
 	t = 0.0;
	number_of_cycles = 1e3;
	system_dimension = 6;
	system = "rigid";

	double perigee_over_semimajor_axis = 0.9451; //Moon
	double G = 1.0;
	double mass = 1.0;

	double y[system_dimension];
	y[0] = 0.1;
	y[1] = 0.1;
	y[2] = perigee_over_semimajor_axis;
	y[3] = 0.0;
	y[4] = 0.0;
	y[5] = sqrt(G * mass * (1.0 + e) / 
				perigee_over_semimajor_axis);

	// evolve_cycle(y, &e, cycle_period, &t);
	// printf("%1.10e %1.10e\n", y[0] , y[1]);

	FILE *tst_ev_orbit
			= fopen("output/test_evolve_orbit.dat", "w");

	evolve_orbit(&e, y, cycle_period, number_of_cycles, 
				&orbit, &orbit_size, system);

	for (int i = 0; i < orbit_size; i++)
	{
		// fprintf(tst_ev_orbit, "%1.15e %1.15e\n", 
		// 	fmod(orbit[i][0], 2.0*M_PI), orbit[i][1]);
		fprintf(tst_ev_orbit, "%1.15e %1.15e\n", 
			angle_mod(orbit[i][0]), orbit[i][1]);
	}

	dealloc_2d_double(&orbit, number_of_cycles);

    fclose(tst_ev_orbit);

	/////////////////////////////////////////////////////////
	/*				   	Kepler equation	    			  */
	/////////////////////////////////////////////////////////

    // FILE *tst_kep_plane 
	// 		= fopen("output/test_kepler_plane.dat", "w");
    // FILE *tst_kep_time 
	// 		= fopen("output/test_kepler_time.dat", "w");

    // e = 0.0549; //Moon
    // for (t = 0.0; t < 1.0 * M_PI; t += M_PI / 1000.0)
    // {
    //     double u = kepler_equation(e,t);
    //     double r = 1 - e * cos(u);
    //     double f = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) 
	// 			* tan(0.5 * u));
    //     double x = r * cos(f);
    //     double y = r * sin(f);
    //     fprintf(tst_kep_plane, "%1.10e %1.10e\n", x , y);
    //     fprintf(tst_kep_time, "%1.10e %1.10e\n", t , r);
    // }

    // fclose(tst_kep_plane);
    // fclose(tst_kep_time);

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
