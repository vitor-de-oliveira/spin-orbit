#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
	double y[2], cycle_period;
	double **orbit;

    e = 0.1;
	y[0] = 0.1;
	y[1] = 0.05;
	cycle_period = 2.0 * M_PI;
 	t = 0.0;
	number_of_cycles = 1e3;
	
	// evolve_cycle(y, &e, cycle_period, &t);
	// printf("%1.10e %1.10e\n", y[0] , y[1]);

	FILE *tst_ev_orbit
			= fopen("output/test_evolve_orbit.dat", "w");

	evolve_orbit(&e, y, cycle_period, number_of_cycles, 
				&orbit);

	for (int i = 0; i < number_of_cycles; i++)
	{
		fprintf(tst_ev_orbit, "%1.15e %1.15e\n", 
			fmod(orbit[i][0], 2.0*M_PI), orbit[i][1]);
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

    // e = 0.99;
    // for (t = 0.0; t < 10.0 * M_PI; t += M_PI / 1000.0)
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
