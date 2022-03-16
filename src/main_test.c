#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "kepler.h"

int main(int argc, char **argv)
{
	/********************** Start clock ************************/

	clock_t begin = clock(), end;

	/******************* Declared variables ********************/

    double e;           // eccentricity
    double t;           // time

	/***********************************************************/
	
	////////////////////////////////////////////////////////////
	/*				    	Kepler equation	    			  */
	////////////////////////////////////////////////////////////

    FILE *tst_kep_plane = fopen("output/test_kepler_plane.dat", "w");
    FILE *tst_kep_time = fopen("output/test_kepler_time.dat", "w");

    // obs: -> cant cconverge on all points for e = 0.5 and e = 0.99
    //      -> negative time behave very similarly

    e = 0.99;
    for (t = 0.0; t > -10.0 * M_PI; t -= M_PI / 1000.0)
    {
        double u = kepler_equation(e,t);
        double r = 1 - e * cos(u);
        double f = 2.0 * atan(sqrt((1.0 + e)/(1.0 - e)) * tan(0.5 * u));
        double x = r * cos(f);
        double y = r * sin(f);
        fprintf(tst_kep_plane, "%1.10e %1.10e\n", x , y);
        fprintf(tst_kep_time, "%1.10e %1.10e\n", t , r);
    }

    fclose(tst_kep_plane);
    fclose(tst_kep_time);


	/********************** Stop clock *************************/

	end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	if (time_spent < 60.0)
	{
		printf("time spent = %1.2e seconds\n", time_spent);
	}
	else if (time_spent < 3600.0)
	{
		printf("time spent = %1.2e minutes\n", time_spent/60.0);
	}
	else
	{
		printf("time spent = %1.2e hours\n", time_spent/3600.0);
	}
	

	/***********************************************************/

	return 0;
}
