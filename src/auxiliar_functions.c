#include "auxiliar_functions.h"

int evolve_cycle(double *y, void *params, double cycle_period, double *t)
{
	// declare variables
	int status, system_dimension;
	double mu = *(double *)params, t0, h;

	// initialiaze control variables
	h = 1e-3 * sign(cycle_period);
	system_dimension = 2;
	t0 = *t;

	// set driver
	gsl_odeiv2_system sys = 
	{field_rigid, jacobian_rigid, system_dimension, &mu};
	gsl_odeiv2_driver *d = 
	gsl_odeiv2_driver_alloc_standard_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-14, 0.0, 0.0, 0.0);
	gsl_odeiv2_driver_set_hmax(d, 1e-3);
	gsl_odeiv2_driver_set_hmin(d, 1e-11);

	// cycle evolution
	status = gsl_odeiv2_driver_apply (d, t, t0 + cycle_period * sign(cycle_period), y);

	// check if integration was successfull
	if (status != GSL_SUCCESS)
	{
		printf("Warning: %s\n", gsl_strerror(status));
		exit(2);
	}

	// free memory
	gsl_odeiv2_driver_free(d);

	return 0;
}

int evolve_orbit(void *params, double *ic, double cycle_period, int number_of_cycles, double ***orbit)
{
	int system_dimension = 2;

	// declare variables
	double y[system_dimension];
	double box = 1e3;
	double t = 0.0;

	// allocate memory and initializes exit data
	if (orbit != NULL)
	{
		alloc_2d_double(orbit, number_of_cycles, system_dimension);
		copy((*orbit)[0], ic, system_dimension);
	}
	else
	{
		printf("Error: NULL orbit.");
		exit(2);
	}

	// orbit evolution
	copy(y, ic, system_dimension);
	for (int i = 1; i < number_of_cycles; i++)
	{
		evolve_cycle(y, params, cycle_period, &t);
	
		// check if orbit diverges
		if (fabs(y[0]) > box || fabs(y[1]) > box)
		{
			printf("Warning: box limit reached\n");
			break;
		}

		// write orbit element
		copy((*orbit)[i], y, system_dimension);
	}

	return 0;
}
