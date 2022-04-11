#include "auxiliar_functions.h"

int evolve_cycle(double *y, void *params, 
				double cycle_period, double *t, 
				char *system)
{
	// declare variables
	int status;
	double h, h_max, h_min;
	double error_abs, error_rel;
	char *method, *control;

	// initialiaze control variables
	h = 1e-3 * sign(cycle_period);
	error_abs = 1e-14;
	error_rel = 0.0;
	h_max = 1e-1;
	h_min = 1e-11;
	method = "rk8pd";
	control = "adaptive";

	// set system, integrator and driver
	gsl_odeiv2_system sys;
	set_system(&sys, params, system);

	const gsl_odeiv2_step_type *T;
	set_integrator(&T, method);

	gsl_odeiv2_driver *d;
	set_driver(&d, &sys, T, h, h_max, h_min,
			   error_abs, error_rel, control);

	// cycle evolution
	status = gsl_odeiv2_driver_apply (d, t, 
		*t + cycle_period, y);

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

int evolve_orbit(void *params, double *ic, 
				double cycle_period, int number_of_cycles, 
				double ***orbit, int *orbit_size,
				char *system)
{
	int system_dimension;

	if (strcmp(system, "rigid") == 0)
	{
		system_dimension = 6;
	}
	else if (strcmp(system, "rigid_kepler") == 0)
	{
		system_dimension = 2;
	}
	else if (strcmp(system, "two_body") == 0)
	{
		system_dimension = 4;
	}
	else
	{
		printf("Warning: undefined system\n");
		exit(2);
	}

	// declare variables
	double y[system_dimension];
	double box = 1e5;
	double t = 0.0;

	// allocate memory and initializes exit data
	if (orbit != NULL)
	{
		// takes into consideration initial condition
		alloc_2d_double(orbit, number_of_cycles + 1, 
			system_dimension);
		copy((*orbit)[0], ic, system_dimension);
	}
	else
	{
		printf("Error: NULL orbit.");
		exit(2);
	}
	
	// takes into consideration initial condition
	int counter = 1;

	// orbit evolution
	copy(y, ic, system_dimension);
	for (int i = 0; i < number_of_cycles; i++)
	{
		evolve_cycle(y, params, cycle_period, &t, system);
	
		// check if orbit diverges
		for (int j = 0; j < system_dimension; j++)
		{
			if (fabs(y[j]) > box)
			{
				printf("Warning: box limit reached\n");
				goto out;
			}
		}

		counter++;

		// write orbit element
		copy((*orbit)[i + 1], y, system_dimension);
	}

	out:;

	*orbit_size = counter;

	return 0;
}
