#include "dynamical_system.h"

dynsys copy_dynsys(dynsys system)
{
	dynsys cp_system;
	cp_system.name = system.name;
	cp_system.field = system.field;
	cp_system.jac = system.jac;
	cp_system.dim = system.dim;
	cp_system.params = system.params;
	return cp_system;
}

anlsis copy_anlsis(anlsis analysis)
{
	anlsis cp_anlsis;
	cp_anlsis.nc = analysis.nc;
	cp_anlsis.nv = analysis.nv;
    cp_anlsis.grid_resolution = analysis.grid_resolution;
	cp_anlsis.number_of_cycles = analysis.number_of_cycles;
	cp_anlsis.cycle_period = analysis.cycle_period;
	cp_anlsis.coordinate_min = analysis.coordinate_min;
	cp_anlsis.coordinate_max = analysis.coordinate_max;
	cp_anlsis.velocity_min = analysis.velocity_min;
	cp_anlsis.velocity_max = analysis.velocity_max;
    cp_anlsis.grid_coordinate_min = analysis.grid_coordinate_min;
	cp_anlsis.grid_coordinate_max = analysis.grid_coordinate_max;
	cp_anlsis.grid_velocity_min = analysis.grid_velocity_min;
	cp_anlsis.grid_velocity_max = analysis.grid_velocity_max;
    cp_anlsis.evolve_box_size = analysis.evolve_box_size;
    cp_anlsis.evolve_basin_eps = analysis.evolve_basin_eps;
    cp_anlsis.time_series_delta = analysis.time_series_delta;
	return cp_anlsis;
}

int evolve_cycle(double *y, double *t,
				 dynsys system, anlsis analysis)
{
	// declare variables
	int status;
	double h, h_max, h_min;
	double error_abs, error_rel;
	char *method, *control;

	// initialize control variables
	h = 1e-3 * sign(analysis.cycle_period);
	error_abs = 1e-14;
	error_rel = 0.0;
	h_max = 1e-1;
	h_min = 1e-11;
	method = "rk8pd";
	control = "adaptive";

	// set system, integrator and driver
	gsl_odeiv2_system sys;
	set_system(&sys, system);

	const gsl_odeiv2_step_type *T;
	set_integrator(&T, method);

	gsl_odeiv2_driver *d;
	set_driver(&d, &sys, T, h, h_max, h_min,
			   error_abs, error_rel, control);

	// cycle evolution
	status = gsl_odeiv2_driver_apply (d, t, 
		*t + analysis.cycle_period, y);

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

int evolve_orbit(double *ic, double ***orbit, int *orbit_size,
				dynsys system, anlsis analysis)
{
	// declare variables
	double y[system.dim];
	double t = 0.0;

	// allocate memory and initializes exit data
	if (orbit != NULL)
	{
		// takes into consideration initial condition
		alloc_2d_double(orbit, analysis.number_of_cycles + 1, 
			system.dim);
		copy((*orbit)[0], ic, system.dim);
	}
	else
	{
		printf("Error: NULL orbit.");
		exit(2);
	}
	
	// takes into consideration initial condition
	int counter = 1;

	// orbit evolution
	copy(y, ic, system.dim);
	for (int i = 0; i < analysis.number_of_cycles; i++)
	{
		evolve_cycle(y, &t, system, analysis);
	
		// check if orbit diverges
		for (int j = 0; j < system.dim; j++)
		{
			if (fabs(y[j]) > analysis.evolve_box_size)
			{
				printf("Warning: box limit reached\n");
				printf("y[%d] = %1.10e\n", j, y[j]);
				printf("Box size = %1.2e\n", 
						analysis.evolve_box_size);
				goto out;
			}
		}

		counter++;

		// write orbit element
		copy((*orbit)[i + 1], y, system.dim);
	}

	out:;

	*orbit_size = counter;

	return 0;
}

int double_to_grid	(int grid[2], 
					double x[2],
					anlsis analysis)
{
	double res_i = (analysis.grid_coordinate_max -
					analysis.grid_coordinate_min) / 
					(double)(analysis.grid_resolution - 1);
	
	double res_j = (analysis.grid_velocity_max -
					analysis.grid_velocity_min) / 
					(double)(analysis.grid_resolution - 1);

	double pos_i = 
			(x[0] - analysis.grid_coordinate_min) / res_i;

	double pos_j = 
			(x[1] - analysis.grid_velocity_min) / res_j;

	grid[0] = (int)floor(pos_i);
	grid[1] = (int)floor(pos_j);

	// grid[0] = (int)pos_i;
	// grid[1] = (int)pos_j;

	return 0;
}

int grid_to_double	(int grid[2],
					double x[2],
					anlsis analysis)
{
	double res_i = (analysis.grid_coordinate_max -
					analysis.grid_coordinate_min) / 
					(double)(analysis.grid_resolution - 1);
	
	double res_j = (analysis.grid_velocity_max -
					analysis.grid_velocity_min) / 
					(double)(analysis.grid_resolution - 1);

	// x[0] = analysis.grid_coordinate_min + 
	// 		(double)grid[0] * res_i;
	// x[1] = analysis.grid_velocity_min + 
	// 		(double)grid[1] * res_j;

	x[0] = analysis.grid_coordinate_min + 
			((double)grid[0] + 0.5) * res_i;
	x[1] = analysis.grid_velocity_min + 
			((double)grid[1] + 0.5) * res_j;

}

int set_driver(gsl_odeiv2_driver **d, 
			   gsl_odeiv2_system *sys, 
			   const gsl_odeiv2_step_type *T, 
			   double h, double h_max, double h_min,
			   double error_abs, double error_rel, 
			   char *control)
{
	*d = gsl_odeiv2_driver_alloc_y_new(sys, 
		T, h, error_abs, error_rel);
		
	if (strcmp(control, "adaptive") == 0)
	{
		gsl_odeiv2_driver_set_hmax(*d, h_max);
		gsl_odeiv2_driver_set_hmin(*d, h_min);
	}
	else if (strcmp(control, "fixed") != 0)
	{
		printf("Warning: invalid GSL error control\n");
		exit(2);
	}
	return 0;
}

int set_integrator(const gsl_odeiv2_step_type **T, 
				   char *integrator)
{
    if (strcmp(integrator, "rk2") == 0)
	{
		*T = gsl_odeiv2_step_rk2;
	}
	else if (strcmp(integrator, "rk4") == 0)
	{
		*T = gsl_odeiv2_step_rk4;
	}
	else if (strcmp(integrator, "rk45") == 0)
	{
		*T = gsl_odeiv2_step_rkf45;
	}
	else if (strcmp(integrator, "rkck") == 0)
	{
		*T = gsl_odeiv2_step_rkck;
	}
	else if (strcmp(integrator, "rk8pd") == 0)
	{
		*T = gsl_odeiv2_step_rk8pd;
	}
	else if (strcmp(integrator, "rk4imp") == 0)
	{
		*T = gsl_odeiv2_step_rk4imp;
	}
	else if (strcmp(integrator, "bsimp") == 0)
	{
		*T = gsl_odeiv2_step_bsimp;
	}
	else if (strcmp(integrator, "msadams") == 0)
	{
		*T = gsl_odeiv2_step_msadams;
	}
	else if (strcmp(integrator, "msbdf") == 0)
	{
		*T = gsl_odeiv2_step_msbdf;	
	}
	else
	{
		printf("Warning: invalid GSL integrator\n");
		exit(2);
	}
	return 0;
}

int set_system(gsl_odeiv2_system *sys, 
			   dynsys system)
{
	sys->function = system.field;
	sys->jacobian = system.jac;
	sys->dimension = system.dim;
	sys->params = system.params;
	return 0;
}
