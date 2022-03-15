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

// int field_zero_velociy_curves_rotational(double t, const double y[], double f[], void *params)
// {
// 	(void)(t);
// 	double mu = *(double *)params;
// 	double mu_1 = 1. - mu;
// 	double mu_2 = mu;
// 	double r_1 = sqrt((y[0] + mu_2) * (y[0] + mu_2) + y[1] * y[1]);
// 	double r_2 = sqrt((y[0] - mu_1) * (y[0] - mu_1) + y[1] * y[1]);
// 	double omega_x = y[0] - mu_1 * (y[0] + mu_2) / (r_1 * r_1 * r_1) - mu_2 * (y[0] - mu_1) / (r_2 * r_2 * r_2);
// 	double omega_y = y[1] - mu_1 * y[1] / (r_1 * r_1 * r_1) - mu_2 * y[1] / (r_2 * r_2 * r_2);
// 	f[0] = 2. * omega_y;
// 	f[1] = -2. * omega_x;
// 	return GSL_SUCCESS;
// }

// double root_function_L3(double x, void *params)
// {
// 	struct root_params *p = (struct root_params *)params;
// 	double mu = p->mu;
// 	double J = p->J;
// 	return x * x + 2. * (1. - mu) / fabs(x + mu) + 2. * mu / fabs(x + mu - 1.) - J;
// }

// double root_derivative_L3(double x, void *params)
// {
// 	struct root_params *p = (struct root_params *)params;
// 	double mu = p->mu;
// 	return 2. * x - 2. * (1. - mu) * (x + mu) / (fabs(x + mu) * fabs(x + mu) * fabs(x + mu)) - 2. * mu * (x + mu - 1.) / (fabs(x + mu - 1.) * fabs(x + mu - 1.) * fabs(x + mu - 1.));
// }

// void root_fdf_L3(double x, void *params, double *y, double *dy)
// {
// 	struct root_params *p = (struct root_params *)params;
// 	double mu = p->mu;
// 	double J = p->J;
// 	*y = x * x + 2. * (1. - mu) / fabs(x + mu) + 2. * mu / fabs(x + mu - 1.) - J;
// 	*dy = 2. * x - 2. * (1. - mu) * (x + mu) / (fabs(x + mu) * fabs(x + mu) * fabs(x + mu)) - 2. * mu * (x + mu - 1.) / (fabs(x + mu - 1.) * fabs(x + mu - 1.) * fabs(x + mu - 1.));
// }

// double root_function_triangular_points(double y, void *params)
// {
// 	struct root_params *p = (struct root_params *)params;
// 	double mu = p->mu;
// 	double J = p->J;
// 	double pos[4];
// 	pos[0] = 0.5 - mu;
// 	pos[1] = y;
// 	pos[2] = 0.0;
// 	pos[3] = 0.0;
// 	return 2. * pseudo_potential(pos, &mu) - J; 
// }

// double root_derivative_triangular_points(double y, void *params)
// {
// 	struct root_params *p = (struct root_params *)params;
// 	double mu = p->mu;
// 	double pos[4];
// 	pos[0] = 0.5 - mu;
// 	pos[1] = y;
// 	pos[2] = 0.0;
// 	pos[3] = 0.0;
// 	return 2. * pseudo_potential_y_derivative(pos, &mu); 
// }

// void root_fdf_triangular_points(double y, void *params, double *z, double *dz)
// {
// 	struct root_params *p = (struct root_params *)params;
// 	double mu = p->mu;
// 	double J = p->J;
// 	double pos[4];
// 	pos[0] = 0.5 - mu;
// 	pos[1] = y;
// 	pos[2] = 0.0;
// 	pos[3] = 0.0;
// 	*z = 2. * pseudo_potential(pos, &mu) - J; 
// 	*dz = 2. * pseudo_potential_y_derivative(pos, &mu);
// }

// int root_zvc(double J, void *params, double *y, double x)
// {
// 	int iter, max_iter, status;
// 	double mu, x0;
// 	mu = *(double *)params;
// 	struct root_params params_root = {mu, J};
// 	const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_steffenson;
// 	gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
// 	gsl_function_fdf FDF;
// 	if (J >= lagrangian_point_jacobi_constant(&mu, 3))
// 	{
// 		FDF.f = &root_function_L3;
// 		FDF.df = &root_derivative_L3;
// 		FDF.fdf = &root_fdf_L3;
// 	}
// 	else if (J >= lagrangian_point_jacobi_constant(&mu, 4))
// 	{
// 		FDF.f = &root_function_triangular_points;
// 		FDF.df = &root_derivative_triangular_points;
// 		FDF.fdf = &root_fdf_triangular_points;
// 	}
// 	FDF.params = &params_root;
// 	max_iter = 100;
// 	iter = 0;
// 	gsl_root_fdfsolver_set(s, &FDF, x);
// 	do
// 	{
// 		iter++;
// 		gsl_root_fdfsolver_iterate(s);
// 		x0 = x;
// 		x = gsl_root_fdfsolver_root(s);
// 		status = gsl_root_test_delta(x, x0, 1e-15, 0);
// 		if (iter == max_iter)
// 		{
// 			printf("Warning: maximum iterate number reached at ZVC root\n");
// 			return 1;
// 		}
// 	} while (status == GSL_CONTINUE);
// 	if (J >= lagrangian_point_jacobi_constant(&mu, 3))
// 	{
// 		y[0] = x;
// 		y[1] = 0.;
// 	}
// 	else if (J >= lagrangian_point_jacobi_constant(&mu, 4))
// 	{
// 		y[0] = 0.5 - mu;
// 		y[1] = x;
// 	}
// 	gsl_root_fdfsolver_free(s);
// 	return 0;
// }
